%======================================================================
%> @file @Collocation/trackMarker.m
%> @brief Collocation function to track marker data 
%> @details
%> Details: Collocation::trackMarker()
%>
%> @author Marlies Nitschke
%> @date June, 2021
%======================================================================

%======================================================================
%> @brief Computes difference between simulated and measured marker data
%>
%> @details
%> Supports data in millimeter ('mm') and meter ('m').
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' of the model
%> @param   data    TrackingData: All rows with type 'marker' will be tracked
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient
%======================================================================
function output = trackMarker(obj,option,X,data)

fctname = 'trackMarker';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states') % check whether model states are stored in X
        error('Model states are not stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    variables = data.variables(strcmp(data.variables.type, 'marker'), :); % data with the correct type
    obj.objectiveInit.(fctname).variables = variables;
    obj.objectiveInit.(fctname).mMeasMean = cell2mat(variables.mean'); % mean measured marker data
    if data.nTrials == 1
        obj.objectiveInit.(fctname).mMeasVar = ones(size(obj.objectiveInit.(fctname).mMeasMean)); % use 1 to operate as if we are not dividing by the variance
    else
        obj.objectiveInit.(fctname).mMeasVar = cell2mat(variables.var'); % variance of measured data
    end

    % get factor to convert simulated data to unit of tracking data
    % (assuming that all markers have the same unit)
    if strcmp(variables.unit{1}, 'mm')
        obj.objectiveInit.(fctname).factorToMeas = 1000;
    elseif strcmp(variables.unit{1}, 'm')
        obj.objectiveInit.(fctname).factorToMeas = 1;
    end

    % get index of segment, position and direction for each marker
    % => It is faster to do this here instead of doing it in each call of simuMarker()
    [~, obj.objectiveInit.(fctname).idxSegment] = ismember(variables.segment, obj.model.segments.Properties.RowNames);
    obj.objectiveInit.(fctname).plocalAll = variables.position;
    obj.objectiveInit.(fctname).dlocalAll = variables.direction;

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
% get variables from initalization (faster)
variables = obj.objectiveInit.(fctname).variables;
mMeasMean = obj.objectiveInit.(fctname).mMeasMean;
mMeasVar = obj.objectiveInit.(fctname).mMeasVar;
factorToMeas = obj.objectiveInit.(fctname).factorToMeas;
idxSegment = obj.objectiveInit.(fctname).idxSegment;
plocalAll = obj.objectiveInit.(fctname).plocalAll;
dlocalAll = obj.objectiveInit.(fctname).dlocalAll;

% initialize some more variables
x = X(obj.idx.states); % extract states
idxq = obj.model.extractState('q'); % indices of x in state vector
nNodes = obj.nNodes; % number of nodes
nVars = height(variables); % number of tracking variables

if strcmp(option,'objval') % objective
    
    % simulate marker data
    mSim = zeros(nNodes,nVars);
    for iNode = 1:nNodes
        
        % get q
        q = x(idxq, iNode);
        
        % simulate data
        mSim(iNode, :) = obj.model.simuMarker(variables, q, idxSegment, plocalAll, dlocalAll);
        
    end
    output = sum(sum((mSim*factorToMeas - mMeasMean).^2./mMeasVar)) / nVars / nNodes;
    
elseif strcmp(option,'gradient') % gradient
    output = zeros(size(X));
    
    % get derivatives and compile gradient
    for iNode = 1:nNodes
        
        % get q
        q = x(idxq, iNode);
        
        % get marker positions and derivatives
        [mSim, dmSim_dq] = obj.model.simuMarker(variables, q, idxSegment, plocalAll, dlocalAll);
        dmSim_dq = dmSim_dq * factorToMeas;
        
        % compile gradient
        output(obj.idx.states(idxq,iNode)) = 2*( mSim'* factorToMeas - mMeasMean(iNode, :)) ./ mMeasVar(iNode, :) * dmSim_dq  /nVars /nNodes;
    end

else
    error('Unknown option.');
end


end