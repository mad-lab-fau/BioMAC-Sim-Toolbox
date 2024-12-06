%======================================================================
%> @file trackAngles.m
%> @brief Collocation function to track angles 
%> @details
%> Details: Collocation::trackAngles()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes difference between simulated and measured angles
%>
%> @details
%> Supports data in radian ('rad') and degree ('deg').
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' of the model
%> @param   data    TrackingData: All rows with type 'angle' will be tracked
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient
%======================================================================
function output = trackAngles(obj,option,X,data)

fctname = 'trackAngles';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states') % check whether model states are stored in X
        error('Model states are not stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    variables = data.variables(strcmp(data.variables.type, 'angle'), :); % data with the correct type
    obj.objectiveInit.(fctname).nVars = height(variables); % total number of variables
    obj.objectiveInit.(fctname).measMean = cell2mat(variables.mean'); % mean measured data
    if data.nTrials == 1
        obj.objectiveInit.(fctname).measVar = ones(size(obj.objectiveInit.(fctname).measMean)); % use 1 to operate as if we are not dividing by the variance
    else
        obj.objectiveInit.(fctname).measVar = cell2mat(variables.var'); % variance of measured data
    end

    obj.objectiveInit.(fctname).idxInX = obj.idx.states(obj.model.extractState('q', variables.name), 1:obj.nNodes); % indices in X of tracking variables

    % get factor to convert simulated data to unit of tracking data
    % (assuming that all angles have the same unit)
    if strcmp(variables.unit{1}, 'rad')
        obj.objectiveInit.(fctname).factorToMeas = 1;
    elseif strcmp(variables.unit{1}, 'deg')
        obj.objectiveInit.(fctname).factorToMeas = 180/pi;
    end

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
% get variables from initalization (faster)
nVars = obj.objectiveInit.(fctname).nVars;
measMean = obj.objectiveInit.(fctname).measMean;
measVar = obj.objectiveInit.(fctname).measVar;
factorToMeas = obj.objectiveInit.(fctname).factorToMeas;
idxInX = obj.objectiveInit.(fctname).idxInX;

if strcmp(option,'objval') %objective value for tracking angles  
    output = 0;
    for iVar = 1:nVars
        simVar = X(idxInX(iVar, :)); %simulated tracking variables

        output = output + ...
            sum((simVar*factorToMeas - measMean(:, iVar)).^2./measVar(:, iVar)) / nVars /obj.nNodes;
    end
    
elseif strcmp(option,'gradient') %gradient for tracking angles  
    output = zeros(size(X));
    for iVar = 1:nVars
        simVar = X(idxInX(iVar, :)); %simulated tracking variables

        output(idxInX(iVar, :)) = output(idxInX(iVar, :)) +...
            2*(simVar*factorToMeas - measMean(:, iVar))./measVar(:, iVar) * factorToMeas / nVars / obj.nNodes;
        
    end
    
else
    error('Unknown option.');
end

end