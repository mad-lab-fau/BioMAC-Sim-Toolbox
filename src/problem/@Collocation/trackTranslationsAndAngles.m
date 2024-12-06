%======================================================================
%> @file trackTranslationsAndAngles.m
%> @brief Collocation function to track translations and angles in one term
%> @details
%> Details: Collocation::trackTranslationsAndAngles()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Computes difference between simulated and measured translations and angles in one term
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' of the model
%> @param   data    TrackingData: All rows with type 'angle' or 'translation' will be tracked
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient for input option 'gradient'
%======================================================================
function output = trackTranslationsAndAngles(obj,option,X,data)

fctname = 'trackTranslationsAndAngles';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states') % check whether model states are stored in X
        error('Model states are not stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    variables = data.variables(strcmp(data.variables.type, 'angle') | strcmp(data.variables.type, 'translation'), :); % data with the correct type
    obj.objectiveInit.(fctname).nVars = height(variables); % total number of variables
    obj.objectiveInit.(fctname).measMean = cell2mat(variables.mean'); % mean measured data
    if data.nTrials == 1
        obj.objectiveInit.(fctname).measVar = ones(size(obj.objectiveInit.(fctname).measMean)); % use 1 to operate as if we are not dividing by the variance
    else
        obj.objectiveInit.(fctname).measVar = cell2mat(variables.var'); % variance of measured data
    end

    obj.objectiveInit.(fctname).idxInX = obj.idx.states(obj.model.extractState('q', variables.name), 1:obj.nNodes); % indices in X of tracking variables

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
% get variables from initalization (faster)
nVars = obj.objectiveInit.(fctname).nVars;
measMean = obj.objectiveInit.(fctname).measMean;
measVar = obj.objectiveInit.(fctname).measVar;
idxInX = obj.objectiveInit.(fctname).idxInX;

if strcmp(option,'objval') %objective value for tracking angles and translations
    output = 0;
    for iVar = 1:nVars
        simVar = X(idxInX(iVar, :)); %simulated tracking variables

        output = output + ...
            sum((simVar - measMean(:, iVar)).^2./measVar(:, iVar)) / nVars /obj.nNodes;
    end
    
elseif strcmp(option,'gradient') %gradient for tracking angles and translations
    output = zeros(size(X));
    for iVar = 1:nVars
        simVar = X(idxInX(iVar, :)); %simulated tracking variables
        
        output(idxInX(iVar, :)) = output(idxInX(iVar, :)) +...
            2*(simVar - measMean(:, iVar))./measVar(:, iVar) / nVars / obj.nNodes;
    end
    
else
    error('Unknown option.');
end

end