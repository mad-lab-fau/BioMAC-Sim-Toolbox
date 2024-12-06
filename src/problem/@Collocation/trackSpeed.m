%======================================================================
%> @file trackSpeed.m
%> @brief Collocation function to track speed of the gait cycle
%> @details
%> Details: Collocation::trackSpeed()
%>
%> @author Anne Koelewijn
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Matlab function to track speed of the gait cycle
%> @details Computes difference between simulated and measured duration
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'speed' of the model
%> @param   data    TrackingData: All rows with type 'speed' will be tracked
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient for input option 'gradient'
%======================================================================
function output = trackSpeed(obj,option,X,data)

fctname = 'trackSpeed';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'speed') % check whether speed is stored in X
        error('Model speed is not stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    idxVar = ismember(data.variables.type, 'speed'); %index of tracking variable
    obj.objectiveInit.(fctname).measMean = data.variables.mean{idxVar};
    if data.nTrials == 1
        obj.objectiveInit.(fctname).measVar = 1; % use 1 to operate as if we are not dividing by the variance
    else
        obj.objectiveInit.(fctname).measVar = data.variables.var{idxVar}; % variance of measured data
    end

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
% get variables from initalization (faster)
measMean = obj.objectiveInit.(fctname).measMean;
measVar = obj.objectiveInit.(fctname).measVar;

if strcmp(option,'objval') %objective value for tracking speed
    simVar = X(obj.idx.speed); %state variable
    output = (simVar - measMean).^2./measVar;
    
elseif strcmp(option,'gradient') %gradient for tracking speed
    output = zeros(size(X));
    simVar = X(obj.idx.speed); %state variable
    output(obj.idx.speed) = 2*(simVar - measMean)/measVar;
else
    error('Unknown option');
end
end