%======================================================================
%> @file trackDuration.m
%> @brief Collocation function to track duration of the gait cycle
%> @details
%> Details: Collocation::trackDuration()
%>
%> @author Anne Koelewijn
%> @date November, 2017
%======================================================================


%======================================================================
%> @brief Matlab function to track duration of the gait cycle
%> @details 
%> Computes difference between simulated and measured duration
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'duration' of the model
%> @param   data    TrackingData: The row with type 'duration' will be tracked.
%>                  It assumes that there is only one row.
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient for input option 'gradient'
%======================================================================

function output = trackDuration(obj,option,X,data)
fctname = 'trackDuration';

% For tracking just the numerical duration
if isscalar(data)
    if strcmp(option,'objval') %objective value for tracking duration  
        simVar = X(obj.idx.dur); %state variable
        output = (simVar - data).^2;
    elseif strcmp(option,'gradient') %gradient for tracking duration
        output = zeros(size(X));

        simVar = X(obj.idx.dur); %state variable
        output(obj.idx.dur) = 2*(simVar - data);
    elseif strcmp(option,'init')
    else
        error('Unknown option');
    end
    return;
end

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'dur') % check whether duration is stored in X
        error('Duration is not stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    idxVar = strcmp(data.variables.type, 'duration'); % assume that it is only one row
    obj.objectiveInit.(fctname).measMean = data.variables.mean{idxVar}; % mean measured data
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

if strcmp(option,'objval') %objective value for tracking duration  
    simVar = X(obj.idx.dur); %state variable
    output = (simVar - measMean).^2./measVar;
    
elseif strcmp(option,'gradient') %gradient for tracking duration 
    output = zeros(size(X));
    
    simVar = X(obj.idx.dur); %state variable
    output(obj.idx.dur) = 2*(simVar - measMean)/measVar;
else
    error('Unknown option');
end
end