%======================================================================
%> @file trackAcc.m
%> @brief Collocation function to track accelerometer data
%> @details
%> Details: Collocation::trackAcc()
%>
%> @author Eva Dorschky
%> @date July, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to track accelerometer data
%>
%> @details
%> Supports data in millimeter per second squared ('mm/s^2') and
%> meter per second squared ('m/s^2').
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' and 'dur' of
%>                  the model
%> @param   data    TrackingData: All rows with type 'acc' will be tracked
%> @param   weight  (optional) String: Defines the weighting of the squared difference.
%>                  The following options exist:
%>                  - 'varOverTrials' (default) for nTrials > 1: \n 
%>                    \f$ \frac{\left(y_{i}[k]- \hat{y}_{i}[k] \right)^2}{\sigma^2_{y,i}[k]} \f$, 
%>                    where \f$ \sigma^2_{y,i}[k] \f$ is the variance over all trials of signal i for the time point k (given in the data table).
%>                  - 'varOverTrials' (default) for nTrials = 1: \n 
%>                    \f$ \frac{\left(y_{i}[k]- \hat{y}_{i}[k] \right)^2}{1} \f$ , 
%>                    where no weighting is performed as the variance over one trial would be zero.
%>                  - 'varOverTime': \n 
%>                    \f$ \frac{\left(y_{i}[k]- \hat{y}_{i}[k] \right)^2}{\sigma^2_{y,i}} \f$ ,
%>                    where \f$ \sigma^2_{y,i} \f$ is the variance over all time points of signal i.
%>                  - 'minMaxSquared': \n 
%>                    \f$ \frac{\left(y_{i}[k]- \hat{y}_{i}[k] \right)^2}{\left(max(y_i)-min(y_i) \right)^2 } \f$ ,
%>                    where the range of signal i over all time points is computed.
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient for input option 'gradient'
%======================================================================
function output = trackAcc(obj,option,X,data,weight)

fctname = 'trackAcc';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'dur')% check whether states and duration are stored in X
        error('Model states and duration need to be stored in state vector X.')
    end

    assert(ismethod(obj.model,'simuAccGyro'), ['Model ' class(obj.model) ' does not support accelerometer and gyroscope tracking.']);

    % initialize some variables (faster to get it once)
    variables = data.variables(strcmp(data.variables.type, 'acc'), :); % data with the correct type
    obj.objectiveInit.(fctname).variables = variables;
    obj.objectiveInit.(fctname).measMean = cell2mat(variables.mean'); % mean measured data
    if (nargin < 5 || strcmp(weight, 'varOverTrials')) && data.nTrials > 1
        % variance over trials of measured data
        obj.objectiveInit.(fctname).weighting = cell2mat(variables.var');
    elseif (nargin < 5 || strcmp(weight, 'varOverTrials')) && data.nTrials == 1
        % no weighting
        measMean = obj.objectiveInit.(fctname).measMean;
        obj.objectiveInit.(fctname).weighting = ones(size(measMean));
    elseif strcmp(weight, 'varOverTime')
        % variance over time
        measMean = obj.objectiveInit.(fctname).measMean;
        obj.objectiveInit.(fctname).weighting = repmat(var(measMean, 0, 1), size(measMean, 1), 1);
    elseif strcmp(weight, 'minMaxSquared')
        % (max-min)^2
        measMean = obj.objectiveInit.(fctname).measMean;
        obj.objectiveInit.(fctname).weighting = repmat((max(measMean, [], 1)-min(measMean, [], 1)).^2, size(measMean, 1), 1);
    else
        error('Weight method ''%s'' is unknown.', weight)
    end

    % get factor to convert simulated data to unit of tracking data
    % (assuming that all markers have the same unit)
    if strcmp(variables.unit{1}, 'mm/s^2')
        obj.objectiveInit.(fctname).factorToMeas = 1000;
    elseif strcmp(variables.unit{1}, 'm/s^2')
        obj.objectiveInit.(fctname).factorToMeas = 1;
    end

    % get index of segment, position and direction for each IMU
    % => It is faster to do this here instead of doing it in each call of simuAccGyro()
    [~, obj.objectiveInit.(fctname).idxSegment] = ismember(variables.segment, obj.model.segments.Properties.RowNames);
    obj.objectiveInit.(fctname).idxAcc = find(ismember(variables.type, 'acc'))';
    obj.objectiveInit.(fctname).dlocalAll = variables.direction;
    obj.objectiveInit.(fctname).plocalAll = variables.position;

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
% get variables from initalization (faster)
variables = obj.objectiveInit.(fctname).variables;
measMean = obj.objectiveInit.(fctname).measMean;
weighting = obj.objectiveInit.(fctname).weighting;
factorToMeas = obj.objectiveInit.(fctname).factorToMeas;
idxSegment = obj.objectiveInit.(fctname).idxSegment;
idxAcc = obj.objectiveInit.(fctname).idxAcc;
dlocalAll = obj.objectiveInit.(fctname).dlocalAll;
plocalAll = obj.objectiveInit.(fctname).plocalAll;

% initialize some more variables
duration = X(obj.idx.dur); % duration T of the simulation
x = X(obj.idx.states); % extract states
iq = obj.model.extractState('q');
iqd = obj.model.extractState('qdot');
nNodesDur = obj.nNodesDur;
h = duration/(nNodesDur-1);
nVars = height(variables); % total number of variables

if strcmp(option,'objval')

    % calculate objfun for accelerometer signals
    simVarAll = zeros(nNodesDur-1,nVars);
    
    for iNode = 1:nNodesDur-1
        % determine the q, qd, and qdd
        q = x(iq, iNode);
        qd = x(iqd, iNode);
        qdd = (x(iqd, iNode+1) - x(iqd, iNode) ) / h;
       
        % simulate accelerometer signals
        simVarAll(iNode, :) = obj.model.simuAccGyro(variables, q, qd, qdd, idxSegment, idxAcc, [], dlocalAll, plocalAll);
        
    end
    output = sum(sum((simVarAll*factorToMeas - measMean).^2./weighting)) / nVars / (nNodesDur-1);
    
elseif strcmp(option,'gradient')
    output = zeros(size(X));
    for iNode = 1:nNodesDur-1
        % determine the q, qd, and qdd
        q = x(iq, iNode);
        qd = x(iqd, iNode);
        qdd = (x(iqd, iNode+1) - x(iqd, iNode) ) / h;
        
        % calculate gradient for accelerometer signals
        [s, ds_dq, ds_dqd, ds_dqdd] = obj.model.simuAccGyro(variables, q, qd, qdd, idxSegment, idxAcc, [], dlocalAll, plocalAll);
        % consider computing ds_dT like this:
        % ds/dT = ds/dqdd * dqdd/dh * dh/dT (chain rule, note ds/dqdd is a Ns x Ndof matrix)
        %       = ds/dqdd * (qd(iNode+1) - qd(iNode)) * (-1/h^2) * 1/(nNodesDur-1)
        %       = -ds/dqdd * qdd / T 
        ds_dT = - ds_dqdd * qdd / duration;
        
        % get gradients
        df_dq   = 2*(s'* factorToMeas - measMean(iNode, :)) ./ weighting(iNode, :) * ds_dq*factorToMeas   /nVars /(nNodesDur-1);
        df_dqd  = 2*(s'* factorToMeas - measMean(iNode, :)) ./ weighting(iNode, :) * ds_dqd*factorToMeas  /nVars /(nNodesDur-1);
        df_dqdd = 2*(s'* factorToMeas - measMean(iNode, :)) ./ weighting(iNode, :) * ds_dqdd*factorToMeas /nVars /(nNodesDur-1);
        df_dT   = 2*(s'* factorToMeas - measMean(iNode, :)) ./ weighting(iNode, :) * ds_dT*factorToMeas   /nVars /(nNodesDur-1);

        % compile gradients
        output(obj.idx.states(iq,  iNode))   = df_dq';
        output(obj.idx.states(iqd, iNode))   = output(obj.idx.states(iqd, iNode))   + df_dqd' - df_dqdd'/h;
        output(obj.idx.states(iqd, iNode+1)) = output(obj.idx.states(iqd, iNode+1)) + df_dqdd'/h;
        output(obj.idx.dur)                  = output(obj.idx.dur) + df_dT;
    end

else
    error('Unknown option.')
end

