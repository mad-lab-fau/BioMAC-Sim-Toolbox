%======================================================================
%> @file trackGyro.m
%> @brief Collocation function to track gyroscope data
%> @details
%> Details: Collocation::trackGyro()
%>
%> @author Eva Dorschky
%> @date July, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to track gyroscope data
%>
%> @details
%> Supports data in rad per second ('rad/s') and degree per second ('deg/s').
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' of the model
%> @param   data    TrackingData: All rows with type 'gyro' will be tracked
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
function output = trackGyro(obj,option,X,data,weight)

fctname = 'trackGyro';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states') % check whether states are stored in X
        error('Model states need to be stored in state vector X.')
    end

    assert(ismethod(obj.model,'simuAccGyro'), ['Model ' class(obj.model) ' does not support accelerometer and gyroscope tracking.']);

    % initialize some variables (faster to get it once)
    variables = data.variables(strcmp(data.variables.type, 'gyro'), :); % data with the correct type
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
    if strcmp(variables.unit{1}, 'rad/s')
        obj.objectiveInit.(fctname).factorToMeas = 1;
    elseif strcmp(variables.unit{1}, 'deg/s')
        obj.objectiveInit.(fctname).factorToMeas = 180/pi;
    end

    % get index of segment and direction for each IMU
    % => It is faster to do this here instead of doing it in each call of simuAccGyro()
    [~, obj.objectiveInit.(fctname).idxSegment] = ismember(variables.segment, obj.model.segments.Properties.RowNames);
    obj.objectiveInit.(fctname).idxGyro = find(ismember(variables.type, 'gyro'))';
    obj.objectiveInit.(fctname).dlocalAll =  variables.direction;

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
idxGyro = obj.objectiveInit.(fctname).idxGyro;
dlocalAll = obj.objectiveInit.(fctname).dlocalAll;

% initialize some more variables
x = X(obj.idx.states); % extract states
iq = obj.model.extractState('q');
iqd = obj.model.extractState('qdot');
nVars = height(variables); % total number of variables

if strcmp(option,'objval')

    % calculate objfun for gyroscope signals
    simVarAll = zeros(obj.nNodes,nVars);
    
    for iNode = 1:obj.nNodes
        % determine the q and qd
        q = x(iq, iNode);
        qd = x(iqd, iNode);
        
        % simulate gyroscope signals
        simVarAll(iNode, :) = obj.model.simuAccGyro(variables, q, qd, NaN, idxSegment, [], idxGyro, dlocalAll);
        
    end
    output = sum(sum((simVarAll*factorToMeas - measMean).^2./weighting)) / nVars / obj.nNodes;
    
elseif strcmp(option,'gradient')
    output = zeros(size(X));
    for iNode = 1:obj.nNodes
        % determine the q and qd
        q = x(iq, iNode);
        qd = x(iqd, iNode);
        
        % calculate gradient for gyroscope signals
        [s, ds_dq, ds_dqd] = obj.model.simuAccGyro(variables, q, qd, NaN, idxSegment, [], idxGyro, dlocalAll);
        
        % get and compile gradient
        output(obj.idx.states(iq,  iNode)) = 2*(s'* factorToMeas - measMean(iNode, :)) ./ weighting(iNode, :) * ds_dq*factorToMeas  /nVars /obj.nNodes;
        output(obj.idx.states(iqd, iNode)) = 2*(s'* factorToMeas - measMean(iNode, :)) ./ weighting(iNode, :) * ds_dqd*factorToMeas /nVars /obj.nNodes;

    end

else
    error('Unknown option.')
end

