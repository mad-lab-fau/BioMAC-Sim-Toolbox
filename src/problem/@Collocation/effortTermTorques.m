%======================================================================
%> @file effortTermTorques.m
%> @brief Collocation function to compute effort of the torques from states vector
%> @details
%> Details: Collocation::effortTermTorques()
%>
%> @author Marlies Nitschke
%> @date April, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to compute effort of the torques from states vector
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'controls' of the model
%> @param exponent       (optional) Positive integer: Exponent of torque controls (default: 2)
%> @param speedWeighting (optional) Boolean: If true, the objective will be divided by 
%>                       norm(speed)^exponent. (default: false)
%======================================================================
function output = effortTermTorques(obj,option,X,exponent,speedWeighting)

fctname = 'effortTermTorques';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'controls') % check whether controls are stored in X
        error('Model controls are not stored in state vector X.')
    end

    % initialize some variables (faster to get it once; even though they are not so expensive)
    obj.objectiveInit.(fctname).idxTorAllNodes = obj.idx.controls(obj.model.extractControl('torque'), 1:obj.nNodes);

    nTor = size(obj.objectiveInit.(fctname).idxTorAllNodes, 1);
    obj.objectiveInit.(fctname).weights = ones(nTor,1)/nTor;

    if nargin < 4
        exponent = 2;
    elseif round(exponent) ~= exponent || exponent < 1
        error('Exponent must be a positive integer greater than 1.'); % For exponent = 1, the objective would not be differentiable
    elseif mod(exponent, 2) == 1
        error('Exponent must be an even number as the torques can be negative.');
    end
    obj.objectiveInit.(fctname).exponent = exponent;

    if nargin < 5
        speedWeighting = 0;
    end
    obj.objectiveInit.(fctname).speedWeighting = speedWeighting;
    if speedWeighting && ~isfield(obj.idx,'speed') % check whether controls are stored in X
        error('Model speed is not stored in state vector X.')
    end

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
% get variables from initalization (faster)
idxTorAllNodes = obj.objectiveInit.(fctname).idxTorAllNodes;
weights = obj.objectiveInit.(fctname).weights;
exponent = obj.objectiveInit.(fctname).exponent;
speedWeighting = obj.objectiveInit.(fctname).speedWeighting;

% initialize some more variables
M = X(idxTorAllNodes);
if speedWeighting
    % Use the norm of the speed for weighting
    speed = X(obj.idx.speed);
    speedNorm = norm(speed);
end

if strcmp(option,'objval') % objective value
    output = sum(weights' * abs(M).^exponent)/ obj.nNodes;
    
    if speedWeighting
       output = output/speedNorm^exponent; 
    end
elseif strcmp(option,'gradient') % gradient 
    output = zeros(size(X));
    output(idxTorAllNodes) = exponent*repmat(weights,1,obj.nNodes) .* abs(M).^(exponent-1) .* sign(M)/ obj.nNodes;
    
    if speedWeighting
       output(idxTorAllNodes)  = output(idxTorAllNodes) /speedNorm^exponent;
       output(obj.idx.speed) = -exponent* sum(weights' * abs(M).^exponent)/ obj.nNodes/ speedNorm.^(exponent+2).*speed;
    end
    
else
    error('Unknown option');
end
end
