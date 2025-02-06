%======================================================================
%> @file effortTermMuscles.m
%> @brief Collocation function to compute effort of the muscles from states vector
%> @details
%> Details: Collocation::effortTermMuscles()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Matlab function to compute effort of the muscles from states vector
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'controls' of the model
%> @param weigthsType    (optional) String parsing weights
%> @param exponent       (optional) Positive integer: Exponent of muscle controls (default: 3)
%> @param speedWeighting (optional) Boolean: If true, the objective will be divided by 
%>                       norm(speed)^exponent. (default: false)
%======================================================================
function output = effortTermMuscles(obj,option,X,weigthsType,exponent,speedWeighting)

fctname = 'effortTermMuscles';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'controls') % check whether controls are stored in X
        error('Model controls are not stored in state vector X.')
    end

    % initialize some variables (faster to get it once; even though they are not so expensive)
    obj.objectiveInit.(fctname).idxNeuralExAllNodes = obj.idx.controls(obj.model.extractControl('u'), 1:obj.nNodes);

    nNeuralEx = size(obj.objectiveInit.(fctname).idxNeuralExAllNodes, 1);
    if nargin < 4
        weights = ones(nNeuralEx,1)/nNeuralEx;  % this is the default
    elseif strcmp(weigthsType,'equal')
        weights = ones(nNeuralEx,1)/nNeuralEx;
    elseif strcmp(weigthsType,'volumeweighted')
        if sum(ismember(obj.model.muscles.Properties.VariableNames, 'weight')) == 0 %weight is not stored for Gait2dc model
            weights = obj.model.muscles.fmax.*obj.model.muscles.lceopt; %musMass = (fmax/sigma)*rho.*l_ceopt;
            weights = weights/sum(weights);
        else
            weights = obj.model.muscles.weight;
            weights = weights/sum(weights);
        end
    else
        error('Unknown type of weights.');
    end
    obj.objectiveInit.(fctname).weights = weights;

    if nargin < 5
        obj.objectiveInit.(fctname).exponent = 3;
    elseif round(exponent) ~= exponent || exponent < 1
        error('Exponent must be a positive integer');
    end
    obj.objectiveInit.(fctname).exponent = exponent;

    if nargin < 6
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
idxNeuralExAllNodes = obj.objectiveInit.(fctname).idxNeuralExAllNodes;
weights = obj.objectiveInit.(fctname).weights;
exponent = obj.objectiveInit.(fctname).exponent;
speedWeighting = obj.objectiveInit.(fctname).speedWeighting;

% initialize some more variables
if speedWeighting
    % Use the norm of the speed for weighting
    speed = X(obj.idx.speed);
    speedNorm = norm(speed); 
end

if strcmp(option,'objval') % objective value
    output = sum(weights'*X(idxNeuralExAllNodes).^exponent)/ obj.nNodes;
    
    if speedWeighting
        output = output / speedNorm.^exponent;
    end
elseif strcmp(option,'gradient') % gradient 
    output = zeros(size(X));
    output(idxNeuralExAllNodes) = exponent*repmat(weights,1,obj.nNodes).*X(idxNeuralExAllNodes).^(exponent-1)/ obj.nNodes;
    
    if speedWeighting
        output(idxNeuralExAllNodes) = output(idxNeuralExAllNodes) / speedNorm.^exponent;
        output(obj.idx.speed) = -exponent* sum(weights'*X(idxNeuralExAllNodes).^exponent)/ obj.nNodes/ speedNorm.^(exponent+2).*speed;
    end
else
    error('Unknown option');
end
end