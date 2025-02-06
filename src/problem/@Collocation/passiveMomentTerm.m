%======================================================================
%> @file passiveMomentTerm.m
%> @brief Collocation function to compute term minimizing passive torques
%> @details
%> Details: Collocation::passiveMomentTerm()
%>
%> @author Anne Koelewijn
%> @date May, 2020
%======================================================================

%======================================================================
%> @brief Collocation function to compute effort of the torques from states vector
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'states' of the model
%> @param exponent       (optional) Positive integer: Exponent of passive moments (default: 2)
%======================================================================
function output = passiveMomentTerm(obj,option,X,exponent)

%% check input parameter
% => Here in the metabolicRateTerm, it is not really worth initializing variables
if strcmp(option,'init')
    if ~isfield(obj.idx,'states') % check whether states are stored in X
        error('Model states are not stored in state vector X.')
    end

    % Return a dummy value
    output = NaN;
    return;
end

if nargin < 4
    exponent = 2;
elseif round(exponent) ~= exponent || exponent < 1
    error('Exponent must be a positive integer greater than 1.'); % For exponent = 1, the objective would not be differentiable
end

%% compute demanded output
if strcmp(option,'objval') % objective value
    output = 0;
    for i = 1:obj.nNodes
        x = X(obj.idx.states(:,i));
        M = obj.model.getPassiveJointmoments(x);
        output = output + sum(M.^exponent)/ obj.nNodes;
    
    end
elseif strcmp(option,'gradient') % gradient 
    output = zeros(size(X));
    for i = 1:obj.nNodes
        x = X(obj.idx.states(:,i));
        [M,dMdx] = obj.model.getPassiveJointmoments(x);
        output(obj.idx.states(1:obj.model.nDofs*2,i)) = exponent*dMdx*M.^(exponent-1)/ obj.nNodes;
    
    end  
else
    error('Unknown option');
end
end
