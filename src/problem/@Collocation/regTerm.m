%======================================================================
%> @file regTerm.m
%> @brief Collocation function to compute regularization term
%> @details
%> Details: Collocation::regTerm()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes differences between collocation nodes to ensure
%> smoothness
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output: 'objval' or 'gradient'
%>                      (or 'init' for initialization)
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model 
%======================================================================
function output = regTerm(obj,option,X)

%% check input parameter
% => Here in the regTerm, it is not really worth initializing variables
if strcmp(option,'init')
    if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur')% check whether controls are stored in X
        error('Model states and controls and dur need to be stored in state vector X.')
    end

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
xd.states = diff(X(obj.idx.states),1,2);
xd.controls = diff(X(obj.idx.controls),1,2);
nvarpernode = size(obj.idx.states,1)+size(obj.idx.controls,1);
duration = X(obj.idx.dur);
nNodesDur = obj.nNodesDur;
factor = (nNodesDur-1)/nvarpernode/duration^2; % see 1/(N-1)/nVars * sum(sum( ((x2-x1)/(dur/(N-1)))^2 )) = (N-1)/nVars/dur^2 * sum(sum( (x2-x1)^2 ))

if strcmp(option,'objval') 
    output = factor*(sum(sum(xd.states.^2)) + sum(sum(xd.controls.^2))); % make it the average
elseif strcmp(option,'gradient')
    output = zeros(size(X));
    output(obj.idx.states(:,1:end-1))   = output(obj.idx.states(:,1:end-1))   - 2*factor*xd.states;
    output(obj.idx.states(:,2:end))     = output(obj.idx.states(:,2:end))     + 2*factor*xd.states;
    output(obj.idx.controls(:,1:end-1)) = output(obj.idx.controls(:,1:end-1)) - 2*factor*xd.controls;
    output(obj.idx.controls(:,2:end))   = output(obj.idx.controls(:,2:end))   + 2*factor*xd.controls;
    output(obj.idx.dur)                 = -2*(nNodesDur-1)/nvarpernode/duration^3*(sum(sum(xd.states.^2)) + sum(sum(xd.controls.^2)));
else
    error('Unknown option.');
end

end