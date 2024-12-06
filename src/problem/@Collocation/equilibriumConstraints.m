%======================================================================
%> @file equilibriumConstraints.m
%> @brief Collocation function to compute dynamic constraint for standing
%> @details
%> Details: Collocation::equilibriumConstraints()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Matlab function to compute dynamic constraint for standing
%>
%> @details
%> It ensures equilibrium at node 1.
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = equilibriumConstraints(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls')  % check whether controls are stored in X
    error('Model states and controls need to be stored in state vector X.')
end

%% compute demanded output
nconstraintspernode = obj.model.nConstraints;

% Get indices
iNode = 1; % always do it for node 1
ix = obj.idx.states(:,iNode);
iu = obj.idx.controls(:,iNode);

% Get input for dynamics
x = X(ix);
xd = zeros(size(x)); % state derivatives (zero in equilibrium)
u = X(iu);

if strcmp(option,'confun')

    % dynamic equations must be zero
    output = obj.model.getDynamics(x,xd,u);

elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode,length(X),obj.Jnnz);

    [~, dfdx, ~, dfdu] = obj.model.getDynamics(x,xd,u);
    output(:,ix) = dfdx';
    output(:,iu) = dfdu';
    
else
    error('Unknown option.');
end
end