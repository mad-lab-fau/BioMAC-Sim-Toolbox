%======================================================================
%> @file equilibriumConstraintsCPOffset.m
%> @brief Collocation function to compute dynamic constraint for standing
%> with a vertical offset of the contact points
%> @details
%> Details: Collocation::equilibriumConstraintsCPOffset()
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================

%======================================================================
%> @brief Matlab function to compute dynamic constraint for standing
%> with a vertical offset of the contact points
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = equilibriumConstraintsCPOffset(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls')  || ~isfield(obj.idx,'CPYOffset')
    error('Model states and controls and the CPYOffset need to be stored in  vector X.')
end

%% compute demanded output
nconstraintspernode = obj.model.nConstraints;

if strcmp(option,'confun')
    output = zeros(nconstraintspernode*obj.nNodes,1);
    % dynamic equations must be zero
    for iNode=1:obj.nNodes
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        x = X(obj.idx.states(:,iNode));
        xd = zeros(size(x)); % state derivatives (zero in equilibrium)
        u = X(obj.idx.controls(:,iNode));

        output(ic) = obj.model.getDynamicsCPOffset(x,xd,u,X(obj.idx.CPYOffset));
    end
elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode*obj.nNodes,length(X),obj.Jnnz);
    for iNode = 1:obj.nNodes
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c

        x = X(obj.idx.states(:,iNode));
        xd = zeros(size(x)); % state derivatives (zero in equilibrium)
        ix = obj.idx.states(:,iNode);
        u = X(obj.idx.controls(:,iNode));
        iu = obj.idx.controls(:,iNode);

        [~, dfdx, ~, dfdu, dfdCPYOffset] = obj.model.getDynamicsCPOffset(x,xd,u,X(obj.idx.CPYOffset));
        output(ic,ix) = dfdx';
        output(ic,iu) = dfdu';
        output(ic, obj.idx.CPYOffset) = dfdCPYOffset;
    end

else
    error('Unknown option.');
end
end