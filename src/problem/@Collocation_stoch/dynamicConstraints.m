%======================================================================
%> @file @Collocation_stoch/dynamicConstraints.m
%> @brief Collocation_stoch function to compute dynamic constraint
%> @details
%> Details: Collocation_stoch::dynamicConstraints()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding dynamic equilibrium
%>
%> @param obj           Collocation_Exo class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = dynamicConstraints(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('Model states and controls and duration need to be stored in state vector X.')
end
assert(ismember(obj.Euler, {'BE', 'ME'}), 'Only BE and ME discretization is implemented.')

%% compute demanded output
h = X(obj.idx.dur)/obj.nNodes;
nconstraintspernode = obj.model.nConstraints;

if strcmp(option,'confun')
    output = zeros(nconstraintspernode*obj.nNodes,1);
    
    % dynamic equations must be zero
    for iNode=1:obj.nNodes
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        x1 = X(obj.idx.states(:,iNode));
        x2 = X(obj.idx.states(:,iNode+1));
        xd =(x2-x1)/h;
        u2 = X(obj.idx.controls(:,iNode+1));
        
        if strcmp(obj.Euler,'BE')
            output(ic) = obj.model.getDynamics(x2,xd,u2);	% backward Euler discretization
        elseif strcmp(obj.Euler,'ME')
            % we're using u2 instead of (u1+u2)/2 because u is
            % open loop and it makes no difference except u2
            % will converge more easily
            output(ic) = obj.model.getDynamics((x1+x2)/2,xd,u2);
        end
        output(ic) = output(ic) + obj.randvals(:,iNode);
    end
elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode*obj.nNodes,length(X),obj.Jnnz);
    
    for iNode = 1:obj.nNodes
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        
        x1 = X(obj.idx.states(:,iNode));
        ix = obj.idx.states(:,iNode);
        x2 = X(obj.idx.states(:,iNode+1));
        ix2 = obj.idx.states(:,iNode+1);
        xd =(x2-x1)/h;
        u2 = X(obj.idx.controls(:,iNode+1));
        iu2 = obj.idx.controls(:,iNode+1);
        
        if strcmp(obj.Euler,'BE')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x2,xd,u2);
            output(ic,ix) = -dfdxdot'/h;
            output(ic,ix2) = dfdx' + dfdxdot'/h;
        elseif strcmp(obj.Euler,'ME')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics((x1+x2)/2,xd,u2);
            output(ic,ix) = dfdx'/2 - dfdxdot'/h;
            output(ic,ix2) = dfdx'/2 + dfdxdot'/h;
        end
        output(ic,iu2) = dfdu';
        
        % derivative of constraints with respect to duration (because h is duration/N)
        output(ic,obj.idx.dur) = -dfdxdot' * (x2-x1) / h^2 / obj.nNodes;
        
    end
else
    error('Unknown option.');
end
end