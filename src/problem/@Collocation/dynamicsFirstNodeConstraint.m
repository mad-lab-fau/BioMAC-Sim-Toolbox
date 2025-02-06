%======================================================================
%> @file dynamicsFirstNodeConstraint.m
%> @brief Collocation function to constrain the first node
%> @details
%> Details: Collocation::dynamicsFirstNodeConstraint()
%>
%> @author Marlies Nitschke
%> @date January, 2022
%======================================================================

%======================================================================
%> @brief Matlab function to constrain the first node
%>
%> @details
%> It ensures the dynamics at the first node:
%> f(x[1], (x[2]-x[1])/h, u[1]) = 0 without the constraints for qdot for
%> the pelvis
%>
%> This constraint can be used when no periodicity constraint is used.
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and 'duration' of movement
%======================================================================
function output = dynamicsFirstNodeConstraint(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls')   || ~isfield(obj.idx,'dur')
    error('Model states and controls and duration need to be stored in state vector X.')
end

%% compute demanded output
nNodesDur = obj.nNodesDur;
h = X(obj.idx.dur)/(nNodesDur-1);
nconstraintspernode = obj.model.nConstraints;

% Get indices
iNode = 1; % always do it for node 1
ix1 = obj.idx.states(:,iNode);
ix2 = obj.idx.states(:,iNode+1);
iu1 = obj.idx.controls(:,iNode);

% Get input for dynamics (forward Euler)
x1 = X(ix1);
x2 = X(ix2);
xd = (x2-x1)/h;
u1 = X(iu1);

if strcmp(option,'confun')

    output = obj.model.getDynamics(x1,xd,u1);
    output(1:6) = 0; % Do not constrain pelvis

elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode,length(X),obj.Jnnz);

    [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x1,xd,u1);
    output(:,ix1) = dfdx' -dfdxdot' /h;
    output(:,ix2) =  dfdxdot' /h;
    output(:,iu1) = dfdu';
    output(:,obj.idx.dur) = -dfdxdot' * (x2-x1) / h^2 / (nNodesDur-1);
    output(1:6, :) = 0; % Do not constrain pelvis

else
    error('Unknown option.');
end
end