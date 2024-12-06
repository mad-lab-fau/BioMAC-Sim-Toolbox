%======================================================================
%> @file @Collocation/dynamicConstraints.m
%> @brief Collocation function to compute dynamic constraint
%> @details
%> Details: Collocation::dynamicConstraints()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding dynamic equilibrium
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = dynamicConstraints(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('Model states and controls and duration need to be stored in state vector X.')
end

%% state variable indices for semi-implicit Euler method
% https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
% Discretization is defined by ixSIE1 (state variables from the node before the time step)
% and ixSIE2 (state variables from the node after the time step).
% Two versions: SIE A and SIE B. Both conserve energy.
% SIE A appeared easier to solve in tests.  This makes sense because it uses backward
% Euler for muscle equations. SIE B code is commented out and provided for comparison.
% SIE A: 
ixSIE2 = sort([obj.model.extractState('qdot'); obj.model.extractState('s'); obj.model.extractState('a')]);
ixSIE1 = setdiff(1:obj.model.nStates, ixSIE2);
% SIE B
% ixSIE1 = sort([obj.model.extractState('q')]);
% ixSIE2 = setdiff(1:obj.model.nStates, ixSIE1);


%% compute demanded output
nNodesDur = obj.nNodesDur;
h = X(obj.idx.dur)/(nNodesDur-1);
nconstraintspernode = obj.model.nConstraints;

if strcmp(option,'confun')
    output = zeros(nconstraintspernode*(nNodesDur-1),1);
    
    % dynamic equations must be zero
    for iNode=1:(nNodesDur-1)
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
        elseif strcmp(obj.Euler,'SIE')
            % semi-implicit Euler is energy neutral
			% take specified states ixSIE2 from x2, all other states from x1
			x1(ixSIE2) = x2(ixSIE2);
			output(ic) = obj.model.getDynamics(x1,xd,u2);
        end
    end
elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode*(nNodesDur-1),length(X),obj.Jnnz);
    
    for iNode = 1:(nNodesDur-1)
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        
        ix1 = obj.idx.states(:,iNode);
        ix2 = obj.idx.states(:,iNode+1);
        iu2 = obj.idx.controls(:,iNode+1);
        x1 = X(ix1);
        x2 = X(ix2);
        xd =(x2-x1)/h;
        u2 = X(obj.idx.controls(:,iNode+1));
        
        if strcmp(obj.Euler,'BE')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x2,xd,u2);
            output(ic,ix1) = -dfdxdot'/h;
            output(ic,ix2) = dfdx' + dfdxdot'/h;
        elseif strcmp(obj.Euler,'ME')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics((x1+x2)/2,xd,u2);
            output(ic,ix1) = dfdx'/2 - dfdxdot'/h;
            output(ic,ix2) = dfdx'/2 + dfdxdot'/h;
        elseif strcmp(obj.Euler,'SIE')
			x1(ixSIE2) = x2(ixSIE2);
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x1,xd,u2);
            output(ic,ix1) = -dfdxdot'/h;
            output(ic,ix2) = dfdxdot'/h;
            output(ic,ix1(ixSIE1)) = output(ic,ix1(ixSIE1)) + dfdx(ixSIE1,:)';
            output(ic,ix2(ixSIE2)) = output(ic,ix2(ixSIE2)) + dfdx(ixSIE2,:)';
		end
        output(ic,iu2) = dfdu';
        
        % derivative of constraints with respect to duration (because h is duration/(N-1))
        output(ic,obj.idx.dur) = -dfdxdot' * (x2-x1) / h^2 / (nNodesDur-1);
        
    end
else
    error('Unknown option.');
end
end
