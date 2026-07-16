%======================================================================
%> @file dynamicPeriodicConstraints.m
%> @brief Collocation function to compute combined dynamics and periodicity constraint
%> @details
%> Details: Collocation::dynamicPeriodicConstraints()
%>
%> @author Antigravity
%> @date July, 2026
%======================================================================

%======================================================================
%> @brief Computes combined dynamics and periodicity constraint violation
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output: 'confun' or 'jacobian'
%> @param X             Double array: State vector containing states, controls, speed,
%>                      and duration of the movement (size N nodes, without N+1)
%> @param sym           Boolean if movement is symmetric (half period is optimized) or not
%======================================================================
function output = dynamicPeriodicConstraints(obj,option,X,sym)
%% check input parameters
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur') || ~isfield(obj.idx,'speed')
    error('Model states, controls, duration, and speed need to be stored in state vector X.')
end

if nargin < 4 || isempty(sym)
    sym = 1; % default to symmetric periodic
end

%% state variable indices for semi-implicit Euler method
if strcmp(obj.Euler,'SIE')
    ixSIE2 = sort([obj.model.extractState('qdot'); obj.model.extractState('s'); obj.model.extractState('a')]);
    ixSIE1 = setdiff(1:obj.model.nStates, ixSIE2);
end

%% compute demanded output
N = obj.nNodes; % Number of nodes (optimization nodes, N)
h = X(obj.idx.dur)/N;
nconstraintspernode = obj.model.nConstraints;
nStates = obj.model.nStates;
nControls = obj.model.nControls;

% Extract displacement parameters
dur = X(obj.idx.dur);
speed = X(obj.idx.speed);

unitdisplacementx = zeros(nStates,1);
unitdisplacementx(obj.model.idxForward) = 1;
unitdisplacementz = zeros(nStates,1);
unitdisplacementz(obj.model.idxSideward) = 1;

displacementx = unitdisplacementx*dur*speed(1);
if numel(speed) > 1
    displacementz = unitdisplacementz*dur*speed(2);
else
    displacementz = 0;
end

% Derivatives of displacement
displacementxddur   = unitdisplacementx*speed(1);
displacementxdspeed = unitdisplacementx*dur;
if numel(speed) > 1
    displacementzddur   = unitdisplacementz*speed(2);
    displacementzdspeed = unitdisplacementz*dur;
else
    displacementzddur   = 0;
    displacementzdspeed = [];
end
ddisp_ddur = displacementxddur + displacementzddur;

if strcmp(option,'confun')
    output = zeros(nconstraintspernode*N,1);
    
    for iNode=1:N
        ic = (1:nconstraintspernode) + (iNode-1)*nconstraintspernode;
        x1 = X(obj.idx.states(:,iNode));
        
        if iNode < N
            x2 = X(obj.idx.states(:,iNode+1));
            u2 = X(obj.idx.controls(:,iNode+1));
        else
            % Last interval goes from node N to virtual node N+1 (mirrored node 1)
            x1_first = X(obj.idx.states(:,1));
            u1_first = X(obj.idx.controls(:,1));
            if sym
                x2 = obj.model.idxSymmetry.xsign .* x1_first(obj.model.idxSymmetry.xindex) + displacementx + displacementz;
                u2 = obj.model.idxSymmetry.usign .* u1_first(obj.model.idxSymmetry.uindex);
            else
                x2 = x1_first + displacementx + displacementz;
                u2 = u1_first;
            end
        end
        
        xd = (x2-x1)/h;
        
        if strcmp(obj.Euler,'BE')
            output(ic) = obj.model.getDynamics(x2,xd,u2);
        elseif strcmp(obj.Euler,'ME')
            output(ic) = obj.model.getDynamics((x1+x2)/2,xd,u2);
        elseif strcmp(obj.Euler,'SIE')
            x1(ixSIE2) = x2(ixSIE2);
            output(ic) = obj.model.getDynamics(x1,xd,u2);
        end
    end
    
elseif strcmp(option,'jacobian')
    if isempty(obj.Jnnz)
        Jnnz = 1;
    else
        Jnnz = obj.Jnnz;
    end
    output = spalloc(nconstraintspernode*N,length(X),Jnnz);
    
    for iNode = 1:N
        ic = (1:nconstraintspernode) + (iNode-1)*nconstraintspernode;
        ix1 = obj.idx.states(:,iNode);
        x1 = X(ix1);
        
        if iNode < N
            ix2 = obj.idx.states(:,iNode+1);
            iu2 = obj.idx.controls(:,iNode+1);
            x2 = X(ix2);
            u2 = X(iu2);
            xd = (x2-x1)/h;
            
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
            
            % Derivative w.r.t duration for intervals 1 to N-1 (h = T/N, dh/dT = 1/N)
            output(ic,obj.idx.dur) = -dfdxdot' * (x2-x1) / (h^2 * N);
            
        else
            % Last interval N -> N+1 (virtual node)
            ix1_first = obj.idx.states(:,1);
            iu1_first = obj.idx.controls(:,1);
            x1_first = X(ix1_first);
            u1_first = X(iu1_first);
            
            if sym
                x2 = obj.model.idxSymmetry.xsign .* x1_first(obj.model.idxSymmetry.xindex) + displacementx + displacementz;
                u2 = obj.model.idxSymmetry.usign .* u1_first(obj.model.idxSymmetry.uindex);
            else
                x2 = x1_first + displacementx + displacementz;
                u2 = u1_first;
            end
            xd = (x2-x1)/h;
            
            if strcmp(obj.Euler,'BE')
                [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x2,xd,u2);
                dfdx_xN = -dfdxdot'/h;
                dfdx_x2 = dfdx' + dfdxdot'/h;
                dfdu_u2 = dfdu';
            elseif strcmp(obj.Euler,'ME')
                [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics((x1+x2)/2,xd,u2);
                dfdx_xN = dfdx'/2 - dfdxdot'/h;
                dfdx_x2 = dfdx'/2 + dfdxdot'/h;
                dfdu_u2 = dfdu';
            elseif strcmp(obj.Euler,'SIE')
                x1_sie = x1;
                x1_sie(ixSIE2) = x2(ixSIE2);
                [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x1_sie,xd,u2);
                dfdx_xN = -dfdxdot'/h;
                dfdx_xN(ixSIE1,:) = dfdx_xN(ixSIE1,:) + dfdx(ixSIE1,:)';
                dfdx_x2 = dfdxdot'/h;
                dfdx_x2(ixSIE2,:) = dfdx_x2(ixSIE2,:) + dfdx(ixSIE2,:)';
                dfdu_u2 = dfdu';
            end
            
            % Derivative w.r.t X(ix1), which is node N state variables
            output(ic,ix1) = dfdx_xN;
            
            % Derivative w.r.t first node states (x1_first) via chain rule on x2 = S(x1_first)
            if sym
                dfdx_x1_first = zeros(nconstraintspernode, nStates);
                dfdx_x1_first(:, obj.model.idxSymmetry.xindex) = dfdx_x2 .* (obj.model.idxSymmetry.xsign');
            else
                dfdx_x1_first = dfdx_x2;
            end
            % Accumulate into node 1 state variables
            output(ic,ix1_first) = output(ic,ix1_first) + dfdx_x1_first;
            
            % Derivative w.r.t first node controls (u1_first) via chain rule on u2 = S_u(u1_first)
            if sym
                dfdu_u1_first = zeros(nconstraintspernode, nControls);
                dfdu_u1_first(:, obj.model.idxSymmetry.uindex) = dfdu_u2 .* (obj.model.idxSymmetry.usign');
            else
                dfdu_u1_first = dfdu_u2;
            end
            % Accumulate into node 1 control variables
            output(ic,iu1_first) = output(ic,iu1_first) + dfdu_u1_first;
            
            % Derivative w.r.t duration
            % h = T/N, so:
            % dC/dT = dC/dxd * dxd/dT + dC/dx2 * dx2/dT
            % dxd/dT = - (x2 - x1) / (h^2 * N)
            % dx2/dT = ddisp_ddur
            dxd_dT = - (x2 - x1) / (h^2 * N);
            dC_dT = dfdxdot' * dxd_dT + dfdx_x2 * ddisp_ddur;
            output(ic,obj.idx.dur) = dC_dT;
            
            % Derivative w.r.t speed
            % dx2/dspeed = ddisp_dspeed
            ddisp_dspeed = [displacementxdspeed, displacementzdspeed];
            dC_dspeed = dfdx_x2 * ddisp_dspeed;
            output(ic,obj.idx.speed) = dC_dspeed;
        end
    end
else
    error('Unknown option.');
end
end
