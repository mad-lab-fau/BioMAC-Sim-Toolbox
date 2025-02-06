% ======================================================================
%> @file @Gait3d/simuAccGyro.m
%> @brief Gait3d function to simulates sensors placed on 3D model
%> @details
%> Details: Gait3d::simuAccGyro()
%>
%> @author Ton,  Marlies
%> @date July, 2017
% ======================================================================

%======================================================================
%> @brief Matlab function to simulates sensors placed on 3D model
%> @public
%>
%> @details
%> Do not change "obj" or "variables" in this function!  This will allow
%> Matlab to "pass by reference" and avoid function call overhead.
%>
%> Position has to given as X, Y and Z position of sensor in segment coordinate 
%> system in m.
%>
%> Direction has to be given in axis definition, e.g [1, 0, 0].
%>
%> @param obj               Gait3d class object
%> @param variables         Table: Variable table containing the accelerometer and gyroscope data
%>                          with at least the columns: type, name, segment, position, direction
%> @param q                 Generalized coordinates (Gait3d.nDofs x 1)
%> @param qd                First derivatives of generalized coordinates (Gait3d.nDofs x 1)
%> @param qdd               Second derivatives of generalized coordinates (Gait3d.nDofs x 1)
%>                          (use NaN to skip)
%> @param idxSegment        (optional) Double array: Indices of segments fitting to respective IMU (variables.nVars.all x 1)
%>                          => It is faster to pass them than to recalculate each time.
%> @param idxAcc            (optional) Double array: Indices in variables that have type 'acc' (nAcc x 1)
%>                          => It is faster to pass them than to recalculate each time.
%> @param idxGyro           (optional) Double array: Indices in variables that have type 'gyro' (nGyro x 1)
%>                          => It is faster to pass them than to recalculate each time.
%> @param dlocalAll         (optional) Double matrix: Local direction of the IMU saved in variables (variables.nVars.all x 3)
%>                          => It is faster to pass them than to extract them each time.
%> @param plocalAll         (optional) Double matrix: Local position of the IMU saved in variables (variables.nVars.all x 3)
%>                          => It is faster to pass them than to extract them each time.
%>
%> @retval s                Simulated sensor signals (variables.nVars.all x 1)
%> @retval ds_dq            Sparse Jacobian ds/dq (variables.nVars.all x Gait3d.nDofs)
%> @retval ds_dqd           Sparse Jacobian ds/dqd (variables.nVars.all x Gait3d.nDofs)
%> @retval ds_dqdd          Sparse Jacobian ds/dqdd (variables.nVars.all x Gait3d.nDofs)
%======================================================================
function [s, ds_dq, ds_dqd, ds_dqdd] = simuAccGyro(obj, variables, q, qd, qdd, idxSegment, idxAcc, idxGyro, dlocalAll, plocalAll)

%% General
% Get segment indices and check if all segments sepecified in the variables table exist
if nargin < 6
    [~, idxSegment] = ismember(variables.segment, obj.segments.Properties.RowNames);
end
if any(idxSegment==0)
    nonExistingSegments = unique(variables.segment(idxSegment==0));
    error('The following segments are not defined in the model: %s',  sprintf('%s ', nonExistingSegments{:}));
end

% Get indices of acc and gyro in variables tabale
if nargin < 7
    idxAcc = find(ismember(variables.type, 'acc'))';
end
if nargin < 8
    idxGyro = find(ismember(variables.type, 'gyro'))';
end
nAcc = numel(idxAcc); %number of accelerometer variables
nGyro = numel(idxGyro); %number of gyroscope variables
nVars = length(idxSegment); % much faster then size(variables, 1) or height(variables)

% Call the MEX function to help us compute the sensor signals
if (nAcc > 0 ) || (nGyro > 0)
    [FK, dFKdq, dFKdotdq] = obj.getFkin(q, qd);
    FKdot = dFKdq * qd;

    if nargin < 9
        dlocalAll = variables.direction; % get it for later
    end
end

% If we have accelerometers, we will also need FKdotdot and its
% Jacobians.  dFKdotdot/dq will be done with finite difference.
% Autolev can do it symbolically but this makes the C code 20%
% longer (slower) and we are already at the compiler limit.
if (nAcc > 0 )
    FKdotdot = dFKdotdq * qd + dFKdq * qdd;
    if (nargout > 1)
        h = 1e-7;
        q_new = q + h*qd;
        qd_new = qd + h*qdd;
        [~,~, dFKdotdq_new] = obj.getFkin(q_new, qd_new);
        dFKdotdotdq = (dFKdotdq_new - dFKdotdq) / h;
    end

    if nargin < 10
        plocalAll = variables.position; % get it for later
    end
end

% Initialize the outputs
s = zeros(nVars, 1);
if (nargout > 1)
    ds_dq = zeros(nVars, obj.nDofs);
    ds_dqd = zeros(nVars, obj.nDofs);
    ds_dqdd = zeros(nVars, obj.nDofs);
end

%% Accelerometer
% simulate accelerometer signal
if ~isempty(idxAcc)
    g = obj.gravity';     % use transpose because Opensim stores it as a row vector
    for iVar = idxAcc
        
        % Specific iVar
        plocal = plocalAll(iVar,:);
        if isrow(plocal); plocal = plocal'; end
        dlocal = dlocalAll(iVar,:);
        if isrow(dlocal); dlocal = dlocal'; end
        
        curSegment = idxSegment(iVar); %indices of tracking variables
        ip = (curSegment-2)*12 + (1:3);	% where the vector p is stored in FK
        iR = (curSegment-2)*12 + (4:12);	% where the matrix R is stored in FK
        pdd = FKdotdot(ip);
        R = reshape(FK(iR), 3, 3)';     % segment orientation matrix
        % 	Rd = reshape(FKdot(iR), 3, 3)';
        Rdd = reshape(FKdotdot(iR), 3, 3)';
        
        % Accelerometer model: s = (R*direction) . (rdd - g),
        % where rdd is the global acceleration of the sensor and g = (0,-9.81,0)
        % Rigid body model: r   = p   + R*plocal
        %                   rdd = pdd + Rdd*plocal
        % Therefore: a = (R*direction . (pdd + Rdd*plocal - g)
        R_dlocal = (R*dlocal)'; % don't calcuate it multiple times
        s(iVar) = R_dlocal * (pdd + Rdd*plocal - g);
        
        % Jacobians
        if (nargout > 1)
            % derivatives with respect to q, qd, and qdd
            for iDoF = find(any(dFKdq([ip,iR],:)))
                % ds/dq = (dR/dq*d)' * (pdd + Rdd*p - g) + (R * d)' * (dpdd/dq + dRdd/dq*p)
                ds_dq(iVar,iDoF)   = (reshape(dFKdq(iR,iDoF), 3, 3)'*dlocal)' * (pdd + Rdd*plocal - g) + R_dlocal * (dFKdotdotdq(ip,iDoF) + reshape(dFKdotdotdq(iR,iDoF), 3, 3)'*plocal);
                % ds/dqd = (dR/dq*d)' * (dpdd/dqd * Rdd/dqd*p) using dFKdd/dqd = dFKd/dq
                ds_dqd(iVar,iDoF)  = R_dlocal * (2*dFKdotdq(ip,iDoF) + reshape(2*dFKdotdq(iR,iDoF), 3, 3)'*plocal);
                % ds/dqdd = (dR/dq*d)' * (dpdd/dqdd * Rdd/dqdd*p) using dFKdd/dqdd = dFK/dq
                ds_dqdd(iVar,iDoF) = R_dlocal * (dFKdq(ip,iDoF) + reshape(dFKdq(iR,iDoF), 3, 3)'*plocal);
            end
        end
        
    end
end
%% Gyroscope
% simulate gyroscope signal
if ~isempty(idxGyro)
    for iVar = idxGyro
        
        % Specific iVar
        dlocal = dlocalAll(iVar,:);
        if isrow(dlocal); dlocal = dlocal'; end
        
        curSegment = idxSegment(iVar); %indices of tracking variables
        iR = (curSegment-2)*12 + (4:12);	% where the matrix R is stored in FK
        R = reshape(FK(iR), 3, 3)';
        Rd = reshape(FKdot(iR), 3, 3)';
        
        % Rate gyro model: s = direction . [W(3,2); W(1,3); W(2,1)] where W is the
        % angular velocity tensor in local coordinates: W = R' * dR/dt
        % e.g. https://shiyuzhao.wordpress.com/2011/06/08/rotation-matrix-angle-axis-angular-velocity/
        W = R' * Rd;
        s(iVar) = dlocal' * [W(3,2); W(1,3); W(2,1)];
        
        if (nargout > 1)
            % gyro signal only depends on q and qd
            for iDoF = find(any(dFKdq(iR,:)))
                dR_dq = reshape(dFKdq(iR,iDoF), 3, 3)';
                dRd_dq = reshape(dFKdotdq(iR,iDoF), 3, 3)';
                dW_dq = dR_dq' * Rd + R' * dRd_dq;
                dWd_qd = R' * dR_dq;
                % derivatives with respect to q, qd, and qdd
                ds_dq(iVar,iDoF)  = dlocal' * [dW_dq(3,2); dW_dq(1,3); dW_dq(2,1)];
                ds_dqd(iVar,iDoF) = dlocal' * [dWd_qd(3,2); dWd_qd(1,3); dWd_qd(2,1)];
                %ds_dqdd is zero
            end
        end
        
    end
end


end
