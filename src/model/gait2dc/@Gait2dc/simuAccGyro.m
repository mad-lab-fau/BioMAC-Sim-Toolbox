%======================================================================
%> @file @Gait2dc/simuAccGyro.m
%> @brief Gait2dc function to simulate acceleration and gyroscope signals
%> @details
%> Details: Gait2dc::simuAccGyro()
%>
%> @author Ton, Eva, Iris, Ann-Kristin, Marlies
%> @date July, 2017
%======================================================================

%======================================================================
%> @brief Function to simulate acceleration and gyroscope signals
%> @public
%> @details
%> Do not change "obj" or "variables" in this function!  This will allow
%> Matlab to "pass by reference" and avoid function call overhead.
%>
%> Position has to given as X and Y position of sensor in segment coordinate 
%> system in m.
%>
%> Direction has to be given in axis definition, e.g [1, 0, 0].
%>
%> @param obj               Gait2dc class object
%> @param variables         Table: Variable table containing the accelerometer and gyroscope data
%>                          with at least the columns: type, name, segment, position, direction
%> @param q                 Generalized coordinates (Gait2dc.nDofs x 1)
%> @param qd                First derivatives of generalized coordinates (Gait2dc.nDofs x 1)
%> @param qdd               Second derivatives of generalized coordinates (Gait2dc.nDofs x 1)
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
%> @retval s                Simulated sensor signals (data.nVars.all x 1)
%> @retval ds_dq            Sparse Jacobian ds/dq (data.nVars.all x Gait2dc.nDofs)
%> @retval ds_dqd           Sparse Jacobian ds/dqd (data.nVars.all x Gait2dc.nDofs)
%> @retval ds_dqdd          Sparse Jacobian ds/dqdd (data.nVars.all x Gait2dc.nDofs)
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


% Initialize the outputs
s = zeros(nVars, 1);
if (nargout > 1)
    ds_dq = zeros(nVars, obj.nDofs);
    ds_dqd = zeros(nVars, obj.nDofs);
    ds_dqdd = zeros(nVars, obj.nDofs);
end

%% Accelerometer
% determine further values for Accelerometer tracking
if nAcc > 0
    hh = 1e-7;
    % determine the 42 accelerometer model coefficients (7 segments x 6 model terms), and their Jacobians
    if (0)
        if nargout > 1
            [a, da_dq, da_dqd, da_dqdd] = gait2dc('Accelerometer',q,qd,qdd);
        else
            a = gait2dc('Accelerometer',q,qd,qdd);
        end
    else
        a = gait2dc('Accelerometer',q,qd,qdd);
        
        if nargout > 1
            % The three Jacobians of a are calculated here by finite differences, because we have a crash if we
            % get them from gait2dc.mexw32.  This needs to be looked
            % at! --> No crash at the moment
            da_dq = zeros(42,obj.nDofs);
            da_dqd = zeros(42,obj.nDofs);
            da_dqdd = zeros(42,obj.nDofs);
            
            for k = 1:obj.nDofs
                tt = q(k);
                q(k) = q(k) + hh;
                da_dq(:,k) = (gait2dc('Accelerometer',q,qd,qdd) - a) / hh;
                q(k) = tt;
                
                tt = qd(k);
                qd(k) = qd(k) + hh;
                da_dqd(:,k) = (gait2dc('Accelerometer',q,qd,qdd) - a) / hh;
                qd(k) = tt;
                
                tt = qdd(k);
                qdd(k) = qdd(k) + hh;
                da_dqdd(:,k) = (gait2dc('Accelerometer',q,qd,qdd) - a) / hh;
                qdd(k) = tt;
            end
            da_dq = sparse(da_dq);
            da_dqd = sparse(da_dqd);
            da_dqdd = sparse(da_dqdd);
        end
    end
    
end


% simulate accelerometer signal
if ~isempty(idxAcc)
    for iVar = idxAcc
        
        % Specific iVar
        curSegment = idxSegment(iVar); %indices of tracking variables
        curPosition = plocalAll(iVar,:);
        curDirection = dlocalAll(iVar,:);
        
        c = [ [1 curPosition(1) curPosition(2)]*curDirection(1), [1 curPosition(1) curPosition(2)]*curDirection(2) ];
        ia = 6*(curSegment-1) + (1:6);
        
        % accelerometer signal s
        value = a(ia);
        s(iVar) = c*value;
        % Jacobians ds/dq, ds/dqd and ds/dh
        if nargout > 1
            ds_dq(iVar, :)  = c*da_dq(ia,:);      % because ds_dq = c*da_dq(ia,:);
            ds_dqd(iVar, :) = c*da_dqd(ia,:);     % because ds_dqd = c*da_dqd(ia,:);
            ds_dqdd(iVar, :)= c*da_dqdd(ia,:);   % because ds_dqdd = c*da_dqdd(ia,:);
        end
    end
end
%% Gyroscope

if nGyro > 0
    % indices for qd for the specific segments (see Iris Kellermanns thesis)
    iSeg_qd{1} = [3];       % trunk
    iSeg_qd{2} = [3:4];     % trunk, right thigh
    iSeg_qd{3} = [3:5];     % trunk, right thigh, right shank
    iSeg_qd{4} = [3:6];     % trunk, right thigh, right shank, right ankle
    iSeg_qd{5} = [3,7];     % trunk, left thigh
    iSeg_qd{6} = [3,7:8];   % trunk, left thigh , left shank
    iSeg_qd{7} = [3,7:9];   % trunk, left thigh , left shank ,left ankle
end


% simulate gyroscope signal
if ~isempty(idxGyro)
    for iVar = idxGyro
        
        % Specific iVar
        curSegment = idxSegment(iVar); %indices of tracking variables
        curDirection = dlocalAll(iVar,:);
        
        % simulated signal in degree
        s(iVar) = curDirection * [0, 0, sum( qd(iSeg_qd{curSegment}))]';
        
        % Jacobians
        if nargout > 1
            % ds_dq(iVar, :) are all zero
            ds_dqd(iVar, iSeg_qd{curSegment}) = curDirection * [0, 0, 1]'; % only one if this qd element is contained in s
            % ds_dqdd(iVar, :) are all zero
        end
    end
    
end
end
