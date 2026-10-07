%======================================================================
%> @file heelstrikeVariancePenalty.m
%> @brief Collocation function to compute smooth heelstrike variance penalty
%>
%> @author Gemini 3.5, Markus
%> @date July, 2026
%======================================================================
function output = heelstrikeVariancePenalty(obj, option, X, force_threshold, std_limit)

fctname = 'heelstrikeVariancePenalty';

%% Initialization
if strcmp(option, 'init')
    if ~isfield(obj.idx, 'states')
        error('Model states are not stored in state vector X.')
    end
    output = NaN;
    return;
end

%% Parse default arguments
if nargin < 4 || isempty(force_threshold)
    force_threshold = 0.05; % in BW
end
if nargin < 5 || isempty(std_limit)
    std_limit = 0.05; % 5.0 cm in m
end

% Constants for smoothing
eps_F = 1e-12;  % smoothing for force softplus
eps_W = 1e-3;   % smoothing for weight division
eps_std = 1e-6; % smoothing for standard deviation sqrt
eps_P = 1e-6;   % smoothing for penalty softplus

%% Extract state indexes
idxFy = obj.model.extractState('Fy');
idxxc = obj.model.extractState('xc');
idxzc = obj.model.extractState('zc');
nCPs = obj.model.nCPs;
N = obj.nNodes;

% Find if the problem is periodic symmetric
sym = 0;
iperiodicityConstraint = find(strcmp({obj.constraintTerms.name},'periodicityConstraint'), 1);
if isempty(iperiodicityConstraint)
    iperiodicityConstraint = find(strcmp({obj.constraintTerms.name},'dynamicPeriodicConstraints'), 1);
end
if ~isempty(iperiodicityConstraint)
    sym = obj.constraintTerms(iperiodicityConstraint).varargin{1};
end

% Get dur, speed, and step displacement dx, dz
dur = X(obj.idx.dur);
speed = X(obj.idx.speed);
dx = dur * speed(1);
dz = 0;
if numel(speed) > 1 && ~isempty(idxzc)
    dz = dur * speed(2);
end

% states at collocation nodes 1:N
states = X(obj.idx.states(:, 1:N));

if strcmp(option, 'objval')
    J = 0;
    
    if sym
        % Symmetric step: left foot touchdown at end of cycle continues as right foot touchdown at beginning
        nCPsHalf = nCPs / 2;
        for i = 1:nCPsHalf
            iR = i;
            iL = i + nCPsHalf;
            
            FR = states(idxFy(iR), :);
            xR = states(idxxc(iR), :);
            FL = states(idxFy(iL), :);
            xL = states(idxxc(iL), :);
            
            F_comb = [FL, FR];
            x_comb = [xL + dx, xR];
            
            F_diff = F_comb - force_threshold;
            F_active = 0.5 * (F_diff + sqrt(F_diff.^2 + eps_F));
            W = sum(F_active) + eps_W;
            
            mean_x = sum(F_active .* x_comb) / W;
            var_x = sum(F_active .* (x_comb - mean_x).^2) / W;
            std_x = sqrt(var_x + eps_std);
            
            y_x = std_x - std_limit;
            smax_x = 0.5 * (y_x + sqrt(y_x.^2 + eps_P));
            J = J + smax_x^2;
            
            if ~isempty(idxzc)
                zR = states(idxzc(iR), :);
                zL = states(idxzc(iL), :);
                z_comb = [-zL + dz, zR];
                
                mean_z = sum(F_active .* z_comb) / W;
                var_z = sum(F_active .* (z_comb - mean_z).^2) / W;
                std_z = sqrt(var_z + eps_std);
                
                y_z = std_z - std_limit;
                smax_z = 0.5 * (y_z + sqrt(y_z.^2 + eps_P));
                J = J + smax_z^2;
            end
        end
    else
        % Non-symmetric: each foot touchdown wraps to itself
        for i = 1:nCPs
            F = states(idxFy(i), :);
            x = states(idxxc(i), :);
            
            F_comb = [F, F];
            x_comb = [x + dx, x];
            
            F_diff = F_comb - force_threshold;
            F_active = 0.5 * (F_diff + sqrt(F_diff.^2 + eps_F));
            W = sum(F_active) + eps_W;
            
            mean_x = sum(F_active .* x_comb) / W;
            var_x = sum(F_active .* (x_comb - mean_x).^2) / W;
            std_x = sqrt(var_x + eps_std);
            
            y_x = std_x - std_limit;
            smax_x = 0.5 * (y_x + sqrt(y_x.^2 + eps_P));
            J = J + smax_x^2;
            
            if ~isempty(idxzc)
                z = states(idxzc(i), :);
                z_comb = [z + dz, z];
                
                mean_z = sum(F_active .* z_comb) / W;
                var_z = sum(F_active .* (z_comb - mean_z).^2) / W;
                std_z = sqrt(var_z + eps_std);
                
                y_z = std_z - std_limit;
                smax_z = 0.5 * (y_z + sqrt(y_z.^2 + eps_P));
                J = J + smax_z^2;
            end
        end
    end
    
    output = J;
    
elseif strcmp(option, 'gradient')
    output = zeros(size(X));
    dJ_dx_states = zeros(size(states));
    dJ_ddur = 0;
    dJ_dspeed = zeros(size(speed));
    
    if sym
        nCPsHalf = nCPs / 2;
        for i = 1:nCPsHalf
            iR = i;
            iL = i + nCPsHalf;
            
            FR = states(idxFy(iR), :);
            xR = states(idxxc(iR), :);
            FL = states(idxFy(iL), :);
            xL = states(idxxc(iL), :);
            
            F_comb = [FL, FR];
            x_comb = [xL + dx, xR];
            
            % --- Forward Pass ---
            F_diff = F_comb - force_threshold;
            sqrt_F_diff = sqrt(F_diff.^2 + eps_F);
            F_active = 0.5 * (F_diff + sqrt_F_diff);
            W = sum(F_active) + eps_W;
            
            mean_x = sum(F_active .* x_comb) / W;
            D_x = x_comb - mean_x;
            var_x = sum(F_active .* D_x.^2) / W;
            std_x = sqrt(var_x + eps_std);
            y_x = std_x - std_limit;
            sqrt_y_x = sqrt(y_x.^2 + eps_P);
            smax_x = 0.5 * (y_x + sqrt_y_x);
            
            % --- Backward Pass ---
            dF_active_dF = 0.5 * (1 + F_diff ./ sqrt_F_diff);
            S_x = (mean_x * eps_W) / W;
            dJ_dy_x = smax_x * (1 + y_x / sqrt_y_x);
            dstd_dvar_x = 0.5 / std_x;
            C_x = dJ_dy_x * dstd_dvar_x / W;
            
            dJ_dx_pos_comb = C_x * 2 * F_active .* (D_x - S_x);
            dJ_dF_active = C_x * (D_x.^2 - 2 * D_x * S_x - var_x);
            
            % Gradient with respect to displacement
            dJ_dxL = dJ_dx_pos_comb(1:N);
            dJ_dxR = dJ_dx_pos_comb(N+1:2*N);
            dJ_ddx = sum(dJ_dxL);
            
            if ~isempty(idxzc)
                zR = states(idxzc(iR), :);
                zL = states(idxzc(iL), :);
                z_comb = [-zL + dz, zR];
                
                mean_z = sum(F_active .* z_comb) / W;
                D_z = z_comb - mean_z;
                var_z = sum(F_active .* D_z.^2) / W;
                std_z = sqrt(var_z + eps_std);
                y_z = std_z - std_limit;
                sqrt_y_z = sqrt(y_z.^2 + eps_P);
                smax_z = 0.5 * (y_z + sqrt_y_z);
                
                S_z = (mean_z * eps_W) / W;
                dJ_dy_z = smax_z * (1 + y_z / sqrt_y_z);
                dstd_dvar_z = 0.5 / std_z;
                C_z = dJ_dy_z * dstd_dvar_z / W;
                
                dJ_dz_pos_comb = C_z * 2 * F_active .* (D_z - S_z);
                dJ_dF_active_z = C_z * (D_z.^2 - 2 * D_z * S_z - var_z);
                
                dJ_dF_active = dJ_dF_active + dJ_dF_active_z;
                
                dJ_dzL = dJ_dz_pos_comb(1:N);
                dJ_dzR = dJ_dz_pos_comb(N+1:2*N);
                dJ_ddz = sum(dJ_dzL);
                
                dJ_dx_states(idxzc(iR), :) = dJ_dx_states(idxzc(iR), :) + dJ_dzR;
                dJ_dx_states(idxzc(iL), :) = dJ_dx_states(idxzc(iL), :) - dJ_dzL; % note negative sign for mirroring
                
                dJ_ddur = dJ_ddur + dJ_ddz * speed(2);
                dJ_dspeed(2) = dJ_dspeed(2) + dJ_ddz * dur;
            end
            
            dJ_dF_comb = dJ_dF_active .* dF_active_dF;
            
            dJ_dx_states(idxxc(iR), :) = dJ_dx_states(idxxc(iR), :) + dJ_dxR;
            dJ_dx_states(idxxc(iL), :) = dJ_dx_states(idxxc(iL), :) + dJ_dxL;
            dJ_dx_states(idxFy(iR), :) = dJ_dx_states(idxFy(iR), :) + dJ_dF_comb(N+1:2*N);
            dJ_dx_states(idxFy(iL), :) = dJ_dx_states(idxFy(iL), :) + dJ_dF_comb(1:N);
            
            dJ_ddur = dJ_ddur + dJ_ddx * speed(1);
            dJ_dspeed(1) = dJ_dspeed(1) + dJ_ddx * dur;
        end
    else
        % Non-symmetric
        for i = 1:nCPs
            F = states(idxFy(i), :);
            x = states(idxxc(i), :);
            
            F_comb = [F, F];
            x_comb = [x + dx, x];
            
            % --- Forward Pass ---
            F_diff = F_comb - force_threshold;
            sqrt_F_diff = sqrt(F_diff.^2 + eps_F);
            F_active = 0.5 * (F_diff + sqrt_F_diff);
            W = sum(F_active) + eps_W;
            
            mean_x = sum(F_active .* x_comb) / W;
            D_x = x_comb - mean_x;
            var_x = sum(F_active .* D_x.^2) / W;
            std_x = sqrt(var_x + eps_std);
            y_x = std_x - std_limit;
            sqrt_y_x = sqrt(y_x.^2 + eps_P);
            smax_x = 0.5 * (y_x + sqrt_y_x);
            
            % --- Backward Pass ---
            dF_active_dF = 0.5 * (1 + F_diff ./ sqrt_F_diff);
            S_x = (mean_x * eps_W) / W;
            dJ_dy_x = smax_x * (1 + y_x / sqrt_y_x);
            dstd_dvar_x = 0.5 / std_x;
            C_x = dJ_dy_x * dstd_dvar_x / W;
            
            dJ_dx_pos_comb = C_x * 2 * F_active .* (D_x - S_x);
            dJ_dF_active = C_x * (D_x.^2 - 2 * D_x * S_x - var_x);
            
            dJ_dxL = dJ_dx_pos_comb(1:N);
            dJ_dxR = dJ_dx_pos_comb(N+1:2*N);
            dJ_ddx = sum(dJ_dxL);
            
            if ~isempty(idxzc)
                z = states(idxzc(i), :);
                z_comb = [z + dz, z];
                
                mean_z = sum(F_active .* z_comb) / W;
                D_z = z_comb - mean_z;
                var_z = sum(F_active .* D_z.^2) / W;
                std_z = sqrt(var_z + eps_std);
                y_z = std_z - std_limit;
                sqrt_y_z = sqrt(y_z.^2 + eps_P);
                smax_z = 0.5 * (y_z + sqrt_y_z);
                
                S_z = (mean_z * eps_W) / W;
                dJ_dy_z = smax_z * (1 + y_z / sqrt_y_z);
                dstd_dvar_z = 0.5 / std_z;
                C_z = dJ_dy_z * dstd_dvar_z / W;
                
                dJ_dz_pos_comb = C_z * 2 * F_active .* (D_z - S_z);
                dJ_dF_active_z = C_z * (D_z.^2 - 2 * D_z * S_z - var_z);
                
                dJ_dF_active = dJ_dF_active + dJ_dF_active_z;
                
                dJ_dzL = dJ_dz_pos_comb(1:N);
                dJ_dzR = dJ_dz_pos_comb(N+1:2*N);
                dJ_ddz = sum(dJ_dzL);
                
                dJ_dx_states(idxzc(i), :) = dJ_dx_states(idxzc(i), :) + dJ_dzL + dJ_dzR;
                
                dJ_ddur = dJ_ddur + dJ_ddz * speed(2);
                dJ_dspeed(2) = dJ_dspeed(2) + dJ_ddz * dur;
            end
            
            dJ_dF_comb = dJ_dF_active .* dF_active_dF;
            
            dJ_dx_states(idxxc(i), :) = dJ_dx_states(idxxc(i), :) + dJ_dxL + dJ_dxR;
            dJ_dx_states(idxFy(i), :) = dJ_dx_states(idxFy(i), :) + dJ_dF_comb(1:N) + dJ_dF_comb(N+1:2*N);
            
            dJ_ddur = dJ_ddur + dJ_ddx * speed(1);
            dJ_dspeed(1) = dJ_dspeed(1) + dJ_ddx * dur;
        end
    end
    
    % Map back to the full state vector X
    output(obj.idx.states(:, 1:N)) = dJ_dx_states;
    output(obj.idx.dur) = dJ_ddur;
    output(obj.idx.speed) = dJ_dspeed;
    
else
    error('Unknown option');
end

end
