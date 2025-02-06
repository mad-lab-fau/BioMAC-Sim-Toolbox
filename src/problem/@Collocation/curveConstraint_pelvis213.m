%======================================================================
%> @file curveConstraint_pelvis213.m
%> @brief Collocation function to compute the periodicity constraint for curved running
%> @details
%> Details: Collocation::curveConstraint_pelvis213()
%>
%> @author Eva Dorschky, Marlies Nitschke
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Function to compute the periodicity constraint for curved running
%>
%> @details 
%> It ensures a perdiodic movement on a circle with a specifc radius. 
%> It was defined for our 3D model with the pelvis rotation sequence 213.
%>
%> It uses the following equations:
%> - \f$ c = \vert\vert v \vert\vert \cdot T \f$
%> - \f$ \theta = 2 \arcsin \left( \frac{c}{2r} \right) \f$
%>
%> @param obj           Collocation class object
%> @param option        String: Parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%> @param radius        Double: Radius of the circle
%======================================================================
function output = curveConstraint_pelvis213(obj,option,X,radius)
%% check input parameter
if ~isa(obj.model,'Gait3d')
    error('Problem:curveConstraint_pelvis213','Model must be a Gait3d model for curve running')
end
if ~strcmp(obj.model.osim.name,'3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence') && ~strcmp(obj.model.osim.name, 'gait14dof22musc and Pelvis Rotation-Obliquity-Tilt Sequence') 
    error('Problem:curveConstraint_pelvis213','OpenSim model name must be: 3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence')
end
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'speed') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('Problem:curveConstraint_pelvis213','State vector X does not contain required fields. It requires states, controls, speed or dur.')
end
if size(X(obj.idx.states),2) ~= (obj.nNodes+1) || size(X(obj.idx.controls),2) ~= (obj.nNodes+1)
    error('Problem:curveConstraint_pelvis213','Model states and controls need to be optimized at collocation nodes + 1')
end


%% compute demanded output
nStates = size(obj.idx.states,1);
nControls =  size(obj.idx.controls,1);
icx = 1:nStates; % first indices are for the states
icu = (nStates+1):(nStates+nControls); % next indices are for the controls

% compute angle theta of circular segment
speed = X(obj.idx.speed);
duration = X(obj.idx.dur);
chord = norm(speed) * duration; % with v_x = deltaX / T and v_z = deltaZ / T
asinArgument = chord/(2*radius);
theta = 2 * asin(asinArgument); % angle in rad = arcLength / radius with arcLength = 2 * radius * asin(chord/radius) (See wikipedia)

% Get the indices for all types in the states
idx_q    = obj.model.extractState('q');
idx_qdot = obj.model.extractState('qdot');
idx_a    = obj.model.extractState('a');
idx_s    = obj.model.extractState('s');
idx_xf    = obj.model.extractState('xf');
idx_yf    = obj.model.extractState('yf');
idx_zf    = obj.model.extractState('zf');
idx_xc    = obj.model.extractState('xc');
idx_yc    = obj.model.extractState('yc');
idx_zc    = obj.model.extractState('zc');
idx_Fx    = obj.model.extractState('Fx');
idx_Fy    = obj.model.extractState('Fy');
idx_Fz    = obj.model.extractState('Fz');

% Get indices for specific rows of the states
idx_q_tx    = obj.model.extractState('q','pelvis_tx');
idx_q_tz    = obj.model.extractState('q','pelvis_tz');
idx_q_rot   = obj.model.extractState('q','pelvis_rotation');
idx_qdot_tx = obj.model.extractState('qdot','pelvis_tx');
idx_qdot_tz = obj.model.extractState('qdot','pelvis_tz');

if strcmp(option,'confun') %constraints of periodicity constraint
    output = zeros(nStates+nControls,1);
    
    % q
    output(icx(idx_q))     = X(obj.idx.states(idx_q,end)) - X(obj.idx.states(idx_q,1)); % in general equal
    output(icx(idx_q_tx))  = X(obj.idx.states(idx_q_tx,end))  - cos(theta) * X(obj.idx.states(idx_q_tx,1)) - sin(theta) * X(obj.idx.states(idx_q_tz,1)); % rotating p = (p_x, p_z)
    output(icx(idx_q_tz))  = X(obj.idx.states(idx_q_tz,end))  + sin(theta) * X(obj.idx.states(idx_q_tx,1)) - cos(theta) * X(obj.idx.states(idx_q_tz,1)); % rotating p = (p_x, p_z)
    output(icx(idx_q_rot)) = X(obj.idx.states(idx_q_rot,end)) - X(obj.idx.states(idx_q_rot,1)) - theta; % pelvis rotation
    
    % qdot
    output(icx(idx_qdot))    = X(obj.idx.states(idx_qdot,end)) - X(obj.idx.states(idx_qdot,1)); % in general equal
    output(icx(idx_qdot_tx)) = X(obj.idx.states(idx_qdot_tx,end)) - cos(theta) * X(obj.idx.states(idx_qdot_tx,1)) - sin(theta) * X(obj.idx.states(idx_qdot_tz,1)); % rotating v = (v_x, v_z)
    output(icx(idx_qdot_tz)) = X(obj.idx.states(idx_qdot_tz,end)) + sin(theta) * X(obj.idx.states(idx_qdot_tx,1)) - cos(theta) * X(obj.idx.states(idx_qdot_tz,1)); % rotating v = (v_x, v_z)
    
    % a, s
    output(icx(idx_a)) = X(obj.idx.states(idx_a,end)) - X(obj.idx.states(idx_a,1)); % equal
    output(icx(idx_s)) = X(obj.idx.states(idx_s,end)) - X(obj.idx.states(idx_s,1)); % equal
    
    % xf, yf, zf
    output(icx(idx_xf)) = X(obj.idx.states(idx_xf,end)) - X(obj.idx.states(idx_xf,1)); % equal
    output(icx(idx_yf)) = X(obj.idx.states(idx_yf,end)) - X(obj.idx.states(idx_yf,1)); % equal
    output(icx(idx_zf)) = X(obj.idx.states(idx_zf,end)) - X(obj.idx.states(idx_zf,1)); % equal
    
    % xc, yc, zc
    output(icx(idx_xc)) = X(obj.idx.states(idx_xc,end)) - X(obj.idx.states(idx_xc,1)); % in general equal
    output(icx(idx_yc)) = X(obj.idx.states(idx_yc,end)) - X(obj.idx.states(idx_yc,1)); % in general equal
    output(icx(idx_zc)) = X(obj.idx.states(idx_zc,end)) - X(obj.idx.states(idx_zc,1)); % in general equal
    output(icx(idx_xc)) = X(obj.idx.states(idx_xc,end))  - cos(theta) * X(obj.idx.states(idx_xc,1)) - sin(theta) * X(obj.idx.states(idx_zc,1)); % rotating p = (p_x, p_z)
    output(icx(idx_zc)) = X(obj.idx.states(idx_zc,end))  + sin(theta) * X(obj.idx.states(idx_xc,1)) - cos(theta) * X(obj.idx.states(idx_zc,1)); % rotating p = (p_x, p_z)
    
    % Fx, Fy, Fz
    output(icx(idx_Fy)) = X(obj.idx.states(idx_Fy,end)) - X(obj.idx.states(idx_Fy,1)); % vertical force is equal
    output(icx(idx_Fx)) = X(obj.idx.states(idx_Fx,end)) - cos(theta) * X(obj.idx.states(idx_Fx,1)) - sin(theta) * X(obj.idx.states(idx_Fz,1)); % rotating F_horizontal = (F_x, F_z)
    output(icx(idx_Fz)) = X(obj.idx.states(idx_Fz,end)) + sin(theta) * X(obj.idx.states(idx_Fx,1)) - cos(theta) * X(obj.idx.states(idx_Fz,1)); % rotating F_horizontal = (F_x, F_z)

    % controls
    output(icu) = X(obj.idx.controls(:,end)) - X(obj.idx.controls(:,1));
    
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    
    output = spalloc(nStates+nControls,length(X),obj.Jnnz);
    
    theta_d_duration = 1 / sqrt(1 - asinArgument^2) / radius * norm(speed);
    theta_d_speed    = 1 / sqrt(1 - asinArgument^2) / radius / norm(speed) * duration .*speed;
    
    % q
    output(icx(idx_q),obj.idx.states(idx_q,end)) =  speye(numel(idx_q));
    output(icx(idx_q),obj.idx.states(idx_q,1))   = -speye(numel(idx_q));
    
    output(icx(idx_q_tx),obj.idx.states(idx_q_tx,1))   = -cos(theta);
    output(icx(idx_q_tx),obj.idx.states(idx_q_tz,1))   = -sin(theta);
    output(icx(idx_q_tx),obj.idx.dur)   = ( sin(theta) * X(obj.idx.states(idx_q_tx,1)) - cos(theta) * X(obj.idx.states(idx_q_tz,1)) ) * theta_d_duration;
    output(icx(idx_q_tx),obj.idx.speed) = ( sin(theta) * X(obj.idx.states(idx_q_tx,1)) - cos(theta) * X(obj.idx.states(idx_q_tz,1)) ) * theta_d_speed;
    
    output(icx(idx_q_tz),obj.idx.states(idx_q_tx,1))   =  sin(theta);
    output(icx(idx_q_tz),obj.idx.states(idx_q_tz,1))   = -cos(theta);
    output(icx(idx_q_tz),obj.idx.dur)   = ( cos(theta) * X(obj.idx.states(idx_q_tx,1)) + sin(theta) * X(obj.idx.states(idx_q_tz,1)) ) * theta_d_duration;
    output(icx(idx_q_tz),obj.idx.speed) = ( cos(theta) * X(obj.idx.states(idx_q_tx,1)) + sin(theta) * X(obj.idx.states(idx_q_tz,1)) ) * theta_d_speed;
    
    output(icx(idx_q_rot),obj.idx.dur)   = -theta_d_duration;
    output(icx(idx_q_rot),obj.idx.speed) = -theta_d_speed;
    
    % qdot
    output(icx(idx_qdot),obj.idx.states(idx_qdot,end)) =  speye(numel(idx_qdot));
    output(icx(idx_qdot),obj.idx.states(idx_qdot,1))   = -speye(numel(idx_qdot));
    
    output(icx(idx_qdot_tx),obj.idx.states(idx_qdot_tx,1))   = -cos(theta);
    output(icx(idx_qdot_tx),obj.idx.states(idx_qdot_tz,1))   = -sin(theta);
    output(icx(idx_qdot_tx),obj.idx.dur)   = ( sin(theta) * X(obj.idx.states(idx_qdot_tx,1)) - cos(theta) * X(obj.idx.states(idx_qdot_tz,1)) ) * theta_d_duration;
    output(icx(idx_qdot_tx),obj.idx.speed) = ( sin(theta) * X(obj.idx.states(idx_qdot_tx,1)) - cos(theta) * X(obj.idx.states(idx_qdot_tz,1)) ) * theta_d_speed;
    
    output(icx(idx_qdot_tz),obj.idx.states(idx_qdot_tx,1))   =  sin(theta);
    output(icx(idx_qdot_tz),obj.idx.states(idx_qdot_tz,1))   = -cos(theta);
    output(icx(idx_qdot_tz),obj.idx.dur)   = ( cos(theta) * X(obj.idx.states(idx_qdot_tx,1)) + sin(theta) * X(obj.idx.states(idx_qdot_tz,1)) ) * theta_d_duration;
    output(icx(idx_qdot_tz),obj.idx.speed) = ( cos(theta) * X(obj.idx.states(idx_qdot_tx,1)) + sin(theta) * X(obj.idx.states(idx_qdot_tz,1)) ) * theta_d_speed;
    
    % a
    output(icx(idx_a),obj.idx.states(idx_a,end)) =  speye(numel(idx_a));
    output(icx(idx_a),obj.idx.states(idx_a,1))   = -speye(numel(idx_a));
    
    % s
    output(icx(idx_s),obj.idx.states(idx_s,end)) =  speye(numel(idx_s));
    output(icx(idx_s),obj.idx.states(idx_s,1))   = -speye(numel(idx_s));
    
    % xf, yf, zf
    output(icx(idx_xf),obj.idx.states(idx_xf,end)) =  speye(numel(idx_xf));
    output(icx(idx_xf),obj.idx.states(idx_xf,1))   = -speye(numel(idx_xf));
    
    output(icx(idx_yf),obj.idx.states(idx_yf,end)) =  speye(numel(idx_yf));
    output(icx(idx_yf),obj.idx.states(idx_yf,1))   = -speye(numel(idx_yf));
   
    output(icx(idx_zf),obj.idx.states(idx_zf,end)) =  speye(numel(idx_zf));
    output(icx(idx_zf),obj.idx.states(idx_zf,1))   = -speye(numel(idx_zf));
    
    % xc, yc, zc   
    output(icx(idx_yc),obj.idx.states(idx_yc,end)) =  speye(numel(idx_yc));
    output(icx(idx_yc),obj.idx.states(idx_yc,1))   = -speye(numel(idx_yc));
   
    output(icx(idx_xc),obj.idx.states(idx_xc,end)) =  speye(numel(idx_xc));
    output(icx(idx_xc),obj.idx.states(idx_xc,1))   = -speye(numel(idx_xc));
    output(icx(idx_xc),obj.idx.states(idx_xc,1))   = -cos(theta)*speye(numel(idx_xc));
    output(icx(idx_xc),obj.idx.states(idx_zc,1))   = -sin(theta)*speye(numel(idx_xc));
    output(icx(idx_xc),obj.idx.dur)      = ( sin(theta) * X(obj.idx.states(idx_xc,1)) - cos(theta) * X(obj.idx.states(idx_zc,1)) ) * theta_d_duration;
    for iSpeed = 1 : numel(obj.idx.speed) % there is only one entry if already the norm is given and not v = [v_x, v_z]
        output(icx(idx_xc),obj.idx.speed(iSpeed)) = ( sin(theta) * X(obj.idx.states(idx_xc,1)) - cos(theta) * X(obj.idx.states(idx_zc,1)) ) * theta_d_speed(iSpeed);
    end
    
    output(icx(idx_zc),obj.idx.states(idx_zc,end)) =  speye(numel(idx_zc));
    output(icx(idx_zc),obj.idx.states(idx_zc,1))   = -speye(numel(idx_zc));
    output(icx(idx_zc),obj.idx.states(idx_xc,1))   =  sin(theta)*speye(numel(idx_zc));
    output(icx(idx_zc),obj.idx.states(idx_zc,1))   = -cos(theta)*speye(numel(idx_zc));
    output(icx(idx_zc),obj.idx.dur)      = ( cos(theta) * X(obj.idx.states(idx_xc,1)) + sin(theta) * X(obj.idx.states(idx_zc,1)) ) * theta_d_duration;
    for iSpeed = 1 : numel(obj.idx.speed) % there is only one entry if already the norm is given and not v = [v_x, v_z]
        output(icx(idx_zc),obj.idx.speed(iSpeed)) = ( cos(theta) * X(obj.idx.states(idx_xc,1)) + sin(theta) * X(obj.idx.states(idx_zc,1)) ) * theta_d_speed(iSpeed);
    end    
    
    % Fx, Fy, Fz
    output(icx(idx_Fy),obj.idx.states(idx_Fy,end)) =  speye(numel(idx_Fy));
    output(icx(idx_Fy),obj.idx.states(idx_Fy,1))   = -speye(numel(idx_Fy));
    
    output(icx(idx_Fx),obj.idx.states(idx_Fx,end)) =  speye(numel(idx_Fx));
    output(icx(idx_Fx),obj.idx.states(idx_Fx,1))   = -speye(numel(idx_Fx));
    output(icx(idx_Fx),obj.idx.states(idx_Fx,1))   = -cos(theta)*speye(numel(idx_Fx));
    output(icx(idx_Fx),obj.idx.states(idx_Fz,1))   = -sin(theta)*speye(numel(idx_Fx));
    output(icx(idx_Fx),obj.idx.dur)      = ( sin(theta) * X(obj.idx.states(idx_Fx,1)) - cos(theta) * X(obj.idx.states(idx_Fz,1)) ) * theta_d_duration;
    for iSpeed = 1 : numel(obj.idx.speed) % there is only one entry if already the norm is given and not v = [v_x, v_z]
        output(icx(idx_Fx),obj.idx.speed(iSpeed)) = ( sin(theta) * X(obj.idx.states(idx_Fx,1)) - cos(theta) * X(obj.idx.states(idx_Fz,1)) ) * theta_d_speed(iSpeed);
    end    
    
    output(icx(idx_Fz),obj.idx.states(idx_Fz,end)) =  speye(numel(idx_Fz));
    output(icx(idx_Fz),obj.idx.states(idx_Fz,1))   = -speye(numel(idx_Fz));
    output(icx(idx_Fz),obj.idx.states(idx_Fx,1))   =  sin(theta)*speye(numel(idx_Fz));
    output(icx(idx_Fz),obj.idx.states(idx_Fz,1))   = -cos(theta)*speye(numel(idx_Fz));
    output(icx(idx_Fz),obj.idx.dur)      = ( cos(theta) * X(obj.idx.states(idx_Fx,1)) + sin(theta) * X(obj.idx.states(idx_Fz,1)) ) * theta_d_duration;
    for iSpeed = 1 : numel(obj.idx.speed) % there is only one entry if already the norm is given and not v = [v_x, v_z]
        output(icx(idx_Fz),obj.idx.speed(iSpeed)) = ( cos(theta) * X(obj.idx.states(idx_Fx,1)) + sin(theta) * X(obj.idx.states(idx_Fz,1)) ) * theta_d_speed(iSpeed);
    end   
    
    % controls
    output(icu,obj.idx.controls(:,end)) =  speye(nControls);
    output(icu,obj.idx.controls(:,1))   = -speye(nControls);
    
else
    error('Unknown option');
end
end

