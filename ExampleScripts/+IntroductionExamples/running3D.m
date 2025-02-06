%======================================================================
%> @file IntroductionExamples/running3D.m
%> @brief Function to specify the optimization problem for 3D running
%>
%> @author Anne Koelewijn, Marlies Nitschke
%> @date July, 2024
%======================================================================

% ======================================================================
%> @brief Function to run standing problem
%>
%> @param   model       Gait3d: Model used for simulation
%> @param   resultFile  String: Filename to log files
%> @retval  problem     Collocation: Generated standing problem
% ======================================================================
% ======================================================================
%> @brief Function to setup running problem
%>
%> @param   model          Gait3d: Model used for simulation
%> @param   dataFile       String: Filename containg data struct for tracking
%> @param   initialFile    String: Result used for initial guess
%> @param   resultFile     String: Filename to log files
%> @param   N              Double: Number of nodes
%> @param   sym            Boolean: Symmetry of the movment
%> @param   Euler          String: Euler discretization method
%> @param   W              Struct: Weights for objective terms
%> @param   idxTracking    Double: Index specifying which variables are tracked
%> @retval  problem        Collocation: Generated standing problem
% ======================================================================
function problem = running3D(model,trackingData,initialFile,resultFile,N,sym,W, targetspeed_x, targetspeed_z, targetdur)

% 1. Create problem
plotLog = 0;
problem = Collocation(model,N,'BE',resultFile,plotLog);

% 2. Change bounds of model
states_min = repmat(model.states.xmin,1,N+1);
states_max = repmat(model.states.xmax,1,N+1);
idxPelvisX = model.extractState('q', 'pelvis_tx');
idxPelvisY = model.extractState('q', 'pelvis_ty');
idxPelvisZ = model.extractState('q', 'pelvis_tz');
% Bounds for global position of pelvis and contact points
p_global_x = [-5, 5];
p_global_y = [-1, 2];
p_global_z = [-5, 5];
states_min(idxPelvisX, :) = p_global_x(1);
states_max(idxPelvisX, :) = p_global_x(2);
states_min(idxPelvisY, :) = p_global_y(1);
states_max(idxPelvisY, :) = p_global_y(2);
states_min(idxPelvisZ, :) = p_global_z(1);
states_max(idxPelvisZ, :) = p_global_z(2);
states_min(model.extractState('xc'), :) = p_global_x(1);
states_max(model.extractState('xc'), :) = p_global_x(2);
states_min(model.extractState('yc'), :) = p_global_y(1);
states_max(model.extractState('yc'), :) = p_global_y(2);
states_min(model.extractState('zc'), :) = p_global_z(1);
states_max(model.extractState('zc'), :) = p_global_z(2);
% Bounds for global coordinates: Start at (x, z) = (0, 0)
states_min(idxPelvisX,1) = 0;
states_max(idxPelvisX,1) = 0;
states_min(idxPelvisZ,1) = 0;
states_max(idxPelvisZ,1) = 0;
% Bounds for the arms
states_min(model.extractState('q', {'arm_flex_r', 'arm_flex_l'}), :)        = -40 /180*pi; % some may be needed
states_max(model.extractState('q', {'arm_flex_r', 'arm_flex_l'}), :)        =  40 /180*pi;
states_min(model.extractState('q', {'arm_add_r', 'arm_add_l'}), :)          = -40 /180*pi; % some may be needed
states_max(model.extractState('q', {'arm_add_r', 'arm_add_l'}), :)          =  40 /180*pi;
states_min(model.extractState('q', {'arm_rot_r', 'arm_rot_l'}), :)          = -40 /180*pi; % some may be needed
states_max(model.extractState('q', {'arm_rot_r', 'arm_rot_l'}), :)          =  40 /180*pi;
states_min(model.extractState('q', {'elbow_flex_r', 'elbow_flex_l'}), :)    = - 0 /180*pi; % some may be needed
states_max(model.extractState('q', {'elbow_flex_r', 'elbow_flex_l'}), :)    = 150 /180*pi;
states_min(model.extractState('q', {'pro_sup_r', 'pro_sup_l'}), :)          = - 0 /180*pi; % some may be needed
states_max(model.extractState('q', {'pro_sup_r', 'pro_sup_l'}), :)          = 150 /180*pi;

% 3. Add states and controls to the collocation problem
% Add states using the adapted bounds from step 2
problem.addOptimVar('states',states_min,states_max);
% Add controls
problem.addOptimVar('controls',model.controls.xmin,model.controls.xmax); 
problem.addOptimVar('controls',repmat(model.controls.xmin,1,N+1), repmat(model.controls.xmax,1,N+1));

% 4. Add duration to the collocation problem 
if isinf(W.dur)
    if sym
        problem.addOptimVar('dur',targetdur/2,targetdur/2);
    else
        problem.addOptimVar('dur',targetdur,targetdur);
    end
else
    if sym
        problem.addOptimVar('dur',0.1,1);
    else
        problem.addOptimVar('dur',0.2,2);
    end
end

% 5. Add speed to the collocation problem
if isinf(W.speed)
    problem.addOptimVar('speed',[targetspeed_x; targetspeed_z],[targetspeed_x; targetspeed_z]);
else
    problem.addOptimVar('speed',[-5;-5],[5;5], [1;1]); % Do not initialize with 0 to avoid errors in effortTermMuscles.m
end

% 6. Initialize the problem with an old result specified in initialFile
problem.makeinitialguess(initialFile); 

% 7. Add objective terms
if W.reg ~= 0
    problem.addObjective(@regTerm,W.reg)
end
if W.effMuscles ~= 0
    speedWeighting = 1;
    problem.addObjective(@effortTermMuscles,W.effMuscles,'volumeweighted',3,speedWeighting);
end
if W.effTorques ~= 0
    speedWeighting = 0;
    problem.addObjective(@effortTermTorques,W.effTorques,2, speedWeighting)
end
if W.track ~= 0
    dofNames = {'pelvis_rotation'    'pelvis_obliquity'    'pelvis_tilt'  'hip_flexion_r'    'hip_adduction_r'    'hip_rotation_r'    'knee_angle_r'    'ankle_angle_r'    'subtalar_angle_r' ...
         'mtp_angle_r'    'hip_flexion_l'    'hip_adduction_l'    'hip_rotation_l'    'knee_angle_l'    'ankle_angle_l'    'subtalar_angle_l'    'mtp_angle_l'    'lumbar_extension'    ...
         'lumbar_bending'    'lumbar_rotation'    'arm_flex_r'    'arm_add_r'    'arm_rot_r'    'elbow_flex_r'    'pro_sup_r'    'arm_flex_l'    'arm_add_l'    'arm_rot_l'    'elbow_flex_l'    'pro_sup_l'};
    grfNames = {'GRF_x_r', 'GRF_y_r', 'GRF_z_r', 'GRF_x_l', 'GRF_y_l', 'GRF_z_l'};
    trackingDataAngles = trackingData.extractData('angle', dofNames);
    trackingDataGRF    = trackingData.extractData('GRF',   grfNames);
    wAng = 1;
    wGRF = 5; % wGRF/wAng is the ratio of the weighting between GRFs and angles
    NAng = trackingDataAngles.nVariables;
    NGRF = trackingDataGRF.nVariables;
    wAngFac = wAng * NAng / (wGRF*NGRF + wAng*NAng); % We do a weighted mean and we have to compensate for 1/NAng which is done in the objctive
    wGRFFac = wGRF * NGRF / (wGRF*NGRF + wAng*NAng); % We do a weighted mean and we have to compensate for 1/NGRF which is done in the objctive
    problem.addObjective(@trackAngles,W.track*wAngFac,trackingDataAngles);
    problem.addObjective(@trackGRF,W.track*wGRFFac,trackingDataGRF);
end

% 8. Add constraints
problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,N),repmat(model.constraints.fmax,1,N))
problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),sym)

end