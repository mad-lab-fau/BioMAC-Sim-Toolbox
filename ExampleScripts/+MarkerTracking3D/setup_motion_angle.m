%======================================================================
%> @file +MarkerTracking3D/setup_motion_angle.m
%> @brief Function to setup motion problem with angle tracking
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================

% ======================================================================
%> @brief Function to setup motion problem with angle tracking
%>
%> @param   model          Gait3d: Model used for simulation
%> @param   dataInvFile    String: Filename containg data struct for angles for tracking
%> @param   dataMeasFile   String: Filename containg data struct for GRFs for tracking
%> @param   initialFile    String: Result used for initial guess
%> @param   resultFile     String: Filename to log files
%> @param   W              Struct: Weights for objective terms
%> @param   nNodesEarlier  Double: Number of samples which we start earlier
%> @retval  problem        Collocation: Generated standing problem
% ======================================================================
function problem = setup_motion_angle(model,dataInvFile,dataMeasFile,initialFile,resultFile,W, nNodesEarlier)

% Load tracking data
trackingDataInv = TrackingData.loadStruct(dataInvFile);
indicesStartEnd = trackingDataInv.movementEvents.index(strcmp(trackingDataInv.movementEvents.name, 'R_IC'));
trackingDataInv.trimData(indicesStartEnd(1)-nNodesEarlier, indicesStartEnd(2)-1);

trackingDataMeas = TrackingData.loadStruct(dataMeasFile);
indicesStartEnd = trackingDataMeas.movementEvents.index(strcmp(trackingDataMeas.movementEvents.name, 'R_IC'));
trackingDataMeas.trimData(indicesStartEnd(1)-nNodesEarlier, indicesStartEnd(2)-1);
N = trackingDataMeas.nSamples;

% Create problem
Euler = 'BE';
plotLog = 0;
problem = Collocation(model,N,Euler,resultFile,plotLog);

% Add variables which are optimized
states_min = repmat(model.states.xmin,1,N);
states_max = repmat(model.states.xmax,1,N);
idxCPxc = model.extractState('xc');
idxCPzc = model.extractState('zc');
p_global_x = [-5, 5];
p_global_z = [-5, 5];
states_min(idxCPxc, :) = p_global_x(1);
states_max(idxCPxc, :) = p_global_x(2);
states_min(idxCPzc, :) = p_global_z(1);
states_max(idxCPzc, :) = p_global_z(2);
problem.addOptimVar('states',states_min,states_max);
problem.addOptimVar('controls',repmat(model.controls.xmin,1,N), repmat(model.controls.xmax,1,N));
h = 1/175; % 175 Hz
targetdur =  h*(N-1);
problem.addOptimVar('dur',targetdur,targetdur);

% Initialize the problem with an old result specified in initialFile
problem.makeinitialguess(initialFile); 

% Add tracking terms
trackingDataTrans = trackingDataInv.extractData('translation', {'pelvis_tx', 'pelvis_ty', 'pelvis_tz'});
trackingDataAng = trackingDataInv.extractData('angle', {'pelvis_rotation'    'pelvis_obliquity'    'pelvis_tilt' ...
    'hip_flexion_r'    'hip_adduction_r'    'hip_rotation_r'    'knee_angle_r'    'ankle_angle_r'    'subtalar_angle_r' ...
    'mtp_angle_r'    'hip_flexion_l'    'hip_adduction_l'    'hip_rotation_l'    'knee_angle_l'    'ankle_angle_l'    ...
    'subtalar_angle_l'    'mtp_angle_l'    'lumbar_extension'   'lumbar_bending'    'lumbar_rotation'    'arm_flex_r' ...
    'arm_add_r'    'arm_rot_r'    'elbow_flex_r'    'pro_sup_r'    'arm_flex_l'    'arm_add_l'    'arm_rot_l'    'elbow_flex_l'    'pro_sup_l'});
trackingDataGRF = trackingDataMeas.extractData('GRF', {'GRF_x_r', 'GRF_y_r', 'GRF_z_r', 'GRF_x_l', 'GRF_y_l', 'GRF_z_l'});
problem.addObjective(@trackTranslations,W.trackTrans,trackingDataTrans);
problem.addObjective(@trackAngles,W.trackAngle,trackingDataAng);
problem.addObjective(@trackGRF,W.trackGRF,trackingDataGRF);
% Add effort terms
speedWeighting = 0;
problem.addObjective(@effortTermMuscles,W.effMuscles,'volumeweighted',3,speedWeighting);
problem.addObjective(@effortTermTorques,W.effTorques,2, speedWeighting)
% add regularization term
problem.addObjective(@regTerm,W.reg)

% Add constraints
problem.addConstraint(@dynamicsFirstNodeConstraint,model.constraints.fmin,model.constraints.fmax)
problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,N-1),repmat(model.constraints.fmax,1,N-1))

end
