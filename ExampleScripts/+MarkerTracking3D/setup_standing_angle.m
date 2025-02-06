%======================================================================
%> @file +MarkerTracking3D/setup_standing_angle.m
%> @brief Function to setup standing problem with angle tracking
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================

% ======================================================================
%> @brief Function to setup standing problem with angle tracking
%>
%> @param   model       Gait3d: Model used for simulation
%> @param   dataFileInv String: Filename containg data struct for angle tracking
%> @param   dataFileMea String: Filename containg data struct for GRF tracking
%> @param   iSample     Int: Index of sample used from measured data
%> @param   resultFile  String: Filename to log files
%> @param   W           Struct: Weights for objective terms
%> @retval  problem     Collocation: Generated standing problem
% ======================================================================
function problem = setup_standing_angle(model, dataFileInv, dataFileMea, iSample, resultFile, W)

% 1. Create problem
problem = Collocation(model,1,'BE',resultFile); % simulate only one collocation node for standing, 'BE' or 'ME' doesn't matter

% 2. Get bounds and initial guess for optimization variables
% Define bounds
states_min = model.states.xmin;
states_max = model.states.xmax;
idxCPxc = model.extractState('xc');
idxCPzc = model.extractState('zc');
p_global_x = [-5, 5];
p_global_z = [-5, 5];
states_min(idxCPxc, :) = p_global_x(1);
states_max(idxCPxc, :) = p_global_x(2);
states_min(idxCPzc, :) = p_global_z(1);
states_max(idxCPzc, :) = p_global_z(2);

% 3. Add states and controls to the collocation problem
% Add states using the adapted bounds from step 2
problem.addOptimVar('states',states_min,states_max);
% Add controls
problem.addOptimVar('controls',model.controls.xmin,model.controls.xmax);
% Add variable for CP offset in y-direction in m
problem.addOptimVar('CPYOffset', -0.1, 0.1, 0);
% Set initial guess
problem.makeinitialguess('random'); % Todo: Do no shuffle in Collocation!

% 4. Add objectives to minimize effort of muscles and torques and difference to data
% Load data and extract specific sample
trackingDataInv = TrackingData.loadStruct(dataFileInv);
trackingDataInv.trimData(iSample, iSample);
trackingDataTrans = trackingDataInv.extractData('translation', {'pelvis_tx', 'pelvis_ty', 'pelvis_tz'});
trackingDataAng = trackingDataInv.extractData('angle', {'pelvis_rotation'    'pelvis_obliquity'    'pelvis_tilt' ...
    'hip_flexion_r'    'hip_adduction_r'    'hip_rotation_r'    'knee_angle_r'    'ankle_angle_r'    'subtalar_angle_r' ...
    'mtp_angle_r'    'hip_flexion_l'    'hip_adduction_l'    'hip_rotation_l'    'knee_angle_l'    'ankle_angle_l'    ...
    'subtalar_angle_l'    'mtp_angle_l'    'lumbar_extension'   'lumbar_bending'    'lumbar_rotation'    'arm_flex_r' ...
    'arm_add_r'    'arm_rot_r'    'elbow_flex_r'    'pro_sup_r'    'arm_flex_l'    'arm_add_l'    'arm_rot_l'    'elbow_flex_l'    'pro_sup_l'});
trackingDataMeas = TrackingData.loadStruct(dataFileMea);
trackingDataMeas.trimData(iSample, iSample);
trackingDataGRF = trackingDataMeas.extractData('GRF', {'GRF_y_r', 'GRF_y_l'});
% Add tracking
problem.addObjective(@trackTranslations,W.trackTrans,trackingDataTrans);
problem.addObjective(@trackAngles,W.trackAngle,trackingDataAng);
problem.addObjective(@trackGRF,W.trackGRF,trackingDataGRF);
% Add effort
problem.addObjective(@effortTermMuscles,W.effMuscles,'volumeweighted',3); 
problem.addObjective(@effortTermTorques,W.effTorques,2); 

% 5. Add constraint to ensure static equilibrium 
problem.addConstraint(@equilibriumConstraintsCPOffset,model.constraints.fmin,model.constraints.fmax); 

end
