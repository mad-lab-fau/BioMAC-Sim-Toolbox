%======================================================================
%> @file +InverseKinematics/setup_standing_IK.m
%> @brief Function to setup inverse kinematics of standing from marker data
%>
%> @author Anne Koelewijn
%> @date September, 2025
%======================================================================

% ======================================================================
%> @brief Function to setup standing problem with marker tracking
%>
%> @param   model       Gait3d: Model used for simulation
%> @param   dataFile    String: Filename containg data struct for tracking
%> @param   iSample     Int: Index of sample used from measured data
%> @param   resultFile  String: Filename to log files
%> @retval  problem     Collocation: Generated standing problem
% ======================================================================
function problem = setup_standing_IK(model, dataFile, iSample, resultFile)

% 1. Create problem
problem = Collocation(model,1,'BE',resultFile); % simulate only one collocation node for standing, 'BE' or 'ME' doesn't matter

% 2. Get bounds and initial guess for optimization variables
% Define bounds
angInds = model.extractState('q');      %extract indices of degrees of freedom in model
states_min = model.states.xmin(angInds);
states_max = model.states.xmax(angInds);

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
% Add variable for CP offset in y-direction in m
problem.addOptimVar('CPYOffset', -0.1, 0.1, 0);
% Set initial guess
problem.makeinitialguess('random'); % Todo: Do no shuffle in Collocation!

% 4. Add objectives to minimize effort of muscles and torques and difference to data
% Load data and extract specific sample
trackingData = TrackingData.loadStruct(dataFile);
trackingData.trimData(iSample, iSample);
trackingDataMar = trackingData.extractData('marker');
% Add tracking
problem.addObjective(@trackMarker,1,trackingDataMar);

end
