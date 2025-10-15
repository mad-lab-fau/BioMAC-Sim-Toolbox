%======================================================================
%> @file +MarkerTracking3D/setup_motion_IK.m
%> @brief Function to setup motion problem with marker tracking
%>
%> @author Anne Koelewijn
%> @date September, 2025
%======================================================================

% ======================================================================
%> @brief Function to setup motion problem with marker tracking
%>
%> @param   model          Gait3d: Model used for simulation
%> @param   dataFile       String: Filename containg data struct for tracking
%> @param   initialFile    String: Result used for initial guess
%> @param   resultFile     String: Filename to log files
%> @param   W              Struct: Weights for objective terms
%> @param   nNodesEarlier  Double: Number of samples which we start earlier
%> @retval  problem        Collocation: Generated standing problem
% ======================================================================
function problem = setup_motion_IK(model,dataFile,initialFile,resultFile)

% Load tracking data
trackingData = TrackingData.loadStruct(dataFile);
indicesStartEnd = trackingData.movementEvents.index(strcmp(trackingData.movementEvents.name, 'R_IC'));
trackingData.trimData(indicesStartEnd(1), indicesStartEnd(2)-1);
N = trackingData.nSamples;

% Create problem
Euler = 'BE';
plotLog = 0;
problem = Collocation(model,N,Euler,resultFile,plotLog);

% Add variables which are optimized
angInds = model.extractState('q');      %extract indices of degrees of freedom in model
states_min = repmat(model.states.xmin(angInds),1,N);
states_max = repmat(model.states.xmax(angInds),1,N);
idxCPxc = model.extractState('xc');
idxCPzc = model.extractState('zc');
p_global_x = [-5, 5];
p_global_z = [-5, 5];
states_min(idxCPxc, :) = p_global_x(1);
states_max(idxCPxc, :) = p_global_x(2);
states_min(idxCPzc, :) = p_global_z(1);
states_max(idxCPzc, :) = p_global_z(2);
problem.addOptimVar('states',states_min,states_max);

% Initialize the problem with an old result specified in initialFile
problem.makeinitialguess(initialFile); 

% Add tracking terms
trackingDataMar = trackingData.extractData('marker');
problem.addObjective(@trackMarker,1,trackingDataMar);

end
