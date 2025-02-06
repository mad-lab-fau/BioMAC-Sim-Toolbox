%======================================================================
%> @file +MarkerTracking3D/run_motion_marker.m
%> @brief Function to do simulation with marker tracking
%>
%> @author Marlies Nitschke
%> @date June, 2021
%======================================================================

%======================================================================
%> @brief Function to do standing simulation with marker tracking
%>
%> @details
%> This function tracks data to simulate a motion without periodicity.
%>
%> @param  workDirectory    String: Current work directory
%> @param  subject          String: Name of subject (e.g. Subject_02)
%> @param  motion           String: Name of motion (e.g. straightrunning)
%> @param  trial            String: Name of triak (e.g. trial0025)
%> @param  WMar             String: Weight of marker tracking term in objective
%> @param  WGRF             String: Weight of GRF tracking term in objective
%> @param  WReg             String: Weight of the regularization term in objective
%> @retval resultFile       String: Filename of result
%======================================================================
function resultFile = run_motion_marker(workDirectory, subject, motion, trial)

% Fixed settings
dataFolder     = ['data' filesep 'MarkerTracking']; % Relative from the work directory
dataTrackFile  = [subject filesep motion filesep trial '_DataStructMeasured.mat'];
modelFile      = [subject filesep subject '.osim'];
resultFolder   = ['results' filesep 'MarkerTracking3D' filesep subject];  % Relative from the path of the repository
resultFileStanding = 'standing_marker';
resultFile     = sprintf('%s_%s_motion_marker', motion, trial);

% Get absolute file names
resultFileStanding = [workDirectory filesep resultFolder filesep resultFileStanding];
resultFile         = [workDirectory filesep resultFolder filesep resultFile];
dataTrackFile      = [workDirectory filesep dataFolder filesep dataTrackFile];
modelFile          = [workDirectory filesep dataFolder filesep modelFile];

% Create resultfolder if it does not exist
if ~exist([workDirectory,filesep,resultFolder], 'dir')
    mkdir([workDirectory,filesep,resultFolder]);
end

% Create an instane of our 3D model class
model = Gait3d(modelFile);

% Adjust CPs based on standing solution
resultStanding = load(resultFileStanding);
model.CPs.position = model.CPs.position + repmat([0 resultStanding.result.X(resultStanding.result.problem.idx.CPYOffset) 0], model.nCPs, 1);

% Specify the optimizaton problem
W.trackMarker= 1e-2;         % Weight of marker tracking term in objective
W.trackGRF   = 1e-3;         % Weight of GRF tracking term in objective
W.effMuscles = 1e+00;        % Weight of effort term for muscles in objective
W.effTorques = 1e-01;        % Weight of effort term for torques in objective
W.reg        = 1e-3;         % Weight of regularization term in objective
initialGuess = resultFileStanding;
nNodesEarlier = 10; %Number of time points before the actual motion of interest to ensure motion artefacts are before the motion of interest
problemRunning = MarkerTracking3D.setup_motion_marker(model,dataTrackFile,initialGuess,resultFile,W,nNodesEarlier);

% Solve
solver = IPOPT();
solver.setOptionField('max_iter', 20000);
solver.setOptionField('tol', 0.0001);
resultRunning = solver.solve(problemRunning);
resultRunning.save(resultFile);

end
