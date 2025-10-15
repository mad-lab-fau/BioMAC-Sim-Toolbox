%======================================================================
%> @file +InverseKinematics/run_motion_IK.m
%> @brief Function to do simulation with marker tracking
%>
%> @author Anne Koelewijn
%> @date September, 2025
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
function resultFile = run_motion_IK(workDirectory, subject, motion, trial)

% Fixed settings
dataFolder     = ['data' filesep 'MarkerTracking']; % Relative from the work directory
dataTrackFile  = [subject filesep motion filesep trial '_DataStructMeasured.mat'];
modelFile      = [subject filesep subject '.osim'];
resultFolder   = ['results' filesep 'InverseKinematics' filesep subject];  % Relative from the path of the repository
resultFileStanding = 'ik_standing';
resultFile     = sprintf('%s_%s_ik_motion', motion, trial);

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
initialGuess = resultFileStanding;
problemRunning = MarkerTracking3D.setup_motion_IK(model,dataTrackFile,initialGuess,resultFile);

% Solve
solver = IPOPT();
solver.setOptionField('max_iter', 20000);
solver.setOptionField('tol', 0.0001);
resultRunning = solver.solve(problemRunning);
resultRunning.save(resultFile);

end
