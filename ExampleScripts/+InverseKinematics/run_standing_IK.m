%======================================================================
%> @file +InverseKinematics/run_standing_IK.m
%> @brief Function to do inverse kinematics of standing motion with marker
%> tracking objective
%>
%> @author Anne Koelewijn
%> @date September, 2025
%======================================================================

%======================================================================
%> @brief Function to do inverse kinematics of standing with marker tracking
%>
%> @details
%> This function tracks data during standing to optimize the CP position.
%>
%> @param  workDirectory    String: Current work directory
%> @param  subject          String: Name of subject (e.g. Participant_02)
%> @param  trial            String: Name of trial (e.g. trial0002)
%> @retval resultFile       String: Filename of result
%======================================================================
function resultFile = run_standing_IK(workDirectory, subject, trial)

% Fixed settings
dataFolder     = ['data' filesep 'MarkerTracking']; % Relative from the work directory
dataTrackFile  = [subject filesep 'N-Pose' filesep trial '_DataStructMeasured.mat'];
modelFile      = [subject filesep subject '.osim'];
resultFolder   = ['results' filesep 'InverseKinematics' filesep subject];  % Relative from the path of the repository

% Get absolute file names
resultFile     = [workDirectory filesep resultFolder filesep 'ik_standing'];
dataTrackFile  = [workDirectory filesep dataFolder filesep dataTrackFile];
modelFile      = [workDirectory filesep dataFolder filesep modelFile];

% Create resultfolder if it does not exist
if ~exist([workDirectory,filesep,resultFolder], 'dir')
    mkdir([workDirectory,filesep,resultFolder]);
end

% Create an instane of our 3D model class
model = Gait3d(modelFile);

% Solve simulation for multiple random initial guesses
rng('default');
nRep = 10;
results = cell(nRep, 1);
objSum = nan(nRep, 1);
for iRep = 1 : nRep
    % Specify the optimizaton problem
    iSample = 10;            % Index of the used sample in the data
    problemStanding = InverseKinematics.setup_standing_IK(model, dataTrackFile, iSample, resultFile);
    
    % Solve
    solver = IPOPT();
    solver.setOptionField('max_iter', 20000);
    solver.setOptionField('tol', 0.0001);
    results{iRep} = solver.solve(problemStanding);
    objSum(iRep) = sum([results{iRep}.problem.objectiveTerms.weightedValue]);
end

% Save simulation result of result with lowest objective value
[~, iMinObj] = min(objSum);
result = results{iMinObj};
result.save(resultFile);

end
