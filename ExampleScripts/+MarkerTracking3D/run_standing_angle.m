%======================================================================
%> @file +MarkerTracking3D/run_standing_angle.m
%> @brief Function to do standing simulation with angle tracking
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================

%======================================================================
%> @brief Function to do standing simulation with angle tracking
%>
%> @details
%> This function tracks data during standing to optimize the CP position.
%>
%> @param  workDirectory    String: Current work directory
%> @param  subject          String: Name of subject (e.g. Subject_02)
%> @param  trial            String: Name of trial (e.g. trial0002)
%> @param  WTra             String: Weight of translation tracking term in objective
%> @param  WAng             String: Weight of angle tracking term in objective
%> @param  WGRF             String: Weight of GRF tracking term in objective
%> @retval resultFile       String: Filename of result
%======================================================================
function resultFile = run_standing_angle(workDirectory, subject, trial)

% Fixed settings
dataFolder     = ['data' filesep 'MarkerTracking']; % Relative from the work directory
dataMeasFile   = [subject filesep 'N-Pose' filesep trial '_DataStructMeasured.mat'];
dataInvFile    = [subject filesep 'N-Pose' filesep trial '_DataStructInverse.mat'];
modelFile      = [subject filesep subject '.osim'];
resultFolder   = ['results' filesep 'MarkerTracking3D' filesep subject];  % Relative from the path of the repository

% Get absolute file names
resultFile   = [workDirectory filesep resultFolder filesep 'standing_angle'];
dataMeasFile = [workDirectory filesep dataFolder filesep dataMeasFile];
dataInvFile  = [workDirectory filesep dataFolder filesep dataInvFile];
modelFile    = [workDirectory filesep dataFolder filesep modelFile];

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
    iSample = 10; % Index of the used sample in the data
    W.trackTrans = 1e-3;        % Weight of translation tracking term in objective
    W.trackAngle = 1e-1;        % Weight of angle tracking term in objective
    W.trackGRF   = 1e-2;        % Weight of GRF tracking term in objective
    W.effMuscles = 1e+00;       % Weight of effort term for muscles in objective
    W.effTorques = 1e-01;       % Weight of effort term for torques in objective
    problemStanding = MarkerTracking3D.setup_standing_angle(model, dataInvFile, dataMeasFile, iSample, resultFile, W);
    
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
