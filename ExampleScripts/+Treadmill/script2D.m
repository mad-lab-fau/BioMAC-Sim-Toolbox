%======================================================================
%> @file Treadmill/script2D.m
%> @brief Script to show how to simulate a treadmill and with the 2D 
%> OpenSim model
%> @example
%> @details
%> This is a script that you can use to see how the treadmill works. 
%>
%> @author Anne Koelewijn
%> @date August, 2024
%======================================================================

clear all 
close all
clc

%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));
% Path to your repository
path2repo = [filePath filesep '..' filesep '..' filesep];

% Fixed settings
dataFolder     = 'data/Walking';                % Relative from the path of the repository
dataFile       = 'Winter_normal.mat';           % Running data from Fukuchi 2017 subject S001 with 3.5 m/s
modelFile      = 'gait2d.osim';                 % Name of the OpenSim model with lumbar joint locked
% modelFile      = 'gait10dof18musc.osim';      % Name of the base model from OpenSim 
resultFolder   = 'results/Treadmill'; 	        % Relative from the path of the repository

%% Initalization
% Get date
dateString = datestr(date, 'yyyy_mm_dd');

% Get absolute file names
resultFileStanding = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_standing'];
resultFileWalking  = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_walking'];
dataFile           = [path2repo,filesep,dataFolder,  filesep,dataFile];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end

%% Standing: Simulate standing with minimal effort without tracking data for one point in time (static)
% Create an instane of the OpenSim 2D model class using the default settings
model = Gait2d_osim(modelFile);
% model = Gait2dc(modelFile);

% Call IntroductionExamples.standing2D() to specify the optimizaton problem
% We use the same function as in IntroductionExamples, since we are solving
% the same problem.
problemStanding = IntroductionExamples.standing2D(model, resultFileStanding);

% Create an object of class solver. We use most of the time the IPOPT here.
solver = IPOPT();

% Change settings of the solver
solver.setOptionField('tol', 0.0000001);
solver.setOptionField('constr_viol_tol', 0.000001);

% Solve the optimization problem
resultStanding = solver.solve(problemStanding);

% Save the result
resultStanding.save(resultFileStanding);

% To plot the result we have to extract the states x from the result vector X
x = resultStanding.X(resultStanding.problem.idx.states);

% Now, we can plot the stick figure visualizing the result
figure();
resultStanding.problem.model.showStick(x);
title('2D Standing');

% If the model is standing on the toes, the optimization ends in a local optimum and not the global one. Rerun this
% section and you should find a different solution, due to a different random
% initial guess. You can run it a couple of times until you find a good
% solution, standing on flat feet.

%% Simulate walking on a single belt or split-belt treadmill
% Load tracking data struct and create a TrackingData object
trackingData = TrackingData.loadStruct(dataFile);

% The global speed is going to be zeros, as the movement now comes from the treadmill
targetSpeed = 0; % m/s

% Create and initialize an instance of the OpenSim 2D model class.
model = Gait2d_osim(modelFile);
singlespeed = 1; %Change to 0 for split-belt treadmill simulation
if singlespeed
    model.setTreadmillSpeed(1.2);
else
    speed.left = 1.15;
    speed.right = 1.25;
    model.setTreadmillSpeed(speed);
end

isSymmetric = 0;
initialGuess = resultFileStanding;
problemWalking = Treadmill.walking2D(model, resultFileWalking, trackingData, targetSpeed, isSymmetric, initialGuess);

% Create solver and change solver settings
solver = IPOPT();
solver.setOptionField('max_iter', 5000);
solver.setOptionField('tol', 0.0005);

% Solve the optimization problem and save the result. 
resultWalking = solver.solve(problemWalking);
resultWalking.save(resultFileWalking); 

% If you want to create plots, take a look at one of the other examples.
resultWalking.problem.writeMovie(resultWalking.X, resultWalking.filename);