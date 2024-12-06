%======================================================================
%> @file IntroductionExamples/script2D.m
%> @brief Script to start working with the Bio-Sim-Toolbox for 2D models
%> @example
%> @details
%> This is a script which you can use if you work with our code for the
%> first time. 
%>
%> @author Marlies Nitschke
%> @date November, 2018
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
dataFolder     = 'data/IntroductionExamples';    % Relative from the path of the repository
dataFile       = 'RBDS001runT35.mat';            % Running data from Fukuchi 2017 subject S001 with 3.5 m/s
modelFile      = 'gait2dc_par.xls';              % Name of our default 2D model settings
resultFolder   = 'results/IntroductionExamples'; % Relative from the path of the repository


%% Initalization
% Get date
dateString = datestr(date, 'yyyy_mm_dd');

% Get absolute file names
resultFileStanding = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_standing'];
resultFileRunning  = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_running'];
dataFile           = [path2repo,filesep,dataFolder,  filesep,dataFile];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end


%% Standing: Simulate standing with minimal effort without tracking data for one point in time (static)
% Create an instane of our 2D model class using the default settings
model = Gait2dc(modelFile);
% => We use a mex function for some functionality of the model. This was
% automatically initialized with the correct settings for the current
% model. (see command line output)

% Call IntroductionExamples.standing2D() to specify the optimizaton problem
% => Take a look at the function ;)
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

%% Running: Simulate running with minimal effort while tracking data
% Load tracking data struct and create a TrackingData object
trackingData = TrackingData.loadStruct(dataFile);

% We set the targetspeed of the simulation. By constraining the speed, we
% ensure that the simulation will use this speed. The speed we choose
% reflects the speed of the trackingData, but speed can also be left free. 
targetSpeed = 3.5; % m/s

% Create and automatically initalize an instance of our 2D model class.  
% To fit the tracking data we have to scale the default model. This is done
% using the height and mass of the participant.
model = Gait2dc(modelFile, trackingData.participantHeight, trackingData.participantMass);

% The data we use in this example was recorded using a treadmill. Therefore,
% we have the change the coeffienct of the air drag in the model. 
model.drag_coefficient = 0;
% => You may have noticed that the mex file was initialized again after
% changing the air drag. This is necessary to update the model parameters
% of the mex file. 

% Call IntroductionExamples.running2D() to specify the optimizaton problem
% => Take a look inside the function ;)
isSymmetric = 1;
initialGuess = resultFileStanding;
problemRunning = IntroductionExamples.running2D(model, resultFileRunning, trackingData, targetSpeed, isSymmetric, initialGuess);

% If you want to check your derivatives, you can use the code below. When
% doing so, it is best to decrease the number of nodes to, e.g., 4, to
% speed up the calculations. If it works for 4 nodes, it should also work
% for 40 or 4000 nodes.
%     problemRunning.derivativetest();

% Create solver and change solver settings
solver = IPOPT();
solver.setOptionField('max_iter', 5000);
solver.setOptionField('tol', 0.0005);

% Solve the optimization problem and save the result. 
resultRunning = solver.solve(problemRunning);
resultRunning.save(resultFileRunning); 

% Now, we can use our report function to look at the results. We can use the
% default settings for extracting variables and plotting by calling the
% function without inputs.
% However, we change the settings now to plot the initial guess and adapt
% the figure size.
settings.plotInitialGuess = 1;
style.figureSize = [0 0 16 26];
% When also giving a filename as input, the function will automatically
% create a pdf summarizing all information and plots. Take a look at it!
% You might have to press Enter a couple of times in the Command Window
% for the saving to continue.
resultRunning.report(settings, style, resultFileRunning);

% If you do not want to plot the variables, but only extract specific
% biomechanical variables, you can instead use the extractData function.
% To obtain the muscle forces of all muscles in addition to the default
% variables, we can adapt the settings struct.
settings.muscleForce = model.muscles.Properties.RowNames;
% For symmetric simulations, we can specify whether the returned
% simVarTable should reflect a full gait cycle or only the simulated
% motion, i.e. half a gait cycle.
getFullCycle = 1;
simVarTable = resultRunning.problem.extractData(resultRunning.X, settings, [], getFullCycle);
% Open simVarTable and inspect it. You will see that is has a similar
% structure as the table in trackingData.variables.

% If you have reference data in the tracking data format, you can use the
% input argument variableTable of extractData():
% variableTable = referenceData.variables;
% This data will be plotted as "measured" data in the report function and
% saved as mean_extra and var_extra in the simVarTable.

% It is worth saving this simVarTable in a sparate file to reuse when
% performing the evaluation of your study. If you also added the reference
% data before calling the report, you might have already everything you
% need for evaluation in one file.

% If you generated multiple simulations, you can use Collocation.plotMultSimVarTables()
% to plot multiple results into one graph and Collocation.plotMeanSimVarTables()
% to plot the mean and variance of multiple results into one graph. The
% example script scriptPlotting.m shows these options.

% We can also write a movie of the movement
resultRunning.problem.writeMovie(resultRunning.X, resultRunning.filename);

% We can also compute the energy expenditure of the movement
% The following command will give the metabolic cost in J/m/kg, but there
% are other options in the function as well.
resultRunning.problem.getMetabolicCost(resultRunning.X)
