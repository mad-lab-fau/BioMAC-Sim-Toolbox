%======================================================================
%> @file ExoPaper/ExoPaper.m
%> @brief Script to show how we can run simulations with exoskeleton.
%> @example
%> @details
%> This script was used for Koelewijn, A. D., & Selinger, J. C. (2022). 
%> Predictive simulations to replicate human gait adaptations and 
%> energetics with exoskeletons. IEEE Transactions on Neural Systems and
%> Rehabilitation Engineering, 30, 1931-1940., to play with and
%> test some settings. It does not exactly reproduce the simulations used
%> there, since here, we only use a single virtual participant.
%>
%> @author Anne Koelewijn
%> @date August, 2022
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
dataFolder       = ['data' filesep 'Walking'];                     			% Relative from the path of the repository
dataFile         = 'Winter_normal.mat';                 			  		% Data from Winter's book (2009)
modelFile        = 'gait2dc_par.xls';                  				 		% Name of our default 2D model settings 
resultFolder     = ['results' filesep 'ExoWalking' filesep 'NoExo'];   		% Relative from the path of the repository
resultFolderHigh = ['results' filesep 'ExoWalking' filesep 'Penalize_High'];% Relative from the path of the repository --> penalize high
resultFolderLow  = ['results' filesep 'ExoWalking' filesep 'Penalize_Low']; % Relative from the path of the repository --> penalize low


%% Initalization
% Get absolute file names
resultFileStanding     = [path2repo,filesep,resultFolder,filesep, mfilename,'_standing'];
resultFileWalking      = [path2repo,filesep,resultFolder,filesep, mfilename,'_walking'];
resultFileWalkingHigh  = [path2repo,filesep,resultFolderHigh,filesep, mfilename,'_walkingexo'];
resultFileWalkingLow   = [path2repo,filesep,resultFolderLow,filesep, mfilename,'_walkingexo'];
dataFile               = [path2repo,filesep,dataFolder,  filesep,dataFile];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolderHigh], 'dir')
    mkdir([path2repo,filesep,resultFolderHigh]);
end

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolderLow], 'dir')
    mkdir([path2repo,filesep,resultFolderLow]);
end


%% Standing: Simulate standing with minimal effort without tracking data for one point in time (static)
% Create an instane of our 2D model class using the default settings
clear model

%Below is an example of how a random virtual participant could be created.
useRandomParticipant = 0; % Change to 1 to use a random virtual participant
if useRandomParticipant
	weight = 65.3+randn(1, 1)*9.8;
	BMI = 23.3+randn(1,1)*2.3;
	body_length = sqrt(weight./BMI);
	SDmusParamMatrix = 10/100*randn(14,16);
	model = Gait2dc(modelFile, body_length, weight);
	model = model.getSameRandomParticipant(SDmusParamMatrix);
else
	model = Gait2dc(modelFile);
end

% Create and solve the optimization problem only when it has not yet been solved
if exist([resultFileStanding '.mat'], 'file') == 0
    problemStanding = ExoPaper.standing2D(model, resultFileStanding);
    solver = IPOPT();

    resultStanding = solver.solve(problemStanding);

    % Save the result
    resultStanding.save(resultFileStanding);
else
    disp('Solution already exists')
    load(resultFileStanding)
    resultStanding = result;
end

%% Walking: Simulate walking with minimal effort while tracking data, first with step rate left free and optimized

% Load tracking data struct and create a TrackingData object
trackingData = TrackingData.loadStruct(dataFile);
targetSpeed = trackingData.extractData('speed', 'speed').variables.mean{1};

Wtrack = 0.1;
isSymmetric = 1;

% Create and solve the optimization problem
if exist([resultFileWalking '.mat'], 'file') == 0
    solver = IPOPT();
    initialGuess = resultFileStanding;
    problemWalking = ExoPaper.walking2D_fixedsteprate(model, Wtrack, resultFileWalking, trackingData, targetSpeed, [], isSymmetric, initialGuess);
    resultWalking = solver.solve(problemWalking);

    % Save the result
    resultWalking.save(resultFileWalking);
else
    disp('Solution already exists')
    load(resultFileWalking)
    resultWalking = result;
end

% Save the optimal steprate to create simulations with fixed step rate
idur = resultWalking.problem.idx.dur;
sr_norm = 60/(resultWalking.X(idur)*2);

% Now create simulations between 85% and 115% of the optimal step rate. We
% define the vector as follows to be able to use the previous % as initial
% guess for the new simulation (e.g. 99% to create 98%)
steprates_rel = [1.0 0.9 0.8 1.0 1.1 1.2];%[1.0:-0.01:0.85 1.0:0.01:1.15];
steprates_min = 0.85*sr_norm;
steprates_max = 1.15*sr_norm;
steprates = steprates_rel.*sr_norm;

for i = 1:length(steprates)
    % trackingData should be loaded again for each simulation
    trackingData = TrackingData.loadStruct(dataFile);
    targetSpeed = trackingData.extractData('speed', 'speed').variables.mean{1};

    isSymmetric = 1;
    if exist('resultFileWalking_now', 'var')
        initialGuess = resultFileWalking_now;
    else
        initialGuess = resultFileWalking;
    end

    resultFileWalking_now = [resultFileWalking '_relsteprate_fromprev_' num2str(steprates_rel(i)*100)];
    
    % Solve the optimization problem
    if exist([resultFileWalking_now '.mat'], 'file') == 0
        problemWalking = ExoPaper.walking2D_fixedsteprate(model, Wtrack, resultFileWalking_now, trackingData, targetSpeed, 1/(steprates(i)/30), isSymmetric, initialGuess);
        solver = IPOPT();
        resultWalking = solver.solve(problemWalking);

        % Save the result
        resultWalking.save(resultFileWalking_now);
    else
        disp('Solution already exists')
        load(resultFileWalking_now)
        resultWalking = result;
    end
end


%% Walking with exo, penalize high
% First with step rate free and optimized
penalizeHigh = 1; %should be 1 for penalize high and 0 for penalize low.

trackingData = TrackingData.loadStruct(dataFile);
targetSpeed = trackingData.extractData('speed', 'speed').variables.mean{1};

if useRandomParticipant
	model_high = Gait2dc_Exo(sr_norm, 'dur', increase, modelFile, body_length, weight); % we use the height and weight that were defined before
	model_high = model1.getSameRandomParticipant(SDmusParamMatrix); % we use the muscle parameters that were defined before
else
	model_high = Gait2dc_Exo(sr_norm, 'dur', penalizeHigh, modelFile);
end

isSymmetric = 1;
initialGuess = resultFileStanding;
           
% Solve the optimization problem
if exist([resultFileWalkingHigh '.mat'], 'file') == 0
    problemWalkingHigh = ExoPaper.walking2D_exo(model_high, Wtrack, resultFileWalkingHigh, trackingData, targetSpeed, [], isSymmetric, initialGuess);
    solver = IPOPT();
    resultWalkingHigh = solver.solve(problemWalkingHigh);
    
    %Save the result
    resultWalkingHigh.save(resultFileWalkingHigh);
else
    disp('Solution already exists')
    load(resultFileWalkingHigh)
    resultWalkingHigh = result;
end

%Now between 85% and 115% of the steprate for normal walking
resultFileWalkingHigh_now = [resultFileWalkingHigh '_relsteprate_fromprev_' num2str(steprates_rel(i)*100)];
for i = 1:length(steprates) 
    trackingData = TrackingData.loadStruct(dataFile);
    targetSpeed = trackingData.extractData('speed', 'speed').variables.mean{1};

    % Solve the optimization problem
    if exist([resultFileWalkingHigh_now '.mat'], 'file') == 0
        isSymmetric = 1;
        if i > 1
            initialGuess = resultFileWalkingHigh_now;
        else
            initialGuess = resultFileWalking;
        end
        
        problemWalkingHigh = ExoPaper.walking2D_exo(model_high, Wtrack, resultFileWalkingHigh_now, trackingData, targetSpeed, 1/(steprates(i)/30), isSymmetric, initialGuess);
        solver = IPOPT();
        resultWalkingHigh = solver.solve(problemWalkingHigh);

        % Save the result
        resultWalkingHigh.save(resultFileWalkingHigh_now);
    else
        disp('Solution already exists')
        load(resultFileWalkingHigh_now)
        resultWalkingHigh = result;
    end
end

%% Walking with exo, penalize low
% First with step rate free and optimized

penalizeLow = 0;

trackingData = TrackingData.loadStruct(dataFile);
targetSpeed = trackingData.extractData('speed', 'speed').variables.mean{1};

if useRandomParticipant
	model_high = Gait2dc_Exo(sr_norm, 'dur', penalizeLow, modelFile, body_length, weight); % we use the height and weight that were defined before
	model_high = model1.getSameRandomParticipant(SDmusParamMatrix); % we use the muscle parameters that were defined before
else
	model_low = Gait2dc_Exo(sr_norm, 'dur', penalizeLow, modelFile);
end
           
% Solve the optimization problem
if exist([resultFileWalkingLow '.mat'], 'file') == 0
    isSymmetric = 1;
    initialGuess = resultFileStanding;
    problemWalkingLow = ExoPaper.walking2D_exo(model_low, Wtrack, resultFileWalkingLow, trackingData, targetSpeed, [], isSymmetric, initialGuess);
    solver = IPOPT();
    resultWalkingLow = solver.solve(problemWalkingLow);
    
    %Save the result
    resultWalkingLow.save(resultFileWalkingLow);
else
    disp('Solution already exists')
    load(resultFileWalkingHigh)
    resultWalkingLow = result;
end

resultFileWalkingLow_now = [resultFileWalkingLow '_relsteprate_fromprev_' num2str(steprates_rel(i)*100)];
for i = 1:length(steprates)
    trackingData = TrackingData.loadStruct(dataFile);
    targetSpeed = trackingData.extractData('speed', 'speed').variables.mean{1};

    % Solve the optimization problem
    if exist([resultFileWalkingLow_now '.mat'], 'file') == 0
        isSymmetric = 1;
        if i > 1
            initialGuess = resultFileWalkingLow_now;
        else
            initialGuess = resultFileWalking;
        end
        
        problemWalkingLow = ExoPaper.walking2D_exo(model_low, Wtrack, resultFileWalkingLow_now, trackingData, targetSpeed, 1/(steprates(i)/30), isSymmetric, initialGuess);
    
        % Create solver and change solver settings
        solver = IPOPT();
        resultWalkingHigh = solver.solve(problemWalkingLow);

        % Save the result
        resultWalkingHigh.save(resultFileWalkingLow_now);
    else
        disp('Solution already exists')
        load(resultFileWalkingLow_now)
        resultWalkingHigh = result;
    end
end