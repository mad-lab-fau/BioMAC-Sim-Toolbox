%======================================================================
%> @file IMUTracking2D/scriptIMU2D.m
%> @brief Script to show how we can run simulations tracking inertial sensor data.
%> 
%> @example
%> @details
%> This is a script that shows how to generate simulations that track data
%> from an inertial measurement unit (IMU), using only angular velocity
%> and acceleration measurements. These scripts are based on the
%> following paper: Dorschky, E., Nitschke, M., Seifer, A. K., van den 
%> Bogert, A. J., & Eskofier, B. M. (2019). Estimation of gait kinematics
%> and kinetics from inertial sensor data using optimal control of
%> musculoskeletal models. Journal of biomechanics, 95, 109278.
%> Additional data can be found here: DOI 10.5281/zenodo.11522049. We only
%> use the normal walking and running trial for participant 1.
%>
%> @author Anne Koelewijn and Eva Dorschky
%> @date November 2024
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
dataFolder           = 'data\IMU2D';                      % Relative from the path of the repository
dataFileWalking      = 'P01_normwalking.mat';                 % Walking data at normal speed from Dorschky et al. (2019)
dataFileRunning      = 'P01_normrunning.mat';                 % Running data at normal speed from Dorschky et al. (2019)
modelFile            = 'gait2dc_par.xls';                   % Name of our default 2D model settings 
resultFolder         = 'results\IMUTracking2D\';            % Relative from the path of the repository

%% Initalization
% Get absolute file names
resultFileStanding     = [path2repo,filesep,resultFolder,filesep, mfilename,'_standing'];
resultFileWalking      = [path2repo,filesep,resultFolder,filesep, mfilename,'_walking'];
resultFileRunning      = [path2repo,filesep,resultFolder,filesep, mfilename,'_running'];
dataFileWalking        = [path2repo,filesep,dataFolder,  filesep,dataFileWalking];
dataFileRunning        = [path2repo,filesep,dataFolder,  filesep,dataFileRunning];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end


%% Standing: Simulate standing with minimal effort without tracking data for one point in time (static)
% Create an instane of our 2D model class using the default settings
model = Gait2dc(modelFile);

% Create and solve the optimization problem only when it has not yet been solved
if exist([resultFileStanding '.mat'], 'file') == 0
    problemStanding = IntroductionExamples.standing2D(model, resultFileStanding); %We use the same standing simulation for simlicity
    solver = IPOPT();

    resultStanding = solver.solve(problemStanding);

    % Save the result
    resultStanding.save(resultFileStanding);
else
    disp('Solution already exists')
    load(resultFileStanding)
    resultStanding = result;
end

%% Walking while tracking inertial sensor data

% Load tracking data struct and create a TrackingData object
trackingDataWalking = TrackingData.loadStruct(dataFileWalking);

W.eff = 100;               % Weight of effort term  in objective
W.reg = 1e-5;              % Weight of regularization term in objective
W.track = 1;               % Weight of tracking term in objective
W.speed = 0;               % Target speed must be exactly achieved
W.dur = Inf;               % No predefined duration
isSymmetric = 0;

% Create ad solve the optimization problem
if exist([resultFileWalking '.mat'], 'file') == 0
    solver = IPOPT();
    initialGuess = resultFileStanding;
    problemWalking = IMUTracking2D.move2D_IMU(model,resultFileWalking, trackingDataWalking, initialGuess,isSymmetric,W);
    resultWalking = solver.solve(problemWalking);

    % Save the result
    resultWalking.save(resultFileWalking);
else
    disp('Solution already exists')
    load(resultFileWalking)
    resultWalking = result;
end

resultWalking.report()

%% Running while tracking inertial sensor data
% Load tracking data struct and create a TrackingData object
trackingDataRunning = TrackingData.loadStruct(dataFileRunning);

W.eff = 100;            % Weight of effort term  in objective
W.reg = 1e-5;           % Weight of regularization term in objective
W.track = 1;            % Weight of tracking term in objective
W.speed = 0;            % Target speed must be exactly achieved
W.dur = Inf;            % No predefined duration
isSymmetric = 0;

% Create ad solve the optimization problem
if exist([resultFileRunning '.mat'], 'file') == 0
    solver = IPOPT();
    initialGuess = resultFileWalking;
    problemRunning = IMUTracking2D.move2D_IMU(model,resultFileRunning, trackingDataRunning, initialGuess,isSymmetric,W);
    resultRunning = solver.solve(problemRunning);

    % Save the result
    resultRunning.save(resultFileRunning);
else
    disp('Solution already exists')
    load(resultFileRunning)
    resultRunning = result;
end

resultRunning.report()