%======================================================================
%> @file IntroductionExamples/script3D.m
%> @brief Script to perform dynamic optimization of muscle states with the
%> BioMAC toolbox.
%> @example
%> @details
%> This is a script that shows how dynamic optimization of muscle states
%> can be done with the Toolbox. This approach was used in the 2019 paper:
%> Koelewijn, A. D., Heinrich, D., & van den Bogert, A. J. (2019).
%> Metabolic cost calculations of gait using musculoskeletal energy models,
%> a comparison study. PloS one, 14(9), e0222037. It assumes that the
%> joint angles and moments are somehow known already, and calculates the
%> muscle states and inputs required to achieve the moments given the joint
%> angles.
%>
%> @author Anne Koelewijn
%> @date September, 2025
%======================================================================

clear all 
close all
clc

%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));
% Path to your repository
inds = strfind(filePath,'BioMAC-Sim-Toolbox');
path2repo =  filePath(1:inds+17);

% Fixed settings
dataFolder     = 'data/MarkerTracking';               % Relative from the path of the repository
dataFile       = 'Participant_02/straightrunning/trial0025_DataStructInverse.mat';   % Straight running data from marker tracking paper
dataFileCurved = 'Participant_02/curvedrunning/trial0100_DataStructInverse.mat';   % Curved running data from marker tracking paper
modelFile      = 'Participant_02/Participant_02.osim';             % Name of our default model. You can add a new model in the osim_files folder of Gait3d, but you can also specify the path to a model similar to how it is done for the dataFile
resultFolder   = 'results/MuscleDynamicOptimization'; % Relative from the path of the repository


%% Initalization
% Get date
dateString = datestr(date, 'yyyy_mm_dd');

% Get absolute file names
resultFileRunning       = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_running'];
resultFileCurvedRunning = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_curvedRunning'];
dataFile                = [path2repo,filesep,dataFolder,  filesep,dataFile];
dataFileCurved          = [path2repo,filesep,dataFolder,  filesep,dataFileCurved];
modelFile               = [path2repo,filesep,dataFolder,  filesep,modelFile];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end


%% Straight Running: Find muscle parameters for periodic movement
trackingData = TrackingData.loadStruct(dataFile);

model = Gait3d(modelFile);

initialGuess = 'mid';
problemRunning = MuscleDynamicOptimization.walkingMuscleStates(model, resultFileRunning, trackingData, initialGuess);
% problemRunning.derivativetest();

% Create solver and change solver settings
solver = IPOPT();
solver.setOptionField('max_iter', 5000);
solver.setOptionField('tol', 0.0005);

% Solve the optimization problem and save the result. 
resultRunning = solver.solve(problemRunning);
resultRunning.save(resultFileRunning); 

%Plotting only relevant things
settings.translation   = {};
settings.angle         = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r', 'hip_flexion_l','knee_angle_l', 'ankle_angle_l'};
settings.moment        = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r', 'hip_flexion_l','knee_angle_l', 'ankle_angle_l'};
settings.GRF           = {};
settings.u             = {'Iliopsoas_r', 'Glutei_r', 'Hamstrings_r', 'Rectus_r', 'Vasti_r', 'Gastroc_r', 'Soleus_r', 'TibialisAnt_r', ...
    'Iliopsoas_l', 'Glutei_l', 'Hamstrings_l', 'Rectus_l', 'Vasti_l', 'Gastroc_l', 'Soleus_l', 'TibialisAnt_l'};
settings.a             = settings.u;
settings.s             = settings.u;
settings.LMTU          = settings.u;
settings.LCE           = settings.u;
settings.LSEE          = settings.u;
settings.muscleForce   = settings.u;
settings.CEForce       = settings.u;
settings.muscleMetRate = settings.u;
settings.duration      = 1;
settings.plotStick     = 0;
table = resultRunning.report(settings);

%% Curved Running: Find muscle parameters for non-periodic movement
trackingData = TrackingData.loadStruct(dataFile_curved);

model = Gait3d(modelFile);

initialGuess = 'mid';
problemCurvedRunning = MuscleDynamicOptimization.movementMuscleStates(model, resultFileRunning, trackingData, initialGuess);
% problemRunning.derivativetest();

% Create solver and change solver settings
solver = IPOPT();
solver.setOptionField('max_iter', 5000);
solver.setOptionField('tol', 0.0005);

% Solve the optimization problem and save the result. 
resultCurvedRunning = solver.solve(problemCurvedRunning);
resultCurvedRunning.save(resultFileCurvedRunning); 

%Plotting only relevant things
settings.translation   = {};
settings.angle         = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r', 'hip_flexion_l','knee_angle_l', 'ankle_angle_l'};
settings.moment        = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r', 'hip_flexion_l','knee_angle_l', 'ankle_angle_l'};
settings.GRF           = {};
settings.u             = {'Iliopsoas_r', 'Glutei_r', 'Hamstrings_r', 'Rectus_r', 'Vasti_r', 'Gastroc_r', 'Soleus_r', 'TibialisAnt_r', ...
    'Iliopsoas_l', 'Glutei_l', 'Hamstrings_l', 'Rectus_l', 'Vasti_l', 'Gastroc_l', 'Soleus_l', 'TibialisAnt_l'};
settings.a             = settings.u;
settings.s             = settings.u;
settings.LMTU          = settings.u;
settings.LCE           = settings.u;
settings.LSEE          = settings.u;
settings.muscleForce   = settings.u;
settings.CEForce       = settings.u;
settings.muscleMetRate = settings.u;
settings.duration      = 1;
settings.plotStick     = 0;
table = resultCurvedRunning.report(settings);