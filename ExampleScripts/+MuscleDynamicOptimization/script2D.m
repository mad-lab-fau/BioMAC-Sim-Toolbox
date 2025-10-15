%======================================================================
%> @file MuscleDynamicOptimization/script2D.m
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
path2repo = [filePath filesep '..' filesep '..' filesep];

% Fixed settings
dataFolder     = 'data/MuscleDynamicOptimization';    % Relative from the path of the repository
dataFile       = 'P01_normwalking.mat';            % Running data from Fukuchi 2017 subject S001 with 3.5 m/s
modelFile      = 'gait2dc_par.xls';              % Name of our default 2D model settings
resultFolder   = 'results/MuscleDynamicOptimization'; % Relative from the path of the repository


%% Initalization
% Get date
dateString = datestr(date, 'yyyy_mm_dd');

% Get absolute file names
resultFileRunning  = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_running'];
dataFile           = [path2repo,filesep,dataFolder,  filesep,dataFile];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end

%% Running: Find muscle parameters for running
trackingData = TrackingData.loadStruct(dataFile);

model = Gait2dc(modelFile, trackingData.participantHeight, trackingData.participantMass);
model.drag_coefficient = 0;

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