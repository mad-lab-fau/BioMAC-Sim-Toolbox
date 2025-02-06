%======================================================================
%> @file MarkerTracking3D/scriptMarkerTracking.m
%> @brief Example script for marker tracking (in 3D)
%> @example
%> @details
%> This is a script that shows how to generate marker tracking simulations
%> in 3D without any task constraints. These scripts are based on the
%> following paper: Nitschke, M., Marzilger, R., Leyendecker, S., Eskofier,
%> B. M., & Koelewijn, A. D. (2023). Change the direction: 3D optimal 
%> control simulation by directly tracking marker and ground reaction 
%> force data. PeerJ, 11, e14852.
%> Additional data can be found here: DOI 10.5281/zenodo.6949011. We only
%> use one trial of each condition for participant 2.
%>
%> @author Marlies Nitschke, Anne Koelewijn
%> @date September, 2024
%======================================================================

clear all 
close all
clc

%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));

% Path to your repository
path2repo = [filePath filesep '..' filesep '..'];
participant = 'Participant_02';

% Trial names
trial_Npose = 'trial0002';
trial_straightrunning = 'trial0025';
trial_curvedslowrunning = 'trial0100';
trial_vcut = 'trial0122';
movements = {'straightrunning', 'curvedrunning', 'vcut'};
trials = {trial_straightrunning, trial_curvedslowrunning, trial_vcut};

% Error measure for analysis of error between tracking and simulated data
measureNames = {'RMSE'};

% Dock all figures in one Window
set(0, 'DefaultFigureWindowStyle', 'docked'); 

% Adapt style (You might have to adapt this to best fit your screen)
style.subFigSettings.nCol = 12;
style.subFigSettings.width = 3;
style.subFigSettings.height = 2;
style.trackColor = 'k--';
style.extraColor = 'r:';
style.xLabelText  = 'Motion in \%';
style.lineWidth = 1.5;

%% Run standing simulations
% Standing with marker tracking
resultFileStandingMarker = MarkerTracking3D.run_standing_marker(path2repo, participant, trial_Npose);

% Standing for joint angle tracking
resultFileStandingAngle = MarkerTracking3D.run_standing_angle(path2repo, participant, trial_Npose);


%% Compare standing simulations (can be easily rerun to inspect the results)
% Load the marker tracking result
load(resultFileStandingMarker);

% Plot the stick figure visualizing the result
figure();
x = result.X(result.problem.idx.states);
result.problem.model.showStick(x);
markerTable = result.problem.objectiveTerms(1).varargin{1}.variables;
markerMean = cell2mat(markerTable.mean')/1000; % extract and convert from mm to meter
result.problem.model.showMarker(x, markerTable, markerMean);
title('Standing Marker Tracking');
view(45, 20);
legend(findobj(gca, 'Type', 'Scatter'), {'Measured', 'Simulated'});

% Extract simulated data
settings.marker = markerTable;
simVarTableStandingMarker = result.problem.extractData(result.X, settings);

% Compute errors to tracked data
simVarTableStandingMarker = calculateMeasures(simVarTableStandingMarker, measureNames);

% Load the angle tracking result
load(resultFileStandingAngle);

% Plot the stick figure visualizing the result
figure();
x = result.X(result.problem.idx.states);
result.problem.model.showStick(x);
result.problem.model.showMarker(x, markerTable, markerMean);
title('Standing Angle Tracking');
view(45, 20);
legend(findobj(gca, 'Type', 'Scatter'), {'Measured', 'Simulated'});

% Extract simulated data
simVarTableStandingAngle = result.problem.extractData(result.X, settings);

% Compute errors to tracked data
simVarTableStandingAngle = calculateMeasures(simVarTableStandingAngle, measureNames);


%% Run movement simulations
% Simulate the different movements as specified before. This will require a
% couple of hours per simulation. 
for iMov = 1:length(movements)
    % Marker tracking
    resultFilesRunningMarker{iMov} = MarkerTracking3D.run_motion_marker(path2repo, participant, movements{iMov}, trials{iMov});

    % Joint angle tracking
    resultFilesRunningAngle{iMov} = MarkerTracking3D.run_motion_angle(path2repo, participant, movements{iMov}, trials{iMov});
end



%% Compare movement simulations (can be easily rerun to inspect the results)
% Get simulated variables
for iMov = 1:length(movements)

    % Load the marker tracking result
    load(resultFilesRunningMarker{iMov});
    
    % Extract simulated data
    markerTable = result.problem.objectiveTerms(1).varargin{1}.variables;
    settings.marker = markerTable;
    simVarTableRunningMarker{iMov} = result.problem.extractData(result.X, settings);
    
    % Compute errors to tracked data
    simVarTableRunningMarker{iMov} = calculateMeasures(simVarTableRunningMarker{iMov}, measureNames);
    
    
    % Load the angle tracking result
    load(resultFilesRunningAngle{iMov});
    
    % Extract simulated data
    simVarTableRunningAngle{iMov} = result.problem.extractData(result.X, settings);
    
    % Compute errors to tracked data
    simVarTableRunningAngle{iMov} = calculateMeasures(simVarTableRunningAngle{iMov}, measureNames);
    
    % Make plots
    Collocation.plotMultSimVarTables({simVarTableRunningMarker{iMov}, simVarTableRunningAngle{iMov}}, style);
    legend({'Marker', 'Angle'})
    
    % Pause loop to look at results
    if iMov < length(movements)
        fprintf('Inspect the results for %s. Press a key to close the figures and continue with the next movement. \n', movements{iMov})
        pause;
        close all;
    end
end



