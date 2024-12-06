%======================================================================
%> @file IntroductionExamples/scriptPlotting.m
%> @brief Script to show how to work with the toolbox's plotting options
%> @example
%> @details
%> This is a script which you can use if you work with our code for the
%> first time. 
%> The script first generates a couple of 2D simulations which are then
%> visualized later. 3D simulations can be visualized in a similar way
%> (might need adaption of input options of some functions).
%>
%> @author Marlies Nitschke
%> @date April, 2021
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
dataFiles      = {'RBDS001runT35.mat', 'RBDS002runT35.mat', 'RBDS003runT35.mat'}; % Running data from Fukuchi 2017 subject S001, S002, S003 with 3.5 m/s
modelFile      = 'gait2dc_par.xls';              % Name of our default 2D model settings
resultFolder   = 'results/IntroductionExamples'; % Relative from the path of the repository

%% Initalization
% Get date
dateString = datestr(date, 'yyyy_mm_dd');
curmfilename = mfilename;

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end

%% Create standing as initial guess 
%  (see script2D.m for detailed comments)
% Ceate model
model = Gait2dc(modelFile); 

% Create problem
resultFileStanding = [path2repo,filesep,resultFolder,filesep,dateString,'_', curmfilename,'_standing'];
problem = IntroductionExamples.standing2D(model, resultFileStanding);

% Create solver
solver = IPOPT();
solver.setOptionField('tol', 0.0000001);
solver.setOptionField('constr_viol_tol', 0.000001);

% Solve the optimization problem and save the result
result = solver.solve(problem);
result.save(resultFileStanding);

%% Create running simulations for different tracking data files 
%  (see script2D.m for detailed comments)
targetSpeed = 3.5; % m/s
for iData = 1 : numel(dataFiles)
    % Load tracking data
    dataFile = [path2repo,filesep,dataFolder,  filesep,dataFiles{iData}];
    trackingData = TrackingData.loadStruct(dataFile);
    
    % Create model
    model = Gait2dc(modelFile, trackingData.subjectHeight, trackingData.subjectMass);
    model.drag_coefficient = 0;
    
    % Create problem
    resultFileRunning  = [path2repo,filesep,resultFolder,filesep,dateString,'_', curmfilename,'_running', sprintf('_S%03i', iData)];
    isSymmetric = 1;
    initialGuess = resultFileStanding;
    problem = IntroductionExamples.running2D(model, resultFileRunning, trackingData, targetSpeed, isSymmetric, initialGuess);
    
    % Create solver 
    solver = IPOPT();
    solver.setOptionField('max_iter', 5000);
    solver.setOptionField('tol', 0.0005);
    
    % Solve the optimization problem and save the result
    result = solver.solve(problem);
    result.save(resultFileRunning);
    
end

%% Plot tracking data
% In general, I like the following command to let all plotted figures appear 
% within one window:
set(0, 'DefaultFigureWindowStyle', 'docked');
% And switch off the warnings this produces
warning('off', 'MATLAB:Figure:SetOuterPosition');

% You can plot all entries of the trackingData object using
% plotVarTable().
plotVarTable(trackingData.variables);


%% Create reports and obtain simulated variables for running simulations
% Close previous fgure to see what is happening now
close all;

% Define what should be included into the report and which simulated
% variables you want to extract first. See the documentation of
% extractData() for details.
% You can for example not extract neural excitation by using an empty cell
settings.u = {};
% If you want to extract something what had an empty cell before, you
% simply have to define in a cell what you want to extract. E.g. by using
% the names of the right muscles which are the first 8 muscles in the
% muscle table of the 2D model.
settings.muscleForce = result.problem.model.muscles.Properties.RowNames(1:8)';

simVarTables = cell(1, numel(dataFiles));
for iData = 1 : numel(dataFiles)
    
    % Load the result object of the single simulations
    resultFileRunning  = [path2repo,filesep,resultFolder,filesep,dateString,'_', curmfilename,'_running', sprintf('_S%03i', iData)];
    load(resultFileRunning);
   
    % Make the plots for each single simulation by calling the report
    % function. If you add a file name to the input variables, a pdf report
    % is saved. The output simVarTable will be a table containing all
    % simulated variables which you defined in settings.
    % Remember: If you only want to get the simVarTable without plotting, use
    % extractData()!
    simVarTables{iData} = result.report(settings);
    
end

%% Create report of S001 with tracking data of S002 as reference
% Close previous fgure to see what is happening now
close all;

% Load the result object of the single simulations of S001
resultFileRunning  = [path2repo,filesep,resultFolder,filesep,dateString,'_', curmfilename,'_running', sprintf('_S%03i', 1)];
load(resultFileRunning);

% Load the tracking data of S002
dataFile = [path2repo,filesep,dataFolder,  filesep,dataFiles{iData}];
trackingData = TrackingData.loadStruct(dataFile);

% Resample the tracking data of S002 to the number of collocation nodes
% used in the simulation
trackingData.resampleData(result.problem.nNodes);

% Use the struct settings from the previouse paragraph and add the table
% containing the variables
settings.variableTable = trackingData.variables;

% Make the report and this time also save the pdf report by using a
% filename as second input. The data specified in settings.variableTable is
% called "measured" in the legends of the report. 
simVarTable = result.report(settings, [], resultFileRunning);
% Off course, you do not have to output the simVarTable if you just want to
% create the pdf report.

%% Compute error measures using the simVarTable
% After extracting a simVarTable with report() or extractData(), you can
% also compute error measures for all the variables in the table.
% Simply specify which error measures you want to have and the columns with
% the result are added to the simVarTable.
measureNames = {'Pearson', 'RMSE'};
simVarTable = calculateMeasures(simVarTable, measureNames);
% Have a look into simVarTable!
% The measures are computed for each combination of the columns "sim",
% "mean" and "mean_extra" if they are provided. In this example, we obtain
% the error measures
% - <measure>_sim_track comparing the columns "sim" and "mean"
%   which is here the comparison between simulation of S001 and measured data of S001.
% - <measure>_sim_extra comparing the columns "sim" and "mean_extra"
%   which is here the comparison between simulation of S001 and measured data of S002.
% - <measure>_track_extra  comparing the columns "mean" and "mean_extra"
%   which is here the comparison between measured data of S001 and S002.
%
% This is super useful for evaluation of your simulations.


%% Create plots showing all simulation results at the same time
% Close previous fgure to see what is happening now
close all;

% We can use a cell vector or cell matrix of simVarTables as input.
Collocation.plotMultSimVarTables(simVarTables);

%% Create plots showing the mean over all simulation results
% Close previous fgure to see what is happening now
close all;

% The style of used in plotMultSimVarTables() and plotMeanSimVarTable() can
% be adjusted using the struct style. We can for example use three
% instead of two columns of subfigures.
style.subFigSettings.nCol = 3;

% We can use a cell vector or cell matrix of simVarTables as input. It takes 
% the mean along the columns of the cell. The function is however currently
% assuming that tracking and reference data is identical over all tables
% within a column and is using only the one from the first row.
Collocation.plotMeanSimVarTable(simVarTables', style);

%% Custom
% It is very likely that you can not directly use the outputs of these
% functions for your thesis itself since not everything can be perfectly
% adapted. However, they are great for data analysis and using report() is 
% the easiest to obtain simulated variables for further processing. 
% Furthermore, their code offers you a good template for writing custom 
% plotting functions. 


