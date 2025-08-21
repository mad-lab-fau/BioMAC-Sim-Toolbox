%======================================================================
%> @file ISB_tutorial/scriptTutorial_final.m
%> @brief Script with all code as used in the ISB tutorial 2025:
%> “In the wild” movement analysis using dynamic simulations
%> 
%> @example
%> @details
%> This script is used for the tutorial '“In the wild” movement analysis
%> using dynamic simulations'. To make the code work, you have to add to
%> the lines that are highlighted with TODO.
%>
%> @author Anne Koelewijn
%> @date July 2025
%======================================================================

clear all 
close all
clc

%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));
% Path to your repository
% path2repo = % TODO 1. Add the path to the location of the BioMAC-Sim-Toolbox here with single quotes. You can use the comment "what('BioMAC-Sim-Toolbox/').path;"

% Fixed settings    
dataFolder           = '';                      % Relative from the path of the repository
dataFileWalking      = 'IMU_walking_P03.mat';   % Walking data at normal speed from Dorschky et al. (2019)
modelFile            = 'gait2dc_par.xls';       % Name of our default 2D model settings 
resultFolder         = 'results\ISBTutorial\';  % Relative from the path of the repository

%% Initalization
% Get absolute file names. We normally use the code as shown in line 35,
% which uses the file name automatically. This doesn't work if you run
% sections independently.
% resultFileStanding     = [path2repo,filesep,resultFolder,filesep, mfilename,'_standing'];
resultFileStanding     = [path2repo,filesep,resultFolder,filesep, 'scriptTutorial_standing'];
resultFileWalking_base = [path2repo,filesep,resultFolder,filesep, 'scriptTutorial_walking'];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end

% Create an instane of our 2D model class using the default settings
model = Gait2dc(modelFile);

%% Standing: Simulate standing with minimal effort without tracking data for one point in time (static)
% Create and solve the optimization problem only when it has not yet been solved
if exist([resultFileStanding '.mat'], 'file') == 0
    problemStanding = ISB_Tutorial.standing2D(model, resultFileStanding); %We use the same standing simulation for simlicity
    solver = IPOPT();
    solver.setOptionField('max_iter', 20000);
    solver.setOptionField('tol', 0.0001);

    resultStanding = solver.solve(problemStanding);

    % To plot the result we have to extract the states x from the result vector X
    x = resultStanding.X(resultStanding.problem.idx.states);

    % Now, we can plot the stick figure visualizing the result
    figure();
    resultStanding.problem.model.showStick(x);
    title('2D Standing');

    % If the model is standing on the toes, the optimization ends in a local optimum and not the global one. Rerun this
    % section and you should find a different solution, due to a different random
    % initial guess. You can run it a couple of times until you find a good
    % solution, standing on flat feet. If it does not. Stop the simulation
    % using the button above or by typing "dbquit" in the Command Window
    % and press play again. If the simulation looks good, press Continue or
    % continue the code using F5.
    keyboard

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

%First create "base" simulation
resultFileWalking = resultFileWalking_base;
if exist([resultFileWalking '.mat'], 'file') == 0 %this way, we avoid overwriting by accident and save computational time and effort. Make sure to change resultFileWalking
    solver = IPOPT();
    initialGuess = resultFileStanding;
    problemWalking = ISB_Tutorial.move2D_IMU(model,resultFileWalking, trackingDataWalking, initialGuess);
    resultWalking = solver.solve(problemWalking);

    % Save the result
    resultWalking.save(resultFileWalking);
    resultWalking.report();
end

%% Objective weighting. These are the results called scriptTutorial_walking_Effortxxxx where xxxx represents the weighting used
efforts = [1 50 100 600 5000 1e4];
for iEff = 1:length(efforts)
    trackingDataWalking = TrackingData.loadStruct(dataFileWalking);
    % resultFileWalking = %TODO 5. create a unique name here
    W.track  = 1;                 % Weight of tracking term in objective
    % W.eff    = %TODO 6. set the effort weight % Weight of effort term  in objective
    W.reg    = 1e-3;              % Weight of regularization term in objective
    if exist([resultFileWalking '.mat'], 'file') == 0 %this way, we avoid overwriting by accident and save computational time and effort. Make sure to change resultFileWalking
        solver = IPOPT();
        initialGuess = resultFileStanding;
        % problemWalking = % TODO 7. create the problem by calling move2D_IMU with the correct inputs
        resultWalking = solver.solve(problemWalking);

        % Save the result
        resultWalking.save(resultFileWalking);
    end
end

%We will create a single plot with the different solutions
simVarTables_weight = cell(1, length(efforts));
for iEff = 1:length(efforts)
    % resultFileWalking =  % TODO 8. Copy from line 101xxx here %unique name
    load(resultFileWalking);

    % Extract the joint angles, joint moments, GRFs, and muscle activations into a table.
    settings.translation      = {};
    settings.angle            = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
    settings.moment           = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
    settings.GRF              = {'GRF_x_r', 'GRF_y_r'};
    settings.u                = {};
    settings.a                = {'Iliopsoas_r', 'Glutei_r', 'Hamstrings_r', 'Rectus_r', 'Vasti_r', 'Gastroc_r', 'Soleus_r', 'TibialisAnt_r'};
    % we plot the accelerometer and gyroscope signals for sensors placed at the same location as the measured data 
    settings.acc              = trackingDataWalking.variables(and(strcmp(trackingDataWalking.variables.type,'acc'),or(trackingDataWalking.variables.direction(:,1) == 1,trackingDataWalking.variables.direction(:,2) == 1)),["type", 'name', 'unit', 'segment', 'position', 'direction']); 
    settings.gyro             = trackingDataWalking.variables(and(strcmp(trackingDataWalking.variables.type,'gyro'),trackingDataWalking.variables.direction(:,3) == 1),["type", 'name', 'unit', 'segment', 'position', 'direction']);
    simVarTables_weight{iEff} = result.problem.extractData(result.X, settings);   
end

close all;
Collocation.plotMultSimVarTables(simVarTables_weight);

%% Initial guess. These results are called scriptTutorial_walkingmidIniGuess and scriptTutorial_walkingwalkIniGuess
resultFileWalking = [resultFileWalking_base 'walkIniGuess'];
trackingDataWalking = TrackingData.loadStruct(dataFileWalking);
if exist([resultFileWalking '.mat'], 'file') == 0 %this way, we avoid overwriting by accident and save computational time and effort. Make sure to change resultFileWalking
    solver = IPOPT();
    % initialGuess = % TODO 9. add the path to the results with effort 5000 here so that it will be used as initial guess
    problemWalking = ISB_Tutorial.move2D_IMU(model,resultFileWalking, trackingDataWalking, initialGuess);
    resultWalking = solver.solve(problemWalking);

    % Save the result
    resultWalking.save(resultFileWalking);
end

resultFileWalking = [resultFileWalking_base 'midIniGuess'];
trackingDataWalking = TrackingData.loadStruct(dataFileWalking);
if exist([resultFileWalking '.mat'], 'file') == 0 %this way, we avoid overwriting by accident and save computational time and effort. Make sure to change resultFileWalking
    solver = IPOPT();
    % initialGuess = % TODO 10. add the correct string to use the midpoint initial guess
    problemWalking = ISB_Tutorial.move2D_IMU(model,resultFileWalking, trackingDataWalking, initialGuess);
    resultWalking = solver.solve(problemWalking);

    % Save the result
    resultWalking.save(resultFileWalking);
end

%We will create a single plot with the different solutions, using default
%settings when using extractData
load(resultFileWalking_base);
simVarTables_IniGuess{1} = result.problem.extractData(result.X);
resultFileWalking = [resultFileWalking_base 'midIniGuess'];
load(resultFileWalking);
simVarTables_IniGuess{2} = result.problem.extractData(result.X);
resultFileWalking = [resultFileWalking_base 'walkIniGuess'];
simVarTables_IniGuess{3} = result.problem.extractData(result.X);

close all;
Collocation.plotMultSimVarTables(simVarTables_IniGuess);

%% Number of nodes. These results are called scriptTutorial_walking_nNodesXXX
nodes = [5 25 50];% cut for speed 75 100 250 500];
for iNode = 1:length(nodes)
    trackingDataWalking = TrackingData.loadStruct(dataFileWalking);

    resultFileWalking = [resultFileWalking_base '_nNodes' num2str(nodes(iNode))]; %unique name
    if exist([resultFileWalking '.mat'], 'file') == 0 %this way, we avoid overwriting by accident and save computational time and effort. Make sure to change resultFileWalking
        solver = IPOPT();
        solver.setOptionField('max_iter', 10000)
        if iNode < 2
            initialGuess = resultFileStanding;
        else %speed up simulations
             initialGuess = [resultFileWalking_base '_nNodes' num2str(nodes(iNode-1))]; %unique name
        end
        % problemWalking = TODO 11. create the problem by calling move2D_IMU with the correct inputs and number of nodes
        resultWalking = solver.solve(problemWalking);

        % Save the result
        resultWalking.save(resultFileWalking);
    end
end

%We will create a single plot with the different solutions
simVarTables_nodes = cell(1, length(nodes));
clear settings
for iNode = 1:length(nodes)
    resultFileWalking = [resultFileWalking_base '_nNodes' num2str(nodes(iNode))]; %unique name
    load(resultFileWalking);

    % Extract the joint angles and moments into a table.
    settings.translation   = {};
    settings.angle         = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
    settings.moment        = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
    settings.GRF           = {};
    settings.u             = {};
    settings.a             = settings.u;
    simVarTables_nodes{iNode} = result.problem.extractData(result.X, settings);
    metCost(iNode) = result.problem.getMetabolicCost(result.X);

end

close all;

% TODO 12: create a plot showing the metabolic cost as a function of the number of nodes
Collocation.plotMultSimVarTables(simVarTables_nodes);

%% Optional: create you own constraint function

% Load tracking data struct and create a TrackingData object
trackingDataWalking = TrackingData.loadStruct(dataFileWalking);

% Test the derivatives
solver = IPOPT();
initialGuess = resultFileStanding;
% problemWalkingDerTest = % TODO optional 8: create a problem with a small number of nodes to run the derivative test efficiently
% TODO optional 9. Run the derivative test function here. problemWalking.derivativetest() 

% Reload tracking data struct and create a TrackingData object
trackingDataWalking = TrackingData.loadStruct(dataFileWalking);

% Recreate "base" simulation with the speed constraint
resultFileWalking = [resultFileWalking_base '_speedConstraint'];
if exist([resultFileWalking '.mat'], 'file') == 0 %this way, we avoid overwriting by accident and save computational time and effort. Make sure to change resultFileWalking
    solver = IPOPT();
    initialGuess = resultFileStanding;
    
    % TODO optional 10. Create the problem using the new function name. problemWalking = 
    resultWalking = solver.solve(problemWalking);

    % Save the result
    resultWalking.save(resultFileWalking);
    resultWalking.report();
end

%We will create a single plot with the original and new solution, using default settings when using extractData 
load(resultFileWalking_base);
simVarTables_Constraint{1} = result.problem.extractData(result.X);
load(resultFileWalking);
simVarTables_Constraint{2} = result.problem.extractData(result.X);

close all;
Collocation.plotMultSimVarTables(simVarTables_IniGuess);