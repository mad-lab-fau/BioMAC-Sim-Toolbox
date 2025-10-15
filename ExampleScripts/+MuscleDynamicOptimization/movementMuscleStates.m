%======================================================================
%> @file MuscleDynamicOptimization/movementMuscleStates.m
%> @brief Function to specify the optimization problem for dynamic
%> optimization of muscle states
%>
%> @author Anne Koelewijn
%> @date September, 2025
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for dynamic
%> optimization of muscle states. The joint angles and joint moments are
%> predefined
%>
%> @param   model          Gait2dc: Model which should be used for the simulation
%> @param   resultfile     String: Name of the resultfile including path
%> @param   trackingData   TrackingData: Tracking Data containing angles and GRFs data 
%> @param   initialGuess   String: Filename with path specifying the initial guess 
%> @retval  problem        Collocation: Optimization problem for 2D running
% ======================================================================
function problem = movementMuscleStates(model, resultfile, trackingData, initialGuess)

%% 1. Fixed settings
% We can choose the number of collocation nodes.
nNodes = 40;   
% Most of the time we use backard euler for discretization which is encoded with 'BE'.
Euler = 'BE';
% We usually use the name of the resultfile for the name of the logfile
logfile = resultfile;
% We want to plot intermediate results during solving of the problem.
plotLog = 1;


%% 2. Create collocation problem
problem = Collocation(model, nNodes, Euler, logfile, plotLog);


%% 3. Add states and controls including their bounds and initial values to the problem
% Get upper and lower bounds of the muscle states and resize it for all
% time nodes
musInds = sort([model.extractState('a'); model.extractState('s')]); %extract indices of muscle states in model
xmin = repmat(model.states.xmin(musInds), 1, nNodes+1); %We use one more than the number of nodes since the last node (the +1) one is then used in the periodicity constraint
xmax = repmat(model.states.xmax(musInds), 1, nNodes+1);

% Add states (initial values will be specified later)
problem.addOptimVar('states_mus', xmin, xmax);

% Add controls to the problem using the default bounds (initial values will be specified later)
problem.addOptimVar('controls',repmat(model.controls.xmin,1,nNodes+1), repmat(model.controls.xmax,1,nNodes+1));

problem.makeinitialguess(initialGuess); 

%% 4. Add objective terms to the problem

Weffort = 1; %We only have one objective
weightsType = 'equal'; 
exponent = 3; 
problem.addObjective(@effortTermMusclesAct, Weffort, weightsType, exponent); 

%% 5. Add constraints to the problem
trackingData.resampleData(nNodes);
problem.addConstraint(@muscleDynamicConstraints,repmat(model.constraints.fmin(strcmp(model.constraints.type, 'muscles')),1,nNodes),repmat(model.constraints.fmax(strcmp(model.constraints.type, 'muscles')),1,nNodes), trackingData.extractData('angle'), trackingData.extractData('dur'),0)

% Constrain joint moments
problem.addConstraint(@jointMomentConstraint,zeros(height(model.joints)*nNodes,1),zeros(height(model.joints)*nNodes,1), trackingData.extractData('moment'), trackingData.extractData('angle'), trackingData.extractData('dur'),0)

end
