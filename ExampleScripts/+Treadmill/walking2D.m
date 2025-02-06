%======================================================================
%> @file TreadMillWalking/walking2D.m
%> @brief Function to specify the optimization problem for 2D running
%>
%> @author Anne Koelewijn
%> @date August 2024
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for 2D running
%>
%> @param   model          Gait2dc: Model which should be used for the simulation
%> @param   resultfile     String: Name of the resultfile including path
%> @param   trackingData   TrackingData: Tracking Data containing angles and GRFs data 
%> @param   targetSpeed    Double: Target speed of the movement in x direction in m/s. 
%>                         This target speed will be enforced.
%> @param   isSymmetric    Bool: Specifys weather we assume symmetry of the
%>                         movement. If we assume symmetry, we simulate only one half 
%>                         of gait cycle. This has to fit to the tracking data.
%> @param   initialGuess   String: Filename with path specifying the initial guess 
%> @retval  problem        Collocation: Optimization problem for 2D running
% ======================================================================
function problem = walking2D(model, resultfile, trackingData, targetSpeed, isSymmetric, initialGuess)

%% Fixed settings
% We can choose the number of collocation nodes.
nNodes = 40;   
% Most of the time we use backard euler for discretization which is encoded with 'BE'.
Euler = 'BE';
% We usually use the name of the resultfile for the name of the logfile
logfile = resultfile;
% We want to plot intermediate results during solving of the problem.
plotLog = 1;


%% Create collocation problem
problem = Collocation(model, nNodes, Euler, logfile, plotLog);


%% Add states and controls including their bounds and initial values to the problem
% Get upper and lower bounds of the model and resize it
xmin = repmat(model.states.xmin, 1, nNodes+1); %We use one more than the number of nodes since the last node (the +1) one is then used in the periodicity constraint
xmax = repmat(model.states.xmax, 1, nNodes+1);

% Adapt the bounds for the first node to start the movement at X = 0
xmax(model.extractState('q', 'pelvis_tx'), 1) = 0;
xmin(model.extractState('q', 'pelvis_tx'), 1) = 0; 

% Add states (initial values will be specified later)
problem.addOptimVar('states', xmin, xmax);

% Add controls to the problem using the default bounds (initial values will be specified later)
problem.addOptimVar('controls',repmat(model.controls.xmin,1,nNodes+1), repmat(model.controls.xmax,1,nNodes+1));

% Add duration of the movement 
problem.addOptimVar('dur',0.2, 2);

% Add speed in x direction of the movement. We choose here targetspeed for
% the lower and upper bound to ensure that we have exactly this speed. You
% could also choose other bounds.
problem.addOptimVar('speed',targetSpeed, targetSpeed);

% After adding all the components to X, we can use a previous solution as
% initial guess. In this case, we use the standing solution we produced
% before.
problem.makeinitialguess(initialGuess); 



%% Add objective terms to the problem
Wtracking = 1;
trackingData.resampleData(nNodes);
GRFSignals   = {'GRF_x_r', 'GRF_y_r', 'GRF_x_l', 'GRF_y_l'};
AngleSignals = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'hip_flexion_l', 'knee_angle_l', 'ankle_angle_l'};
nGRF   = length(GRFSignals);
nAngle = length(AngleSignals);
nAll   = nGRF + nAngle;
problem.addObjective(@trackGRF   , Wtracking*nGRF/nAll  , trackingData.extractData('GRF', GRFSignals));
problem.addObjective(@trackAngles, Wtracking*nAngle/nAll, trackingData.extractData('angle', AngleSignals));

Weffort = 1;
weightsType = 'equal'; 
exponent = 3; 
problem.addObjective(@effortTermMuscles, Weffort, weightsType, exponent); 

Wreg = 0.0001;
problem.addObjective(@regTerm, Wreg);


%% Add constraints to the problem
problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,nNodes),repmat(model.constraints.fmax,1,nNodes))
problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),isSymmetric)


end
