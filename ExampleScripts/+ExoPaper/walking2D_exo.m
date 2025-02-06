%======================================================================
%> @file walking2D_fixedsteprate.m
%> @brief Function to specify the optimization problem for 2D running
%>
%> @author Anne Koelewijn
%> @date October, 2019
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for normal 2D walking
%> 
%> @details
%> This function creates walking simulation with an exoskeleton. If targetDuration
%> is empty, the duration is left free and optimized. If it is given, the targetDuration
%> is added as constraint
%>
%> @param   model          Gait2dc: Model which should be used for the simulation
%> @param   Wtracking      Double: Weight used for the tracking objective
%> @param   resultfile     String: Name of the resultfile including path
%> @param   trackingData   TrackingData: Tracking Data containing angles and GRFs data 
%> @param   targetSpeed    Double: Target speed of the movement in x direction in m/s. 
%>                         This target speed will be enforced.
%> @param   targetDuration Double: Target duration of the movement in steps/min. 
%>                         This duration will be enforced. Duration will be
%>                         optimized if it is empty: [].
%> @param   isSymmetric    Bool: Specifys weather we assume symmetry of the
%>                         movement. If we assume symmetry, we simulate only one half 
%>                         of gait cycle. This has to fit to the tracking data.
%> @param   initialGuess   String: Filename with path specifying the initial guess 
%> @retval  problem        Collocation: Optimization problem for 2D running
% ======================================================================
function problem = walking2D_exo(model, Wtracking, resultfile, trackingData, targetSpeed, targetDuration, isSymmetric, initialGuess)

%% Fixed settings
nNodes = 40;   
Euler = 'BE';
logfile = resultfile;
plotLog = 1;

%% Create collocation problem
problem = Collocation_Exo(model, nNodes, Euler, logfile, plotLog);

%% Add states and controls including their bounds and initial values to the problem
% Get upper and lower bounds of the model and resize it
xmin = repmat(model.states.xmin, 1, nNodes+1); 
xmax = repmat(model.states.xmax, 1, nNodes+1);

% Adapt the bounds for the first node to start the movement at X = 0
xmax(model.extractState('q', 'pelvis_tx'), 1) = 0;
xmin(model.extractState('q', 'pelvis_tx'), 1) = 0; 

problem.addOptimVar('states', xmin, xmax);
problem.addOptimVar('controls',repmat(model.controls.xmin,1,nNodes+1), repmat(model.controls.xmax,1,nNodes+1));
problem.addOptimVar('dur',0.2, 2);
problem.addOptimVar('speed',0, 5); 

problem.makeinitialguess(initialGuess); 

%% Add objective terms to the problem
if isSymmetric
    trackingData.resampleData(nNodes*2);
    trackingData.useHalfGaitCycleData();
else
    trackingData.resampleData(nNodes);
end

GRFSignals   = {'GRF_x_r', 'GRF_y_r', 'GRF_x_l', 'GRF_y_l'};
AngleSignals = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'hip_flexion_l', 'knee_angle_l', 'ankle_angle_l'};

nGRF   = length(GRFSignals);
nAngle = length(AngleSignals);
nAll   = nGRF + nAngle;
problem.addObjective(@trackGRF   , Wtracking*nGRF/nAll  , trackingData.extractData('GRF', GRFSignals));
problem.addObjective(@trackAngles, Wtracking*nAngle/nAll, trackingData.extractData('angle', AngleSignals));

Weffort = 100;
weightsType = 'equal'; 
exponent = 3; 
problem.addObjective(@effortTermMuscles, Weffort, weightsType, exponent); 

Wreg = 0.0001;
problem.addObjective(@regTerm, Wreg);


%% Add constraints to the problem
problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,nNodes),repmat(model.constraints.fmax,1,nNodes))
problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),isSymmetric)
problem.addConstraint(@speedConstraint, 0,0, targetSpeed)

if ~isempty(targetDuration)
    problem.addConstraint(@durationConstraint, 0,0, targetDuration);
end
end