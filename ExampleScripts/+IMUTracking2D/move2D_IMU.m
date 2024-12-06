%======================================================================
%> @file +IMUTracking2D/move2D_IMU.m
%> @brief Function to set up optimal control problem to simulate walking or
%> running from inertial sensor data
%>
%> @author Eva Dorschky, Anne Koelewijn
%> @date November, 2024
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for 2D running
%>
%> @param   model          Gait2dc: Model which should be used for the simulation
%> @param   resultfile     String: Name of the resultfile including path
%> @param   trackingData   TrackingData: Tracking Data containing angles and GRFs data 
%> @param   initialGuess   String: Filename with path specifying the initial guess 
%> @param   isSymmetric    Bool: Specifies if we assume movement symmetry. If so, 
%>                         we simulate only one half gait cycle. This has to fit to the tracking data.
%> @param   W              Struct: objective weigths
%> @param   targetSpeed    Double: Target speed of the movement in x direction in m/s. 
%>                         This speed will be enforced.
%> @param   targetdur      Double: Target durection of the movement in s. 
%>                         This duration will be enforced.
%> @retval  problem        Collocation: Optimization problem for 2D running
% ======================================================================
function problem = move2D_IMU(model,resultfile, trackingData, initialGuess,isSymmetric,W)

%Initialize parameters
if ~isfield(W,'effort')
    W.effort = 1;
end
if ~isfield(W,'track')
    W.track = 0;
    warning('Proceeding without tracking')
end
if ~isfield(W,'reg')
    W.reg = 1;
end
if ~isfield(W,'speed')
    W.speed = 0;
end
if ~isfield(W,'dur')
    W.dur = 0;
end

% Initialize collocation problem
N = 100; %Number of collocation nodes
Euler = 'BE';
plotLog = 1; % Show progress of optimization
problem = Collocation(model,N,Euler,resultfile,plotLog);

% Create optimization variables
xmin = model.states.xmin;
xmax = model.states.xmax;

states_min = repmat(xmin,1,N+1);
states_max = repmat(xmax,1,N+1);

idxPelvisX = model.extractState('q', 'pelvis_tx');
states_min(idxPelvisX,1) = 0; % trunk position in x direction is zero in first node
states_max(idxPelvisX,1) = 0; % trunk position in x direction is zero in first node

problem.addOptimVar('states',states_min,states_max);
problem.addOptimVar('controls',repmat(model.controls.xmin,1,N+1),repmat(model.controls.xmax,1,N+1));

% add duration 
if isSymmetric
    problem.addOptimVar('dur',0.1,1);
else
    problem.addOptimVar('dur',0.2,2);
end

% enforce speed
targetspeed =  trackingData.variables.mean{strcmp(trackingData.variables.type,'speed')}; %prescribe speed
problem.addOptimVar('speed',targetspeed,targetspeed);

% Create objective terms
if W.reg ~= 0
    problem.addObjective(@regTerm,W.reg)
end

if W.effort ~= 0
    problem.addObjective(@effortTermMuscles,W.eff,'equal',2)
    
end

if W.track ~= 0
    if trackingData.nSamples~=N
        trackingData.resampleData(N);
    end

    %By removing sensors in the list below, it is also possible to create
    %simulations with a sparse sensor set.
    problem.addObjective(@trackAcc,W.track*2,trackingData.extractData('acc',{'pelvis','femur_r','tibia_r','foot_r','femur_l','tibia_l','foot_l'},{[1,0,0],[0,1,0]}))
    problem.addObjective(@trackGyro,W.track,trackingData.extractData('gyro',{'pelvis','femur_r','tibia_r','foot_r','femur_l','tibia_l','foot_l'},[0,0,1]))
end

if W.dur ~= 0 && ~isinf(W.dur)
    problem.addObjective(@trackDuration,W.dur,trackingData);
end

% Create constraints
problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,N),repmat(model.constraints.fmax,1,N))
problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),isSymmetric)

% Make initialguess
problem.makeinitialguess(initialGuess); 

end
