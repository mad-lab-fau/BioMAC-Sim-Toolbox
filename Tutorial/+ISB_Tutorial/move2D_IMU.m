%======================================================================
%> @file +ISB_tutorial/move2D_IMU.m
%> @brief Function to set up optimal control problem to simulate walking or
%> running from inertial sensor data
%>
%> @author Eva Dorschky, Anne Koelewijn
%> @date November, 2024
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for 2D running, used
%> in the 2025 ISB tutorial
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
function problem = move2D_IMU(model,resultfile, trackingData, initialGuess, varargin)

if isempty(varargin)
    W.track  = 1;                 % Weight of tracking term in objective
    W.eff    = 600;               % Weight of effort term  in objective
    W.reg    = 1e-3;              % Weight of regularization term in objective
    N = 100;
else
    for i = 1:length(varargin)
        if isa(varargin{i}, 'struct')
            W = varargin{i};
        else
            W.track  = 1;                 % Weight of tracking term in objective
            W.eff    = 600;               % Weight of effort term  in objective
            W.reg    = 1e-3;              % Weight of regularization term in objective
        end
        
        if isa(varargin{i}, 'double')
            N = varargin{i};
        else
            N = 100;
        end
    end
end

Euler = 'BE'; % we use backward Euler for discretization
isSymmetric = 0; % we assume no symmetry
plotLog = 1; % Show progress of optimization

problem = Collocation(model,N,Euler,resultfile,plotLog);

% Tracking data preprocessing
if trackingData.nSamples~=N
    trackingData.resampleData(N);
end
targetspeed =  trackingData.variables.mean{strcmp(trackingData.variables.type,'speed')};

% Create optimization variables
xmin = model.states.xmin;
xmax = model.states.xmax;

states_min = repmat(xmin,1,N+1); %We add N+1 states because we apply the periodicity constraint on the (N+1)th node, which should be the same as the first
states_max = repmat(xmax,1,N+1);

idxPelvisX = model.extractState('q', 'pelvis_tx');
states_min(idxPelvisX,1) = 0; % trunk position in x direction is zero in first node to avoid an infinite amount of possible solutions
states_max(idxPelvisX,1) = 0; % trunk position in x direction is zero in first node to avoid an infinite amount of possible solutions

problem.addOptimVar('states',states_min,states_max);
problem.addOptimVar('controls',repmat(model.controls.xmin,1,N+1),repmat(model.controls.xmax,1,N+1));

% add gait cycle duration as optimization variable
problem.addOptimVar('dur',0.2,2);

% prescribe speed
% problem.addOptimVar('speed', 0, 10); %TODO optional 2: change the optimization variable definition for the speed to allow a range of speeds
problem.addOptimVar('speed', targetspeed, targetspeed);

% Create objective terms
problem.addObjective(@regTerm,W.reg)
problem.addObjective(@effortTermMusclesAct,W.eff,'equal',2)
    
%By removing sensors in the list below, it is also possible to create simulations with a sparse sensor set.
problem.addObjective(@trackAcc,W.track*2,trackingData.extractData('acc',{'pelvis','femur_r','tibia_r','foot_r','femur_l','tibia_l','foot_l'},{[1,0,0],[0,1,0]}))
problem.addObjective(@trackGyro,W.track,trackingData.extractData('gyro',{'pelvis','femur_r','tibia_r','foot_r','femur_l','tibia_l','foot_l'},[0,0,1]))

% Create constraints
problem.addConstraint(@dynamicConstraints,repmat(model.constraints.fmin,1,N),repmat(model.constraints.fmax,1,N))
problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),isSymmetric)
% TODO optional 7. add a line of code to add the speed constraint

% Make initialguess
problem.makeinitialguess(initialGuess); 

end
