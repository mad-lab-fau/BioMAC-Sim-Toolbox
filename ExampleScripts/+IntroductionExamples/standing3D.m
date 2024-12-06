%======================================================================
%> @file IntroductionExamples/standing3D.m
%> @brief Function to specify the optimization problem for 3D standing
%>
%> @author Marlies Nitschke
%> @date February, 2019
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for 3D standing
%>
%> @param   model          Gait3d: Model which should be used for the simulation
%> @param   resultfile     String: Name of the resultfile including path
%> @retval  problem        Collocation: Optimization problem for 3D standing
% ======================================================================
function problem = standing3D(model, resultfile)

%% Fixed settings
% Number of collocation nodes is one since we want to simulate static standing 
% at one time point
nNodes = 1;   
% Discretization method is unimportant here since we do not have to compute deriatives. 
% However, most of the time we use backard euler which is encoded with 'BE'.
Euler = 'BE';
% We usually use the name of the resultfile for the name of the logfile
logfile = resultfile;
% We want to plot intermediate results during solving of the problem.
plotLog = 1;


%% Create collocation problem
problem = Collocation(model, nNodes, Euler, logfile, plotLog);


%% Add states and controls including their bounds and initial values to the problem
% Get upper and lower bounds for the states
xmin = model.states.xmin;
xmax = model.states.xmax;

% Before adding the states, we want to change the bounds of the DOFs to make
% it easier for the solver to find a proper solution.
xmin(model.extractState('q', 'pelvis_rotation'))    =   0 /180*pi; % should be zero due to rotation invariance 
xmax(model.extractState('q', 'pelvis_rotation'))    =   0 /180*pi;
xmin(model.extractState('q', 'pelvis_obliquity'))   = - 5 /180*pi; % should be zero due to symmetry of model
xmax(model.extractState('q', 'pelvis_obliquity'))   =   5 /180*pi;
xmin(model.extractState('q', 'pelvis_tilt'))        = -30 /180*pi; % some forward or backward lean may be needed
xmax(model.extractState('q', 'pelvis_tilt'))        =  30 /180*pi;
xmin(model.extractState('q', 'pelvis_tx'))          =   0        ; % should be zero due to translation invariance in horizontal plane
xmax(model.extractState('q', 'pelvis_tx'))          =   0        ;
xmin(model.extractState('q', 'pelvis_ty'))          =   0.5      ; % pelvis should be about 1 m above ground during standing
xmax(model.extractState('q', 'pelvis_ty'))          =   1.5      ;
xmin(model.extractState('q', 'pelvis_tz'))          =   0        ; % should be zero due to translation invariance in horizontal plane
xmax(model.extractState('q', 'pelvis_tz'))          =   0        ;
xmin(model.extractState('q', {'hip_flexion_r', 'hip_flexion_l'}))         = -30 /180*pi; % some flexion may be needed
xmax(model.extractState('q', {'hip_flexion_r', 'hip_flexion_l'}))         =  30 /180*pi;
xmin(model.extractState('q', {'hip_adduction_r', 'hip_adduction_l'}))     = -10 /180*pi; % some may be needed
xmax(model.extractState('q', {'hip_adduction_r', 'hip_adduction_l'}))     =  10 /180*pi;
xmin(model.extractState('q', {'hip_rotation_r', 'hip_rotation_l'}))       = -30 /180*pi; % should be small
xmax(model.extractState('q', {'hip_rotation_r', 'hip_rotation_l'}))       =  30 /180*pi;
xmin(model.extractState('q', {'knee_angle_r', 'knee_angle_l'}))           = -30 /180*pi; % some flexion may be needed
xmax(model.extractState('q', {'knee_angle_r', 'knee_angle_l'}))           =  10 /180*pi;
xmin(model.extractState('q', {'ankle_angle_r', 'ankle_angle_l'}))         = -30 /180*pi; % some may be needed
xmax(model.extractState('q', {'ankle_angle_r', 'ankle_angle_l'}))         =  30 /180*pi;
xmin(model.extractState('q', {'subtalar_angle_r', 'subtalar_angle_l'}))   = -30 /180*pi; % some may be needed
xmax(model.extractState('q', {'subtalar_angle_r', 'subtalar_angle_l'}))   =  30 /180*pi;
xmin(model.extractState('q', {'mtp_angle_r', 'mtp_angle_l'}))             = -30 /180*pi; % some may be needed
xmax(model.extractState('q', {'mtp_angle_r', 'mtp_angle_l'}))             =  30 /180*pi;
xmin(model.extractState('q', 'lumbar_extension'))   = -10 /180*pi; % some may be needed
xmax(model.extractState('q', 'lumbar_extension'))   =  10 /180*pi;
xmin(model.extractState('q', 'lumbar_bending'))     = -10 /180*pi; % should be zero due symmetry of model
xmax(model.extractState('q', 'lumbar_bending'))     =  10 /180*pi;
xmin(model.extractState('q', 'lumbar_rotation'))    = -10 /180*pi; % should be zero due symmetry of model
xmax(model.extractState('q', 'lumbar_rotation'))    =  10 /180*pi;
xmin(model.extractState('q', {'arm_flex_r', 'arm_flex_l'}))        = -88 /180*pi; % some may be needed
xmax(model.extractState('q', {'arm_flex_r', 'arm_flex_l'}))        =  88 /180*pi;
xmin(model.extractState('q', {'arm_add_r', 'arm_add_l'}))          = -88 /180*pi; % some may be needed
xmax(model.extractState('q', {'arm_add_r', 'arm_add_l'}))          =  88 /180*pi;
xmin(model.extractState('q', {'arm_rot_r', 'arm_rot_l'}))          = -88 /180*pi; % some may be needed
xmax(model.extractState('q', {'arm_rot_r', 'arm_rot_l'}))          =  88 /180*pi;
xmin(model.extractState('q', {'elbow_flex_r', 'elbow_flex_l'}))    = -10 /180*pi; % some may be needed
xmax(model.extractState('q', {'elbow_flex_r', 'elbow_flex_l'}))    =  88 /180*pi;
xmin(model.extractState('q', {'pro_sup_r', 'pro_sup_l'}))          = -10 /180*pi; % some may be needed
xmax(model.extractState('q', {'pro_sup_r', 'pro_sup_l'}))          = 150 /180*pi;

% Get a random initial guess for the states (within [xmin, xmax])
xinit = xmin + (xmax-xmin) .* rand(size(xmin));

% Add states to the problem using the bounds and inital values which were
% adapted above
problem.addOptimVar('states', xmin, xmax, xinit);

% Add controls to the problem using the default bounds and neutral controls
problem.addOptimVar('controls', model.controls.xmin, model.controls.xmax, model.controls.xneutral); 


%% Add objective term(s) to the problem
% For this simulation, we want to minimize muscle effort but we do not want
% to track data. 
% We have one term for muscle effort, for which we use a weight of 1.
Weff = 1; 
% In the effort term, we use volume dependent weighting to account for
% different muscle sizes.
weightsType = 'volumeweighted'; 
% In the effort term, we want to use the squared neural excitation of the muscles.
exponent = 2; 
% Now we can add the objective term by specifiying the function handle 'effortTermMuscles'.
problem.addObjective(@effortTermMuscles, Weff, weightsType, exponent); 

% Additionally, we have arms in our 3D model which are actuated using
% torques. We use here the same weighting and the same exponent as above.
problem.addObjective(@effortTermTorques, Weff, exponent); 

%% Add constraint(s) to the problem
% For this simulation, we are using only one constraint. The function
% 'equilibriumConstraints' computes the dynamic constraints (Netwonian's
% equations, muscle activation dynamics and contact equations) for the
% special case of nNodes = 1 (static equilibrium). We have to specify also
% the bounds which are defined by the model.
problem.addConstraint(@equilibriumConstraints, model.constraints.fmin, model.constraints.fmax); 


end