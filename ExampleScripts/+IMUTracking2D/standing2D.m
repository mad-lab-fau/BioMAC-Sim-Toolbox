%======================================================================
%> @file ExoPaper/standing2D.m
%> @brief Function to specify the optimization problem for 2D standing
%>
%> @author Marlies Nitschke
%> @date November, 2018
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for 2D standing
%>
%> @param   model          Gait2dc: Model which should be used for the simulation
%> @param   resultfile     String: Name of the resultfile including path
%> @retval  problem        Collocation: Optimization problem for 2D standing
% ======================================================================
function problem = standing2D(model, resultfile)

%% Fixed settings
nNodes = 1;   
Euler = 'BE';
logfile = resultfile;
plotLog = 1;


%% Create collocation problem
problem = Collocation(model, nNodes, Euler, logfile, plotLog);

%% Add states and controls including their bounds and initial values to the problem
xmin = model.states.xmin;
xmax = model.states.xmax;

xmin(model.extractState('q', 'pelvis_tilt'))        = -30 /180*pi; % some forward or backward lean may be needed
xmax(model.extractState('q', 'pelvis_tilt'))        =  30 /180*pi;
xmin(model.extractState('q', 'pelvis_tx'))          =   0        ; % should be zero due to translation invariance in horizontal plane
xmax(model.extractState('q', 'pelvis_tx'))          =   0        ;
xmin(model.extractState('q', 'pelvis_ty'))          =   0.5      ; % pelvis should be about 1 m above ground during standing
xmax(model.extractState('q', 'pelvis_ty'))          =   1.5      ;
xmin(model.extractState('q', {'hip_flexion_r', 'hip_flexion_l'}))  = -30 /180*pi; % some hip flexion may be needed
xmax(model.extractState('q', {'hip_flexion_r', 'hip_flexion_l'}))  =  30 /180*pi;
xmin(model.extractState('q', {'knee_angle_r' , 'knee_angle_l'}))   = -30 /180*pi; % some knee flexion may be needed
xmax(model.extractState('q', {'knee_angle_r' , 'knee_angle_l'}))   =  10 /180*pi;
xmin(model.extractState('q', {'ankle_angle_r', 'ankle_angle_l'}))  = -30 /180*pi; % some ankle flexion may be needed
xmax(model.extractState('q', {'ankle_angle_r', 'ankle_angle_l'}))  =  30 /180*pi;

% Get a random initial guess for the states (within [xmin, xmax])
xinit = xmin + (xmax-xmin) .* rand(size(xmin));

problem.addOptimVar('states', xmin, xmax, xinit);
problem.addOptimVar('controls', model.controls.xmin, model.controls.xmax, model.controls.xneutral); 

%% Add objective term(s) to the problem
Weff = 1; 
weightsType = 'equal'; 
exponent = 2; 
problem.addObjective(@effortTermMuscles, Weff, weightsType, exponent); 

%% Add constraint(s) to the problem
problem.addConstraint(@equilibriumConstraints, model.constraints.fmin, model.constraints.fmax);

end