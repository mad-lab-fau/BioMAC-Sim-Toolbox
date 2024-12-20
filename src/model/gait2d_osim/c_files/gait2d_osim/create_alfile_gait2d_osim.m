function create_alfile_gait2d_osim(alfile)
% Create the Autolev source file
% This is done with custom code for this model.  Eventually we may interpret the .osim file
% so we can convert any model from Opensim to Autolev.
% NOTE: We may not continue with Autolev.  See, for example: http://pydy.org/ which could
% be a good alternative.
global nbodies

alfile_FK = strrep(alfile,'.al','_FK.al');      % commands with _FK from "create_alfile_gait3d_pelvis213"
fprintf('Making %s...\n', alfile_FK);
fid = fopen(alfile_FK,'w');
if (fid < 0)
    error('Could not write %s\n', alfile_FK);
end
nbodies = 0;
rawcfile = strrep(alfile_FK,'.al','_raw.c');

% preliminaries
fprintf(fid,'%% This file was generated by %s on %s\n', mfilename,date);
fprintf(fid,'%% Autolev source code for the gait2d_osim multibody model\n\n');
fprintf(fid,'AUTOZ ON\n');
% fprintf(fid,'OVERWRITE ALL\n');
fprintf(fid,'PAUSE 0\n');
fprintf(fid,'UnitSystem  kg,meter,sec\n');
fprintf(fid,'Newtonian ground\n');
fprintf(fid,'\n');

% --------------------------------------------
% Generalized coordinates q
%
% Our name    Name in the original .osim model
% ---------------------------------------------
% q1..........pelvis_tilt
% q2..........pelvis_tx
% q3..........pelvis_ty
% q4..........hip_flexion_r
% q5..........knee_flexion_r
% q6..........ankle_angle_r
% q7..........hip_flexion_l
% q8..........knee_flexion_l
% q9..........ankle_angle_l
% 0........lumbar_extension
ndof = 9;
% ------------------------------------
fprintf(fid,'MotionVariables'' q{%d}''''\n', ndof);   	% this generates: MotionVariables' q{10}''
fprintf(fid,'VARIABLES G{%d}\n', 6);					% ground reaction force variables G1..G6
fprintf(fid,'CONSTANTS par__bodyweight\n');				% body weight (Newtons)

% ------------------------------------
% Pelvis
% ------------------------------------
% pelvis generalized coordinates may seem unusual, but this is what the Opensim model uses, we need the same here!
addbody(fid, 'pelvis', 'ground', 	...
    'q2', 'q3',                     ...		% position
    'SIMPROT', 'q1');                       % orientation

%------------------------------------
% Legs
%------------------------------------
addleg(fid, 'R', 4);			% add a right leg, using DOFs q4..q6
addleg(fid, 'L', 7);			% add a left leg, using DOFs q7..q9

%------------------------------------
% Torso
%------------------------------------
joint = 'back';					% this is how the joint was named in osim
tx = ['par__' joint '_x'];
ty = ['par__' joint '_y'];
fprintf(fid, 'CONSTANTS %s, %s\n', tx, ty);
addbody(fid, 'torso', 'pelvis', ...
    tx, ty,                    ...     % position
    'BODY123', '0');                   % orientation: lumbar extension

% New section according to "create_alfile_gait3d_pelvis213"
%-----------------------------------------------------------------------------------------------------
% Generate symbolic Jacobians
%-----------------------------------------------------------------------------------------------------
% generate matrices so we can generate Jacobians with respect to q, qd, qdd,
% and G (ground reaction forces)
genmatrix(fid, 'q', ndof);
genmatrix(fid, 'qd', ndof);
genmatrix(fid, 'qdd', ndof);
fprintf(fid, 'dfk_dq    = ZEE(D(fk,q))\n');
fprintf(fid, 'fkdot     = ZEE(DT(fk))\n');		% IS THIS OUTPUT NEEDED?  We can do fkdot = dfk/dq * qdot
fprintf(fid, 'dfkdot_dq = ZEE(D(fkdot,q))\n');

%-----------------------------------------------------------
% Generate C code
%-----------------------------------------------------------
fprintf(fid, 'Encode fk, dfk_dq, fkdot, dfkdot_dq\n');
fprintf(fid, 'Code Algebraic() %s\n', rawcfile);

% close the file
fclose(fid);

% Open File for dynamics without derivatives
alfile_dynNoDer = strrep(alfile,'.al','_NoDer.al');
fprintf('Making %s...\n', alfile_dynNoDer);
fid = fopen(alfile_dynNoDer,'w');
if (fid < 0)
    error('Could not write %s\n', alfile_dynNoDer);
end
rawcfile = strrep(alfile_dynNoDer,'.al','_raw.c');

% preliminaries
fprintf(fid,'%% This file was generated by %s on %s\n', mfilename,date);
fprintf(fid,'%% Autolev source code for the gait2d_osim multibody model w/ lumbarlock \n\n');
fprintf(fid,'AUTOZ ON\n');
% fprintf(fid,'OVERWRITE ALL\n');
fprintf(fid,'PAUSE 0\n');
fprintf(fid,'\n');

%--------------------------------------------------------------------
% Call FK
%--------------------------------------------------------------------
fprintf(fid,'RUN %s\n', alfile_FK);
fprintf(fid,'\n');

%--------------------------------------------------------------------
% Apply gravity
%--------------------------------------------------------------------
fprintf(fid, '%%------------------------------------------------------\n');
fprintf(fid, 'CONSTANTS par__gravity_x, par__gravity_y\n');
fprintf(fid, 'Gravity(Vector(ground, par__gravity_x, par__gravity_y, 0))\n');

%--------------------------------------------------------------------
% Apply air drag force to the Trunk CM
%--------------------------------------------------------------------
fprintf(fid, '%%------------------------------------------------------\n');
fprintf(fid, 'CONSTANTS par__wind, par__airdrag\n');
fprintf(fid, 'airspeed> = V_pelvisO_Ground> - par__wind * ground1>\n');		% velocity of trunk CM, relative to the air
fprintf(fid, 'VARIABLES sx, sy, s\n');                                      % xy components of airspeed, and magnitude
fprintf(fid, 'sx = dot(airspeed>,ground1>)\n');
fprintf(fid, 'sy = dot(airspeed>,ground2>)\n');
% define magnitude of airspeed such that it can't be zero
fprintf(fid, 's = sqrt(sx*sx + sy*sy + 1e-6)\n');
fprintf(fid, 'FORCE(groundO/pelvisO, -par__airdrag*s*airspeed>)\n');		% quadratic equation

%--------------------------------------------------------------------
% Apply ground reaction force/moment to each of the foot segments
%--------------------------------------------------------------------
fprintf(fid, '%%------------------------------------------------------\n');
addgrf(fid, 'Rfoot',1);         % 1,2,3 for Rfoot
addgrf(fid, 'Lfoot',4);         % 4,5,6 for Lfoot
ngrfseg = 2;					% number of segments in which GRF are applied

%--------------------------------------------------------------------
% Generate equations of motion, scaled to bodyweight
%--------------------------------------------------------------------
fprintf(fid, '%%------------------------------------------------------\n');
fprintf(fid, 'f = ( Fr() + FrStar() ) / par__bodyweight\n');

%-----------------------------------------------------------
% Generate C code
%-----------------------------------------------------------
fprintf(fid, 'Encode f\n');
fprintf(fid, 'Code Algebraic() %s\n', rawcfile);

% close the file
fclose(fid);

% Open File for dynamics with derivatives
alfile_dynDer = alfile;
fprintf('Making %s...\n', alfile_dynDer);
fid = fopen(alfile_dynDer,'w');
if (fid < 0)
    error('Could not write %s\n', alfile_dynDer);
end
rawcfile = strrep(alfile_dynDer,'.al','_raw.c');

% preliminaries
fprintf(fid,'%% This file was generated by %s on %s\n', mfilename,date);
fprintf(fid,'%% Autolev source code for the gait2d_osim multibody model\n\n');
fprintf(fid,'AUTOZ ON\n');
% fprintf(fid,'OVERWRITE ALL\n');
fprintf(fid,'PAUSE 0\n');
fprintf(fid,'\n');

%--------------------------------------------------------------------
% Call dynamics file without derivatives
%--------------------------------------------------------------------
fprintf(fid,'RUN %s\n', alfile_dynNoDer);
fprintf(fid,'\n');

%-----------------------------------------------------------------------------------------------------
% Generate symbolic Jacobians
%-----------------------------------------------------------------------------------------------------
% generate matrices so we can generate Jacobians with respect to G (ground reaction forces)
genmatrix(fid, 'G', 3*ngrfseg);
fprintf(fid, 'df_dq     = ZEE(D(f,q))\n');
fprintf(fid, 'df_dqd    = ZEE(D(f,qd))\n');
fprintf(fid, 'df_dqdd   = ZEE(D(f,qdd))\n');
fprintf(fid, 'df_dG     = ZEE(D(f,G))\n');

%-----------------------------------------------------------
% Generate C code
%-----------------------------------------------------------
fprintf(fid, 'Encode df_dq, df_dqd, df_dqdd, df_dG\n');
fprintf(fid, 'Code Algebraic() %s\n', rawcfile);
fprintf(fid, 'EXIT\n');

% close the file
fclose(fid);

fprintf('Done.\n');

end
%=======================================================================
function addbody(fid, name, parentname, t1, t2, rotationtype, r1)
% adds a body to the model

global nbodies

fprintf(fid, '%%------------------------------------------------------\n');

% if name contains the string "massless_", it's a massless reference frame
massless = ~isempty(strfind(name,'massless_'));
if (massless)
    % Define a coordinate frame for the massless body
    name = strrep(name, 'massless_', '');		% strip off the "massless_" part from the name
    fprintf(fid,'FRAME %s\n', name);
else
    % define the body and its mass properties
    fprintf(fid,'BODY %s\n', name);
    fprintf(fid,'MASS %s = par__%s_M\n', name, name);
    fprintf(fid,'INERTIA  %s, 0, 0, par__%s_Izz\n', name, name);
end

% point where body is connected to its parent, and its velocity
point = [name 'Joint'];
fprintf(fid,'POINT %s\n', point);
if strcmp(parentname, 'ground')
    fprintf(fid,'P_groundO_%s> = Vector(ground, %s, %s, 0)\n', point, t1, t2);
    fprintf(fid,'V_%s_ground> = DT(P_groundO_%s>, ground)\n', point, point);
else
    fprintf(fid,'P_%sJoint_%s> = Vector(%s, %s, %s, 0)\n', parentname, point, parentname, t1, t2);
    fprintf(fid,'V2PTS(ground, %s, %sJoint, %s)\n', parentname, parentname, point);
end

% orientation and angular velocity of the body that was just added
if strcmp(rotationtype,'SIMPROT')
    if strcmp(r1,'0')
        fprintf(fid,'SIMPROT(%s, %s, 3, %s)\n', parentname, name, r1);
        fprintf(fid,'W_%s_ground> = W_%s_ground> + Ground3> * %s\n', name, parentname, r1);
    else
        fprintf(fid,'SIMPROT(%s, %s, 3, %s)\n', parentname, name, r1);
        fprintf(fid,'W_%s_ground> = W_%s_ground> + Ground3> * %s''\n', name, parentname, r1);
    end
else
    fprintf(fid,'DIRCOS(%s, %s, %s, %s, 0, 0)\n', parentname, name, rotationtype, r1);
    fprintf(fid,'ANGVEL(%s, %s, %s, %s, 0, 0)\n', parentname, name, rotationtype, r1);
end

if (~massless)
    % center of mass position and velocity
    fprintf(fid,'CONSTANTS par__%s_CMx, par__%s_CMy\n', name, name);
    fprintf(fid,'P_%s_%sO> = Vector(%s, par__%s_CMx, par__%s_CMy, 0)\n', point,name,name,name,name);
    fprintf(fid,'V2PTS(ground, %s, %s, %sO)\n', name,point,name); 	% generate velocity of CM relative to ground
end

% Add position and orientation to the forward kinematics array FK. This will contain
% 12*Nsegments values for the whole model. FK does not include ground.
% Forward kinematics is needed for marker tracking and for contact models.
if (nbodies == 0)
    fprintf(fid,'fk = [');
else
    fprintf(fid,'fk := [fk;');
end
nbodies = nbodies + 1;
% stored as a column of 6 values, position of origin and the rotation matrix elements, extracted row-wise
fprintf(fid,'dot(P_groundO_%s>, ground1>);', point);
fprintf(fid,'dot(P_groundO_%s>, ground2>);', point);
%fprintf(fid,'dot(P_groundO_%s>, ground3>);', point);
% Lose n = 5 elements from the rotation matrix that contain a 3
fprintf(fid,'ground_%s[1,1];', name);
fprintf(fid,'ground_%s[1,2];', name);
fprintf(fid,'ground_%s[2,1];', name);
fprintf(fid,'ground_%s[2,2]]\n', name);
fprintf(fid,'\n');

% print name on the screen
% fprintf('%s\n', name);

end
%==========================================================================
function addleg(fid, side, iq)
% adds a leg to the model

% fid.......File handle for writing Autolev code
% side......(string) 'R' or 'L'
% iq........(int) index of the first of the 3 DOF in the leg

if side == 'R'
    sign = '';
elseif side == 'L'
    sign = '-';
else
    error('addleg: side must be R or L');
end

%------------------------------------
% Femur
%------------------------------------
dof1 = ['q' num2str(iq)];				% flexion
joint = [side 'hip'];
femur = [side 'femur'];
tx = ['par__' joint '_x'];
ty = ['par__' joint '_y'];
fprintf(fid, 'constants %s, %s\n', tx, ty);
addbody(fid, femur, 'pelvis', ...
    tx, ty, ...					% position
    'SIMPROT', dof1); 			% orientation

%------------------------------------
% Tibia
%------------------------------------
dof = ['q' num2str(iq+1)];			% knee angle
tibia = [side 'tibia'];
joint = [side 'knee'];

% x and y position of knee in femur are functions of knee angle (moving center of rotation)
% we use a 4th order polynomial (5 coefficients)
% parameters: joint z and polynomial coefficients for joint x and y
tx = [		  'par__' joint '_x1 * ' dof '^4' ...
    '+ par__' joint '_x2 * ' dof '^3' ...
    '+ par__' joint '_x3 * ' dof '^2' ...
    '+ par__' joint '_x4 * ' dof      ...
    '+ par__' joint '_x5'];		% expression for X position (function of dof)
ty = [		  'par__' joint '_y1 * ' dof '^4' ...
    '+ par__' joint '_y2 * ' dof '^3' ...
    '+ par__' joint '_y3 * ' dof '^2' ...
    '+ par__' joint '_y4 * ' dof      ...
    '+ par__' joint '_y5'];		% expression for Y position (function of dof)
fprintf(fid, 'constants par__%s_x{5}, par__%s_y{5}\n', joint, joint);

if side == 'R'
    % Note: Prosthesis alignment does not work with simprot
    addbody(fid, tibia, femur, ...
        tx, ty,  ...							% position
        'SIMPROT', dof);                       % orientation
else
    addbody(fid, tibia, femur, ...
        tx, ty, ...							% position
        'SIMPROT', dof);					% orientation
end

%------------------------------------
% Foot
%------------------------------------
dof = ['q' num2str(iq+2)];
joint = [side 'ankle'];
foot = [side 'foot'];
tx = '0';
ty = ['par__' joint '_y'];
fprintf(fid, 'CONSTANT %s\n', ty);			% location of ankle joint in tibia (x and z are hard zeros)
addbody(fid, foot, tibia, ...
    tx, ty, ...					% position
    'SIMPROT', dof); 			% orientation

end
%============================================================================
function genmatrix(fid, name, n)

qname = name(1);			% first letter is used as qname

% letter 2 and 3 determine if matrix elements must be q, q', or q''
if numel(name) == 1
    suffix = '';
elseif numel(name) == 2
    suffix = '''';
elseif numel(name) == 3
    suffix = '''''';
else
    error('genmatrix: name must have length 1, 2, or 3.');
end

fprintf(fid, '%s = [', name);
for i = 1:n
    if (i<n)
        fprintf(fid, '%s%d%s, ', qname, i, suffix);
    else
        fprintf(fid, '%s%d%s]\n', qname, i, suffix);
    end
end

end
%============================================================================
function addgrf(fid, segment, igrf)
% apply ground reaction force to a segment
% GRF is represented by 3 variables: Fxy, Mz in segment reference frame
%
% segment...........(string) segment where GRF should be applied
% igrf..............(integer) index to the location of the 3 GRF variables in the G array

% Define a point on segment that coincides with the ground origin
point = [segment 'GRF'];
fprintf(fid, 'POINTS %s\n', point);
fprintf(fid, 'P_groundO_%s> = 0>\n', point);
fprintf(fid, 'V2PTS(ground, %s, %sO, %s)\n', segment, segment, point);	% compute its velocity, needed for equations of motion

% Apply the ground reaction force and moment to the segment
% Remember that the forces and moments in G are normalized to body weight, so convert to N and Nm
fprintf(fid, 'FORCE(groundO/%s, par__bodyweight*VECTOR(ground, G%d, G%d, 0))\n', point, igrf, igrf+1);
fprintf(fid, 'TORQUE(ground/%s, par__bodyweight*VECTOR(ground, 0, 0, G%d))\n', segment, igrf+2);

end
