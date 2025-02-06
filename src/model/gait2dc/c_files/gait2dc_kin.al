% gait2dc_kin.al
%
% Generate forward kinematics for the 2D gait model
% (is used in gait2dc_dyn.al, gait2dc_acc.al, gait2dc_FK.al, gait2dc_stick.al)
%
% Author: Ton van den Bogert
% Last revised: 01/25/2025

PAUSE 0
AUTOZ ON
OVERWRITE ALL

% Define model parameters as variables, so they can be modified from Matlab
Constants par__ThighLen, par__ShankLen

%--------------------------------------------------------------------
%	Ground reference frame, units
%--------------------------------------------------------------------
Newtonian   Ground
UnitSystem  kg,meter,sec

%---------------------------------------------------
% Generalized coordinates
%---------------------------------------------------
MotionVariables' q{9}''
% matrices needed for Jacobians
q = [q1, q2, q3, q4, q5, q6, q7, q8, q9]
qd = [q1', q2', q3', q4', q5', q6', q7', q8', q9']
qdd = [q1'', q2'', q3'', q4'', q5'', q6'', q7'', q8'', q9'']

%----------------------------------
% Define bodies and points
%----------------------------------
Bodies Trunk, RThigh, RShank, RFoot, LThigh, LShank, LFoot
Points Hip, RKnee, LKnee, RAnkle, LAnkle

% The next section of code creates, for each segment, the position, velocity, and acceleration
% of the joint connecting it to the parent, as well as the angular position, velocity, and
% acceleration of the segment.

%-----------------------------------
% Trunk
%-----------------------------------
% from ground, translation (x,y)=(q1,q2) and rotation q3
P_GroundO_Hip> = Vector(Ground, q1, q2, 0)		% trunk translation: x and y of hip
V_Hip_Ground> = Vector(Ground, q1', q2', 0)		% and velocity
A_Hip_Ground> = Vector(Ground, q1'', q2'', 0)	% and acceleration
Simprot(Ground, Trunk, 3, q3);					% trunk rotation about z axis
W_Trunk_Ground> = q3' * Ground3>					
ALF_Trunk_Ground> = q3'' * Ground3>					

%-----------------------------------
% RThigh
%-----------------------------------
% from trunk, rotation q4 on z axis
Simprot(Trunk, RThigh, 3, q4)
W_RThigh_Ground> = (q4' + q3') * Ground3>	
ALF_RThigh_Ground> = (q4'' + q3'') * Ground3>	

%-----------------------------------
% RShank
%-----------------------------------
% from thigh: rotation q5 on z axis
Simprot(RThigh, RShank, 3, q5)
P_Hip_RKnee> = -par__ThighLen*RThigh2>;     % define the knee position
V2pts(Ground, RThigh, Hip, RKnee);          % its velocity is needed also
A2pts(Ground, RThigh, Hip, RKnee);          % and acceleration
W_RShank_Ground> = (q5' + q4' + q3') * Ground3>		
ALF_RShank_Ground> = (q5'' + q4'' + q3'') * Ground3>		

%-----------------------------------
% RFoot
%-----------------------------------
% from shank, rotation q6 on z axis
Simprot(RShank, RFoot, 3, q6)
P_RKnee_RAnkle> = -par__ShankLen*RShank2>; % define the ankle position
V2pts(Ground, RShank, RKnee, RAnkle);      % its velocity is needed also
A2pts(Ground, RShank, RKnee, RAnkle);      % and acceleration
W_RFoot_Ground> = (q6' + q5' + q4' + q3') * Ground3>         % angular velocity
ALF_RFoot_Ground> = (q6'' + q5'' + q4'' + q3'') * Ground3>   % angular velocity

%-----------------------------------
% LThigh
%-----------------------------------
% from trunk, rotation q7 on z axis
Simprot(Trunk, LThigh, 3, q7)
W_LThigh_Ground> = (q7' + q3') * Ground3>	
ALF_LThigh_Ground> = (q7'' + q3'') * Ground3>		

%-------------------
% LShank
%-------------------
% from thigh, rotation q8 on z axis
Simprot(LThigh, LShank, 3, q8)
P_Hip_LKnee> = -par__ThighLen*LThigh2>;  % define the knee position
V2pts(Ground, LThigh, Hip, LKnee);       % its velocity is needed also
A2pts(Ground, LThigh, Hip, LKnee);       % and acceleration
W_LShank_Ground> = (q8' + q7' + q3') * Ground3>		
ALF_LShank_Ground> = (q8'' + q7'' + q3'') * Ground3>		

%-----------------------------------
% LFoot
%-----------------------------------
% from shank, rotation q6 on z axis
Simprot(LShank, Lfoot, 3, q9)
P_LKnee_LAnkle> = -par__ShankLen*LShank2>;  % define the ankle position
V2pts(Ground, LShank, LKnee, LAnkle);       % its velocity is needed also
A2pts(Ground, LShank, LKnee, LAnkle);       % and acceleration
W_LFoot_Ground> = (q9' + q8' + q7' + q3') * Ground3>	% angular velocity
ALF_LFoot_Ground> = (q9'' + q8'' + q7'' + q3'') * Ground3>	% angular velocity
