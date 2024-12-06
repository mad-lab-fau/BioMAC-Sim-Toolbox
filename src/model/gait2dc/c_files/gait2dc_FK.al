% gait2dc_FK.al
%
% Equations of motion for the 2D gait model with directional deformable contact
%
% Author: Ton van den Bogert
% Last revised: 02/07/2011

PAUSE 0
AUTOZ ON
OVERWRITE ALL

% Define model parameters as variables, so they can be modified from Matlab
Constants par__TrunkMass, par__TrunkInertia, par__TrunkCMy
Constants par__ThighMass, par__ThighInertia, par__ThighCMy, par__ThighLen
Constants par__ShankMass, par__ShankInertia, par__ShankCMy, par__ShankLen
Constants par__FootMass, par__FootInertia, par__FootCMx, par__FootCMy
Constants par__bodyweight, par__airdrag, par__wind, par__slope

%--------------------------------------------------------------------
%	Ground reference frame, units
%--------------------------------------------------------------------
Newtonian   Ground
UnitSystem  kg,meter,sec
Constants   par__gravity					% Local gravitational acceleration

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

%-----------------------------------
% Trunk
%-----------------------------------
% define kinematics
P_GroundO_Hip> = Vector(Ground, q1, q2, 0)		% trunk translation: x and y of hip
V_Hip_Ground> = Vector(Ground, q1', q2', 0)		% and velocity
A_Hip_Ground> = Vector(Ground, q1'', q2'', 0)	% and acceleration
Simprot(Ground, Trunk, 3, q3);					% trunk rotation about z axis

% define mass properties and CM
Mass     Trunk = par__TrunkMass
Inertia  Trunk, 0, 0, par__TrunkInertia
P_Hip_TrunkO> = par__TrunkCMy*Trunk2>

% define angular velocity, angular acceleration, velocity, acceleration
W_Trunk_Ground> = q3' * Ground3>					
ALF_Trunk_Ground> = q3'' * Ground3>					
V2pts(Ground, Trunk, Hip, TrunkO);		
A2pts(Ground, Trunk, Hip, TrunkO);		

%-----------------------------------
% RThigh
%-----------------------------------
% define kinematics (hinge joint, rotation q4 on z axis)
Simprot(Trunk, RThigh, 3, q4)

% define mass properties and CM
Mass     RThigh = par__ThighMass
Inertia  RThigh, 0, 0, par__ThighInertia
P_Hip_RThighO> = par__ThighCMy*RThigh2>

% define angular velocity, angular acceleration, velocity, acceleration
W_RThigh_Ground> = (q4' + q3') * Ground3>	
ALF_RThigh_Ground> = (q4'' + q3'') * Ground3>	
V2pts(Ground, RThigh, Hip, RThighO);				
A2pts(Ground, RThigh, Hip, RThighO);		

%-----------------------------------
% RShank
%-----------------------------------
% define kinematics (hinge joint, rotation q5 on z axis)
Simprot(RThigh, RShank, 3, q5)
P_Hip_RKnee> = -par__ThighLen*RThigh2>;			% define the knee position
V2pts(Ground, RThigh, Hip, RKnee);			% its velocity is needed also
A2pts(Ground, RThigh, Hip, RKnee);			% and acceleration

% define mass properties and CM
Mass     RShank = par__ShankMass
Inertia  RShank, 0, 0, par__ShankInertia	
P_RKnee_RShankO> = par__ShankCMy*RShank2>

% define angular velocity, angular acceleration, velocity, acceleration
W_RShank_Ground> = (q5' + q4' + q3') * Ground3>		
ALF_RShank_Ground> = (q5'' + q4'' + q3'') * Ground3>		
V2pts(Ground, RShank, RKnee, RShankO);					
A2pts(Ground, RShank, RKnee, RShankO);		

%-----------------------------------
% RFoot
%-----------------------------------
% define kinematics (hinge joint, rotation q6 on z axis)
Simprot(RShank, RFoot, 3, q6)
P_RKnee_RAnkle> = -par__ShankLen*RShank2>;	% define the ankle position
V2pts(Ground, RShank, RKnee, RAnkle);					% its velocity is needed also
A2pts(Ground, RShank, RKnee, RAnkle);					% and acceleration

% define mass properties and CM
Mass     RFoot = par__FootMass
Inertia  RFoot, 0, 0, par__FootInertia
P_RAnkle_RFootO> = par__FootCMx*RFoot1> + par__FootCMy*RFoot2>

% define angular velocity, angular acceleration, velocity, acceleration
W_RFoot_Ground> = (q6' + q5' + q4' + q3') * Ground3>	% angular velocity
ALF_RFoot_Ground> = (q6'' + q5'' + q4'' + q3'') * Ground3>	% angular velocity
V2pts(Ground, RFoot, RAnkle, RFootO);		
A2pts(Ground, RFoot, RAnkle, RFootO);	

%-----------------------------------
% LThigh
%-----------------------------------
% define kinematics (hinge joint, rotation q7 on z axis)
Simprot(Trunk, LThigh, 3, q7)

% define mass properties and CM
Mass     LThigh = par__ThighMass
Inertia  LThigh, 0, 0, par__ThighInertia
P_Hip_LThighO> = par__ThighCMy*LThigh2>

% define angular velocity, angular acceleration, velocity, acceleration
W_LThigh_Ground> = (q7' + q3') * Ground3>	
ALF_LThigh_Ground> = (q7'' + q3'') * Ground3>	
V2pts(Ground, LThigh, Hip, LThighO);				
A2pts(Ground, LThigh, Hip, LThighO);		

%-------------------
% LShank
%-------------------
% define kinematics (hinge joint, rotation q8 on z axis)
Simprot(LThigh, LShank, 3, q8)
P_Hip_LKnee> = -par__ThighLen*LThigh2>;	% define the knee position
V2pts(Ground, LThigh, Hip, LKnee);						% its velocity is needed also
A2pts(Ground, LThigh, Hip, LKnee);						% and acceleration

% define mass properties and CM
Mass     LShank = par__ShankMass
Inertia  LShank, 0, 0, par__ShankInertia	
P_LKnee_LShankO> = par__ShankCMy*LShank2>

% define angular velocity, angular acceleration, velocity, acceleration
W_LShank_Ground> = (q8' + q7' + q3') * Ground3>		
ALF_LShank_Ground> = (q8'' + q7'' + q3'') * Ground3>		
V2pts(Ground, LShank, LKnee, LShankO);					
A2pts(Ground, LShank, LKnee, LShankO);		

%-----------------------------------
% LFoot
%-----------------------------------
% define kinematics (hinge joint, rotation q6 on z axis)
Simprot(LShank, Lfoot, 3, q9)
P_LKnee_LAnkle> = -par__ShankLen*LShank2>;		% define the ankle position
V2pts(Ground, LShank, LKnee, LAnkle);			% its velocity is needed also
A2pts(Ground, LShank, LKnee, LAnkle);			% and acceleration

% define mass properties and CM
Mass     LFoot = par__FootMass
Inertia  LFoot, 0, 0, par__FootInertia
P_LAnkle_LFootO> = par__FootCMx*LFoot1> + par__FootCMy*LFoot2>

% define angular velocity, angular acceleration, velocity, acceleration
W_LFoot_Ground> = (q9' + q8' + q7' + q3') * Ground3>	% angular velocity
ALF_LFoot_Ground> = (q9'' + q8'' + q7'' + q3'') * Ground3>	% angular velocity
V2pts(Ground, LFoot, LAnkle, LFootO);		
A2pts(Ground, LFoot, LAnkle, LFootO);	

%--------------------------------------------------------------------------
% Produce forward kinematics output for the segment positions and orientations
%--------------------------------------------------------------------------
% For each segment, generate a row vector T of 6 elements containing the XY components of orientation and position
TTrunk= [                                              					     	&
	dot(P_GroundO_Hip>, Ground1>), dot(P_GroundO_Hip>, Ground2>),               &   % origin of trunk is hip
	Ground_Trunk[1,1], Ground_Trunk[1,2], Ground_Trunk[2,1], Ground_Trunk[2,2]]

TRThigh= [                                              				    	&
	dot(P_GroundO_Hip>, Ground1>), dot(P_GroundO_Hip>, Ground2>),               &
	Ground_RThigh[1,1], Ground_RThigh[1,2], Ground_RThigh[2,1], Ground_RThigh[2,2]]

TRShank= [                                              				    	&
	dot(P_GroundO_RKnee>, Ground1>), dot(P_GroundO_RKnee>, Ground2>),           &
	Ground_RShank[1,1], Ground_RShank[1,2], Ground_RShank[2,1], Ground_RShank[2,2]]

TRFoot = [                                              						&
	dot(P_GroundO_RAnkle>, Ground1>), dot(P_GroundO_RAnkle>, Ground2>),         &
	Ground_RFoot[1,1], Ground_RFoot[1,2], Ground_RFoot[2,1], Ground_RFoot[2,2]]

TLThigh= [                                              				    	&
	dot(P_GroundO_Hip>, Ground1>), dot(P_GroundO_Hip>, Ground2>),               &
	Ground_LThigh[1,1], Ground_LThigh[1,2], Ground_LThigh[2,1], Ground_LThigh[2,2]]

TLShank= [                                              				    	&
	dot(P_GroundO_LKnee>, Ground1>), dot(P_GroundO_LKnee>, Ground2>),           &
	Ground_LShank[1,1], Ground_LShank[1,2], Ground_LShank[2,1], Ground_LShank[2,2]]

TLFoot = [																		&
	dot(P_GroundO_LAnkle>, Ground1>), dot(P_GroundO_LAnkle>, Ground2>),         &
	Ground_LFoot[1,1], Ground_LFoot[1,2], Ground_LFoot[2,1], Ground_LFoot[2,2]]

% put them all together in one long column vector
fk = [Transpose(TTrunk) ; Transpose(TRThigh) ; Transpose(TRShank) ; Transpose(TRFoot) ; &
      Transpose(TLThigh) ; Transpose(TLShank) ; Transpose(TLFoot)];

% and compute the derivatives
dfk_dq = ZEE(D(fk,q))
fkdot = ZEE(DT(fk))
dfkdot_dq = ZEE(D(fkdot,q))
Encode fk, dfk_dq, fkdot, dfkdot_dq

Code Algebraic() gait2dc_FK_al_raw.c

% EXIT