% gait2dc_dynamics.al
%
% Equations of motion for the 2D gait model with directional deformable contact
%
% Author: Ton van den Bogert
% Last revised: 02/07/2011

PAUSE 0
AUTOZ ON
OVERWRITE ALL

% first, generate the kinematic expressions we need, and also the FK outputs
RUN gait2dc_FK.al

% Define model parameters as variables, so they can be modified from Matlab
Constants par__TrunkMass, par__TrunkInertia, par__TrunkCMy
Constants par__ThighMass, par__ThighInertia, par__ThighCMy
Constants par__ShankMass, par__ShankInertia, par__ShankCMy
Constants par__FootMass, par__FootInertia, par__FootCMx, par__FootCMy
Constants par__bodyweight, par__airdrag, par__wind, par__slope
Constants par__gravity					% Local gravitational acceleration value

% Generate expressions for the position, velocity, and acceleration
% of each segment's origin, which is at the center of mass.
% We also assign mass properties here.

%-----------------------------------
% Trunk
%-----------------------------------
Mass     Trunk = par__TrunkMass
Inertia  Trunk, 0, 0, par__TrunkInertia
P_Hip_TrunkO> = par__TrunkCMy*Trunk2>				
V2pts(Ground, Trunk, Hip, TrunkO);		
A2pts(Ground, Trunk, Hip, TrunkO);		

%-----------------------------------
% RThigh
%-----------------------------------
Mass     RThigh = par__ThighMass
Inertia  RThigh, 0, 0, par__ThighInertia
P_Hip_RThighO> = par__ThighCMy*RThigh2>	
V2pts(Ground, RThigh, Hip, RThighO);				
A2pts(Ground, RThigh, Hip, RThighO);		

%-----------------------------------
% RShank
%-----------------------------------
Mass     RShank = par__ShankMass
Inertia  RShank, 0, 0, par__ShankInertia	
P_RKnee_RShankO> = par__ShankCMy*RShank2>	
V2pts(Ground, RShank, RKnee, RShankO);					
A2pts(Ground, RShank, RKnee, RShankO);		

%-----------------------------------
% RFoot
%-----------------------------------
Mass     RFoot = par__FootMass
Inertia  RFoot, 0, 0, par__FootInertia
P_RAnkle_RFootO> = par__FootCMx*RFoot1> + par__FootCMy*RFoot2>
V2pts(Ground, RFoot, RAnkle, RFootO);		
A2pts(Ground, RFoot, RAnkle, RFootO);	

%-----------------------------------
% LThigh
%-----------------------------------
Mass     LThigh = par__ThighMass
Inertia  LThigh, 0, 0, par__ThighInertia
P_Hip_LThighO> = par__ThighCMy*LThigh2>
V2pts(Ground, LThigh, Hip, LThighO);				
A2pts(Ground, LThigh, Hip, LThighO);		

%-------------------
% LShank
%-------------------
Mass     LShank = par__ShankMass
Inertia  LShank, 0, 0, par__ShankInertia	
P_LKnee_LShankO> = par__ShankCMy*LShank2>		
V2pts(Ground, LShank, LKnee, LShankO);					
A2pts(Ground, LShank, LKnee, LShankO);		

%-----------------------------------
% LFoot
%-----------------------------------
Mass     LFoot = par__FootMass
Inertia  LFoot, 0, 0, par__FootInertia
P_LAnkle_LFootO> = par__FootCMx*LFoot1> + par__FootCMy*LFoot2>
V2pts(Ground, LFoot, LAnkle, LFootO);		
A2pts(Ground, LFoot, LAnkle, LFootO);	

%--------------------------------------------------------------------
% Apply gravity
%--------------------------------------------------------------------
Gravity( par__gravity*(-sin(par__slope)*Ground1> - cos(par__slope)*Ground2>) )	

%--------------------------------------------------------------------
% Apply air drag force to the Trunk CM
%--------------------------------------------------------------------
	airspeed> = V_TrunkO_Ground> - par__wind * Ground1>			% velocity of trunk CM, relative to the air
	Force(GroundO/TrunkO, -par__airdrag*mag(airspeed>)*airspeed>)

%--------------------------------------------------------------------------
% Generate and apply ground reaction force and moment on each foot
%--------------------------------------------------------------------------
extforce(RFoot)
extforce(LFoot)
	
%--------------------------------------------------------------------
% Apply joint moments
%--------------------------------------------------------------------
Variables mom{6}		% right hip, knee, ankle, left hip, knee, ankle

Torque(Trunk/RThigh, 	mom1*Ground3>)
Torque(RThigh/RShank, 	mom2*Ground3>)
Torque(RShank/RFoot,	mom3*Ground3>)
Torque(Trunk/LThigh, 	mom4*Ground3>)
Torque(LThigh/LShank, 	mom5*Ground3>)
Torque(LShank/LFoot,	mom6*Ground3>)

%--------------------------------------------------------------------
% Generate equations of motion
%--------------------------------------------------------------------
Zero = (Fr() + FrStar())/100		% divide by 100 because these are moment imbalances, quite large, and we scale them down

%-----------------------------------------------------------------------------------------------------
% For implicit dynamics: generate symbolic expressions for Zero and Jacobians d(Zero)/d(q,qd,qdd,moments,grf)
%-----------------------------------------------------------------------------------------------------
mom = [mom1, mom2, mom3, mom4, mom5, mom6]
grf = [RfootFx,RFootFy,RFootMz,LfootFx,LFootFy,LFootMz]

dz_dq = ZEE(D(Zero,q))
dz_dqd = ZEE(D(Zero,qd))
dz_dqdd = ZEE(D(Zero,qdd))
dz_dmom = ZEE(D(Zero,mom))
dz_dgrf = ZEE(D(Zero,grf))
Encode Zero, dz_dq, dz_dqd, dz_dqdd, dz_dmom, dz_dgrf
Code Algebraic() gait2dc_dyn_al_raw.c

EXIT
