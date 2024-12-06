% gait2dc_dynamics.al
%
% Equations of motion for the 2D gait model with directional deformable contact
%
% Author: Ton van den Bogert
% Last revised: 02/07/2011

PAUSE 0
AUTOZ ON
OVERWRITE ALL

RUN gait2dc_FK.al

% segment reference frames:
% origin is at center of mass of segment, X anterior while neutral standing

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
Zero = (Fr() + FrStar())/100		% divide by 100 because these are moment imbalances, quite large

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
Code Algebraic() gait2dc_dynamics_al_raw.c

EXIT
