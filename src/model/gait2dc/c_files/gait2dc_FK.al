% gait2dc_FK.al
%
% Generate forward kinematics for the 2D gait model
%
% Author: Ton van den Bogert
% Last revised: 01/25/2025

% generate the kinematic expressions
RUN gait2dc_kin.al

%--------------------------------------------------------------------------
% Produce forward kinematics output for the segment positions and orientations
%--------------------------------------------------------------------------
% For each segment, generate a row vector T of 6 elements containing the XY position
% of the proximal joint, and the rotation matrix.
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
