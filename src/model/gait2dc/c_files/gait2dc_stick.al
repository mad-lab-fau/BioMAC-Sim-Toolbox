% gait2dc_stick.al
%
% Equations of motion for the 2D gait model with directional deformable contact
%
% Author: Ton van den Bogert
% Last revised: 02/07/2011

PAUSE 0
AUTOZ ON
OVERWRITE ALL

RUN gait2dc_FK.al

%---------------------------------------------------------------------------------------
% Produce output for a 6-point stick figure: 
% trunk, hip, rknee, rankle, rheel, rtoe, lknee, lankle
%---------------------------------------------------------------------------------------

Stick = [	dot(P_GroundO_TrunkO>, Ground1>) , dot(P_GroundO_TrunkO>, Ground2>) ; &
			dot(P_GroundO_Hip>,    Ground1>) , dot(P_GroundO_Hip>,    Ground2>) ; &
			dot(P_GroundO_RKnee>,  Ground1>) , dot(P_GroundO_RKnee>,  Ground2>) ; &
			dot(P_GroundO_RAnkle>, Ground1>) , dot(P_GroundO_RAnkle>, Ground2>) ; &
			dot(P_GroundO_LKnee>,  Ground1>) , dot(P_GroundO_LKnee>,  Ground2>) ; &
			dot(P_GroundO_LAnkle>, Ground1>) , dot(P_GroundO_LAnkle>, Ground2>) ]

encode Stick

Code Algebraic() gait2dc_stick_al_raw.c

EXIT
