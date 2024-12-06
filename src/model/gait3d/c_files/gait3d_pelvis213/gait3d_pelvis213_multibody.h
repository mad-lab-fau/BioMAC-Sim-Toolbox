// gait3d_pelvis213_multibody.h
// This file defines the data structures that are needed in Autolev multibody model

#define NDOF 33			// number of kinematic degrees of freedom
#define NFKSEG 20		// number of Opensim body segments (excluding ground)
#define NFK 12*NFKSEG	// number of forward kinematic outputs: 12 for each segment: Ox,Oy,Oz,R11,R12,R13,R21,R22,R23,R31,R32,R33
#define NGRFSEG 4		// number of body segments that have contact points attached to them (Rcalc, Rtoes, Lcalc, Ltoes)
#define NGRF 6*NGRFSEG	// number of GRF inputs to equation of motion (3D force and moment on each segment that has contact points)

// macro to declare inertial properties of a segment S
#define INERTIAL(S)  double S##_M, S##_CMx, S##_CMy, S##_CMz, S##_Ixx, S##_Iyy, S##_Izz, S##_Ixy, S##_Iyz, S##_Ixz

typedef struct {
	// Gravity and body weight
	double gravity_x, gravity_y, gravity_z;
	double bodyweight;
	
	// Air drag
	double airdrag;		// drag coefficient, N/(m/s)^2
	double wind;		// wind speed, m/s (positive value is a tailwind)

	// Body segment parameters
	INERTIAL(pelvis);
	INERTIAL(torso);
	
	INERTIAL(Rfemur);
	INERTIAL(Rtibia);
	INERTIAL(Rcalcaneus);
	INERTIAL(Rtoes);
	INERTIAL(Rhumerus);
	INERTIAL(Rulna);
	INERTIAL(Rradius);
	INERTIAL(Rhand);
	
	INERTIAL(Lfemur);
	INERTIAL(Ltibia);
	INERTIAL(Lcalcaneus);
	INERTIAL(Ltoes);
	INERTIAL(Lhumerus);
	INERTIAL(Lulna);
	INERTIAL(Lradius);
	INERTIAL(Lhand);
	
	// Joint axis orientations (are assumed to be the same left and right, with mirroring relative to local XY plane)
	double Rankle1,    Rankle2,    Rankle3;
	double Rsubtalar1, Rsubtalar2, Rsubtalar3;
	double Rmtp1,                  Rmtp3;			// mtp2 is hard zero (as coded in addleg function in gait3d_pelvis213make.m)
	double Relbow1,    Relbow2,    Relbow3;
	double Rraduln1,   Rraduln2,   Rraduln3;
	
    double Lankle1,    Lankle2,    Lankle3;
	double Lsubtalar1, Lsubtalar2, Lsubtalar3;
	double Lmtp1,                  Lmtp3;			// mtp2 is hard zero (as coded in addleg function in gait3d_pelvis213make.m)
	double Lelbow1,    Lelbow2,    Lelbow3;
	double Lraduln1,   Lraduln2,   Lraduln3;
    
	// Joint axis positions (parameters that were coded as hard zeros in gait3d_pelvis213make.m are not included here)
	double back_x, back_y;
	
	double Rhip_x, Rhip_y, Rhip_z;
	double Rknee_x1, Rknee_x2, Rknee_x3, Rknee_x4, Rknee_x5;	// polynomial coefficients for knee x as function of knee angle
	double Rknee_y1, Rknee_y2, Rknee_y3, Rknee_y4, Rknee_y5;	// polynomial coefficients for knee y as function of knee angle
	double Rankle_y;
	double Rsubtalar_x, Rsubtalar_y, Rsubtalar_z;
	double Rmtp_x, Rmtp_y, Rmtp_z;
	double Racromial_x, Racromial_y, Racromial_z;	
	double Relbow_x, Relbow_y, Relbow_z;
	double Rraduln_x, Rraduln_y, Rraduln_z;
	double Rwrist_x, Rwrist_y, Rwrist_z;	
	
	double Lhip_x, Lhip_y, Lhip_z;
	double Lknee_x1, Lknee_x2, Lknee_x3, Lknee_x4, Lknee_x5;	// polynomial coefficients for knee x as function of knee angle
	double Lknee_y1, Lknee_y2, Lknee_y3, Lknee_y4, Lknee_y5;	// polynomial coefficients for knee y as function of knee angle
	double Lankle_y;
	double Lsubtalar_x, Lsubtalar_y, Lsubtalar_z;
	double Lmtp_x, Lmtp_y, Lmtp_z;
	double Lacromial_x, Lacromial_y, Lacromial_z;	
	double Lelbow_x, Lelbow_y, Lelbow_z;
	double Lraduln_x, Lraduln_y, Lraduln_z;
	double Lwrist_x, Lwrist_y, Lwrist_z;
		
} param_struct;

// prototype for the Autolev C function for multibody dynamics
void gait3d_pelvis213_al(param_struct* par, 
	double q[NDOF], 
	double qd[NDOF], 
	double qdd[NDOF],
	double G[NGRF],
	double f[NDOF],
	double df_dq[NDOF][NDOF], 
	double df_dqd[NDOF][NDOF],
	double df_dqdd[NDOF][NDOF],
	double df_dG[NDOF][NGRF],
	double fk[NFK],
	double dfk_dq[NFK][NDOF],
	double fkdot[NFK],
	double dfkdot_dq[NFK][NDOF]);

void gait3d_pelvis213_NoDer_al(param_struct* par,
	double q[NDOF],
	double qd[NDOF],
	double qdd[NDOF],
	double G[NGRF],
	double f[NDOF],
	double fk[NFK],
	double dfk_dq[NFK][NDOF],
	double fkdot[NFK],
	double dfkdot_dq[NFK][NDOF]);

void gait3d_pelvis213_FK_al(param_struct* par,
	double q[NDOF],
	double qd[NDOF],
	double qdd[NDOF],
	double fk[NFK],
	double dfk_dq[NFK][NDOF],
	double fkdot[NFK],
	double dfkdot_dq[NFK][NDOF]);
