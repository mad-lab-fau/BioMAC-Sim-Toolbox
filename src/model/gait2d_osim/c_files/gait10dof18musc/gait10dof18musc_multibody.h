// gait10dof18musc_multibody.h
// This file defines the data structures that are needed in Autolev multibody model

#define NDOF 10 		// number of kinematic degrees of freedom
#define NFKSEG 8		// number of Opensim body segments (excluding ground)
#define NFK 6*NFKSEG	// number of forward kinematic outputs: 6 for each segment: Ox,Oy,R11,R12,R21,R22
#define NGRFSEG 2		// number of body segments that have contact points attached to them (Rcalc, Rtoes, Lcalc, Ltoes)
#define NGRF 3*NGRFSEG	// number of GRF inputs to equation of motion (2D force and moment on each segment that has contact points)

// macro to declare inertial properties of a segment S
#define INERTIAL(S)  double S##_M, S##_CMx, S##_CMy, S##_Izz

typedef struct {
	// Gravity and body weight
	double gravity_x, gravity_y;
	double bodyweight;
	
	// Air drag
	double airdrag;		// drag coefficient, N/(m/s)^2
	double wind;		// wind speed, m/s (positive value is a tailwind)

	// Body segment parameters
	INERTIAL(pelvis);
	INERTIAL(torso);
    
	INERTIAL(Rfemur);
	INERTIAL(Rtibia);
    INERTIAL(Rfoot);
//	INERTIAL(Rcalcaneus); // Does not have toes or other stuff atm
//	INERTIAL(Rtoes);
	
	INERTIAL(Lfemur);
	INERTIAL(Ltibia);
	INERTIAL(Lfoot);
//  INERTIAL(Lcalcaneus);
//	INERTIAL(Ltoes);
	
	// Alignment parameters for BKA study
	double trans, flex;
    
	// Joint axis positions (parameters that were coded as hard zeros in gait3d_pelvis213make.m are not included here)
	double back_x, back_y;
	
	double Rhip_x, Rhip_y;
	double Rknee_x1, Rknee_x2, Rknee_x3, Rknee_x4, Rknee_x5;	// polynomial coefficients for knee x as function of knee angle
	double Rknee_y1, Rknee_y2, Rknee_y3, Rknee_y4, Rknee_y5;	// polynomial coefficients for knee y as function of knee angle
	double Rankle_y;

	
	double Lhip_x, Lhip_y;
	double Lknee_x1, Lknee_x2, Lknee_x3, Lknee_x4, Lknee_x5;	// polynomial coefficients for knee x as function of knee angle
	double Lknee_y1, Lknee_y2, Lknee_y3, Lknee_y4, Lknee_y5;	// polynomial coefficients for knee y as function of knee angle
	double Lankle_y;
} param_struct;

// prototype for the Autolev C function for multibody dynamics
void gait10dof18musc_al(param_struct* par, 
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

void gait10dof18musc_NoDer_al(param_struct* par,
	double q[NDOF],
	double qd[NDOF],
	double qdd[NDOF],
	double G[NGRF],
	double f[NDOF],
	double fk[NFK],
	double dfk_dq[NFK][NDOF],
	double fkdot[NFK],
	double dfkdot_dq[NFK][NDOF]);

void gait10dof18musc_FK_al(param_struct* par,
	double q[NDOF],
	double qd[NDOF],
	double qdd[NDOF],
	double fk[NFK],
	double dfk_dq[NFK][NDOF],
	double fkdot[NFK],
	double dfkdot_dq[NFK][NDOF]);
