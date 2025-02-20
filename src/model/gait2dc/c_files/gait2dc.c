/**
 *
 * @file gait2dc.c
 * @brief C source code for the MEX function gait2dc.mex32.
 *
 * @author Antonie J. (Ton) van den Bogert, Iris Kellermann, Eva Dorschky
 * @date Mai 31, 2016
 *
 * @copyright 2009-2012 Orchard Kinetics LLC
 *
 * @details Implicit differential equation for 2D musculoskeletal model : f(x,dx/dt,u) = 0
 * with implicit contact formulation.
 * This is the source code for the MEX function gait2cd.mex32
 * The musculoskeletal model is documented in the file gait2dc_reference.odt
 * The function documentation is in gait2dc.m. 
 * The model includes an accelerometer model.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mex.h"
#include "gait2dc.h"

// Compilation settings
#define VERBOSE 0					// set this to 1 for debugging

// size of the model (some other size constants are in gait2dc.h)
#define NMUS 16						// number of muscles 
#define MAXCONTACTS 20				// maximum number of contact points in the model
#define MAXSTATES 2*NDOF+2*NMUS+4*MAXCONTACTS		// max number of system state variables, including contact variables

// M_PI is known in gcc but not in Visual C++ 2008
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

// On some compilers, this macro can detect NaN
#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

// We ignore Jacobian elements smaller than this
#define ALMOSTZERO 1e-15

// define struct that holds muscle properties
typedef struct {
	double Lceopt; 			// Optimal length of CE (m)
	double Width;			// Width of CE force-length relationship (m)
	double Fmax;			// Maximal isometric force of CE (N)
	double Vmax;			// Max. contraction velocity in Lceopt/s
	double Tact, Tdeact;	// activation and deactivation time constants
	double gmax;			// maximum eccentric force
	double SEEslack;		// Slack length of the SEE, in meters
	double PEEslack;		// Slack length of the PEE, relative to Lceopt
	double umax;			// strain of SEE at Fmax load
	double kPEE;			// Stiffness parameter of PEE, in Fmax/Lceopt^2
	double L0;				// Muscle+tendon length when all DOFs are zero
	double MA[NMOM];		// Moment arms (positive when muscle causes ccw rotation of distal joint, model facing +X)
	// the following parameters are derived from other parameters during initialization
	double c1,c2;			// Activation-deactivation rate constants
	double c3;				// Continuity parameter for eccentric force-velocity relationship
	double kSEE;			// Stiffness parameter of SEE, in Fmax/m^2
	double Gmax;			// Max. eccentric force, relative to Fmax
	double Arel;			// Hill's a parameter, relative to Fmax
    double FT;              // Percentage of fast-twitch fibers
} muscleprop;

// global variables
static int initialized = 0;
static mwSize nstates = 0;						// number of state variables
static int ncontacts = 0;					// number of contact points
static double zeros[MAXSTATES];				// contains zeros always
static muscleprop muscles[NMUS];			// contains all muscle properties
static contactprop contacts[MAXCONTACTS];	// contains all contact properties
static param_struct param;					// contains other model parameters
static int nonzeromomentarms;				// number of nonzero moment arms (nonzeros in dL/dq), for memory allocation
static double bodymass;						// total body mass

// ===================================================================================
// ExtractParameters: extract model parameters from the parameter matrix
// ===================================================================================
void ExtractParameters(double par[]) {
	// The parameter matrix is obtained by: xlsread('gait2dc_par.xls').
	// par(1,1) is the first numerical element in the XLS file, normally cell A3 of the sheet.
	// par contains NaN for the text cells in the XLS file	
	int i,j,s;
	#define MAXSECTIONS 10
	int sectionrows[MAXSECTIONS];
	
	// Locate sections of parameter matrix, and store them in proper place here
	int nrows = 0;
	int label;
	while ((label = (int) par[nrows++]) != 999) {
		if (label >= 0 && label < MAXSECTIONS) sectionrows[label] = nrows-1;
	};
	// Read gravity and air drag coefficient from section 1
	s = sectionrows[1] + 1 + nrows;
	param.gravity = par[s++];
	param.airdrag = par[s++];
	param.wind = par[s++];
	param.slope = par[s++]*M_PI/180.0;
	
	// Read body segment parameters from section 1, row 11 and following
	s = sectionrows[1] + 11 + nrows;
	param.TrunkMass 	= par[s++];
	param.ThighMass 	= par[s++];
	param.ShankMass 	= par[s++];
	param.FootMass 		= par[s++];
	s = sectionrows[1] + 11 + 2*nrows;
	param.TrunkInertia 	= par[s++];
	param.ThighInertia 	= par[s++];
	param.ShankInertia 	= par[s++];
	param.FootInertia 	= par[s++];
	s = sectionrows[1] + 14 + 3*nrows;
	param.FootCMx		= par[s++];
	s = sectionrows[1] + 11 + 4*nrows;
	param.TrunkCMy 		= par[s++];
	param.ThighCMy 		= par[s++];
	param.ShankCMy 		= par[s++];
	param.FootCMy 		= par[s++];
	s = sectionrows[1] + 12 + 5*nrows;
	param.ThighLen 		= par[s++];
	param.ShankLen 		= par[s++];
	
	// Read muscle parameters from section 2
	for (i=0; i<NMUS; i++) {
		s = sectionrows[2] + 2 + i + nrows;
		muscles[i].Fmax 	= par[s]; s = s + nrows;
		muscles[i].Lceopt 	= par[s]; s = s + nrows;
		muscles[i].Width 	= par[s]; s = s + nrows;
		muscles[i].PEEslack	= par[s]; s = s + nrows;
		muscles[i].SEEslack	= par[s]; s = s + nrows;
		muscles[i].L0		= par[s]; s = s + nrows;
		for (j=0; j<NMOM; j++) {
			muscles[i].MA[j] = par[s]; s = s + nrows;	
		}
		muscles[i].kPEE		= par[s]; s = s + nrows;
		muscles[i].umax		= par[s]; s = s + nrows;
		muscles[i].Vmax		= par[s]; s = s + nrows;
		muscles[i].Tact		= par[s]; s = s + nrows;
		muscles[i].Tdeact	= par[s]; s = s + nrows;
		muscles[i].Gmax	= par[s]; s = s + nrows;
		muscles[i].Arel	= par[s]; s = s + nrows;
        muscles[i].FT = par[s]; s = s + nrows;
	}

	// Read contact model parameters
	ncontacts = 0;
	s = sectionrows[3] + 2 + nrows;
	while ( (ncontacts < MAXCONTACTS) && !isnan(par[s]) ){
		contacts[ncontacts].foot 	= (int)par[s];
		contacts[ncontacts].x 		= par[s+nrows];
		contacts[ncontacts].y 		= par[s+2*nrows];
		contacts[ncontacts].k1 		= par[s+3*nrows];
		contacts[ncontacts].k2 		= par[s+4*nrows];
		contacts[ncontacts].a 		= par[s+5*nrows]*M_PI/180;
		contacts[ncontacts].c 		= par[s+6*nrows];
		contacts[ncontacts].b 		= par[s+7*nrows];
		ncontacts++;
		s++;
	}
	
	// Read joint range of motion parameters
	for (i=0; i<NMOM; i++) {
		s = sectionrows[4] + 2 + i + nrows;
		param.MinAngle[i] 	= par[s];
		param.MaxAngle[i] 	= par[s+nrows];
		param.JointK2[i]		= par[s+2*nrows];
		param.JointD[i]		= par[s+3*nrows];
		param.JointB[i]		= par[s+4*nrows];
		param.JointK1[i]		= par[s+5*nrows];
		param.JointPhi0[i]	= par[s+6*nrows];
	}	
	
#if VERBOSE
	printf("\nBody segment parameters:\n");
	printf("          Mass          Inertia         CMx         CMy         Length\n");
	printf("----------------------------------------------------------------------\n");
	printf("Trunk   %8.4f      %8.4f          zero      %8.4f        n/a\n", param.TrunkMass,param.TrunkInertia,param.TrunkCMy);
	printf("Thigh   %8.4f      %8.4f          zero      %8.4f     %8.3f\n", param.ThighMass,param.ThighInertia,param.ThighCMy,param.ThighLen);
	printf("Shank   %8.4f      %8.4f          zero      %8.4f     %8.3f\n", param.ShankMass,param.ShankInertia,param.ShankCMy,param.ShankLen);
	printf("Foot    %8.4f      %8.4f      %8.4f      %8.4f        n/a\n", param.FootMass,param.FootInertia,param.FootCMx,param.FootCMy);

	printf("\nMuscle parameters:\n");
	printf("muscle    Fmax    Lceopt   Width PEEslack SEEslack kPEE   umax   Vmax   Tact  Tdeact   Gmax   Arel\n");
	printf("--------------------------------------------------------------------------------------------------\n");
	for (i=0; i<NMUS; i++) {
		printf("%5d    %7.2f %7.4f %7.3f %7.3f %7.4f %7.3f %7.3f %7.2f %7.3f %7.3f %7.3f %7.3f\n", i+1,
			muscles[i].Fmax, muscles[i].Lceopt, muscles[i].Width, muscles[i].PEEslack, 
			muscles[i].SEEslack, muscles[i].kPEE, muscles[i].umax, muscles[i].Vmax, 
			muscles[i].Tact, muscles[i].Tdeact, muscles[i].Gmax, muscles[i].Arel);
	}
	printf("\nMuscle path models:\n");
	printf("muscle      L0	   dRhip	dRknee	dRankle	dLhip	dLknee	dLankle\n");
	printf("-------------------------------------------------------------------\n");
	for (i=0; i<NMUS; i++) {
		printf("%5d    %7.3f", i+1, muscles[i].L0);
		for (j=0; j<NMOM; j++)
			if (muscles[i].MA[j] == 0)
				printf("        ");
			else
				printf(" %7.3f", muscles[i].MA[j]);
		printf("\n");
	}
	printf("\nContact parameters:\n");
	printf("point foot  X       Y       K1      K2      A      C\n");
	printf("----------------------------------------------------\n");
	for (i=0; i<ncontacts; i++) {
		printf("%3d    %d %7.4f %7.4f %e %e %7.4f %7.4f\n", i+1, contacts[i].foot,
			contacts[i].x, contacts[i].y, contacts[i].k1, contacts[i].k2, contacts[i].a, contacts[i].c);
	}
	
	printf("\nJoint range of motion parameters:\n");
	printf("Joint    Min.angle   Max.angle   param.JointK2   param.JointD   param.JointB   param.JointK1   param.JointPhi0\n");
	printf("--------------------------------------------------------------------------------\n");
	for (i=0; i<NMOM; i++) {
		printf("%5d    %8.2f  %8.2f  %8.3f  %8.3f  %8.3f  %8.3f  %8.2f\n", i+1,
			param.MinAngle[i], param.MaxAngle[i], param.JointK2[i], param.JointD[i], 
				param.JointB[i], param.JointK1[i], param.JointPhi0[i]);
	}
#endif

// Preprocessing and error checking
	bodymass = param.TrunkMass + 2*(param.ThighMass + param.ShankMass + param.FootMass);
	param.bodyweight = param.gravity*bodymass;

	// convert range of motion limits to radians and ensure there is a range
	for (i=0; i<NMOM; i++) {
		param.MinAngle[i] 	= param.MinAngle[i]*M_PI/180.0;
		param.MaxAngle[i] 	= param.MaxAngle[i]*M_PI/180.0;
		if (param.MinAngle[i] > param.MaxAngle[i]) {
			printf("Error in joint %d\n", i);
			mexErrMsgTxt("Max angle must be greater than min angle.");
		}
		param.JointPhi0[i] 	= param.JointPhi0[i]*M_PI/180.0;
		param.JointD[i] 		= param.JointD[i]*M_PI/180.0;
	}
	
	// Determine the number of nonzero moment arms, for memory allocation
	// and do other muscle related preprocessing
	nonzeromomentarms = 0;
	for (i=0; i<NMUS; i++) {
		for (j=0; j<NMOM; j++) {
			if (muscles[i].MA[j] != 0) nonzeromomentarms++;
		}
		muscles[i].kSEE = 1.0/(muscles[i].umax*muscles[i].umax*muscles[i].SEEslack*muscles[i].SEEslack);
		muscles[i].c2 = 1.0/muscles[i].Tdeact;
		muscles[i].c1 = 1.0/muscles[i].Tact - muscles[i].c2;
		if (muscles[i].c2 < 0 || muscles[i].c1 < 0) {
			printf("Error in muscle %d\n", i);
			mexErrMsgTxt("Muscle time constants must be positive, and deactivation > activation.");
		}
		muscles[i].c3 = muscles[i].Vmax * muscles[i].Arel * (muscles[i].Gmax - 1.0) / (muscles[i].Arel + 1.0);
	}

	// convert stiffness parameters k1,k2,a into stiffness matrix in local reference frame
	for (i=0; i<ncontacts; i++) {
		// code was derived with Autolev: 
		/*
			variables k1,k2,a
			R = [cos(a),sin(a);-sin(a),cos(a)]
			K = inv(R)*[k1,0;0,k2]*R
			-> (5) K[1,1] = k1*COS(a)^2 + k2*SIN(a)^2
			-> (6) K[1,2] = SIN(a)*COS(a)*(k1-k2)
			-> (7) K[2,1] = SIN(a)*COS(a)*(k1-k2)
			-> (8) K[2,2] = k2 + SIN(a)^2*(k1-k2)
		*/
		double k1 = contacts[i].k1/param.bodyweight;		// internally, all forces will be in units of BW
		double k2 = contacts[i].k2/param.bodyweight;
		double a =  contacts[i].a;
		double sa2 = sin(a)*sin(a);
		double ca2 = cos(a)*cos(a);
		
		contacts[i].kxx = k1*ca2 + k2*sa2;
		contacts[i].kxy = sin(a)*cos(a)*(k1-k2);
		contacts[i].kyy = k2*ca2 + k1*sa2;
		contacts[i].b = contacts[i].b/param.bodyweight;
		contacts[i].LambdaX = 1e-6;					// default value
		contacts[i].LambdaY = 1e-6;					// default value

		// normally, k2 is greater than 1, if not, user wants the old contact model
		if (contacts[i].k2 > 1.0) {
			contacts[i].oldcontactmodel = 0.0;
		}
		else {
			contacts[i].oldcontactmodel = 1.0;
		}
		
		// set the parameters k, d0, bo, v0, gamma for the old contact model
		// (this is an undocumented feature!)
		contacts[i].k = contacts[i].k1/param.bodyweight;
		contacts[i].d0 = contacts[i].k2;
		contacts[i].bo = contacts[i].a*180/M_PI;		// undo the radians conversion
		contacts[i].v0 = contacts[i].c;
		contacts[i].gamma = contacts[i].b*param.bodyweight;        // undo the divide by body weight
	}
}

// ================================================================================
//  Softpos: soft positive function y = Softpos(x) and its derivative
// ================================================================================
double Softpos(double x, double x0, double *dydx) {
	double s,y;
	s = sqrt(x*x + x0*x0);
	y = 0.5*(s + x);
	*dydx = 0.5*(x/s + 1.0); 
	return y;
}

// ===================================================================================
// MusclePath: calculate muscle-tendon length and its derivatives w.r.t. joint angles
// ===================================================================================
double MusclePath(muscleprop *m, double ang[NMOM], double dLm_dang[NMOM]) {

	// this is the 2D model with constant moment arms (which are stored in MA)
	// in the 3D model we will have polynomial models for Lm(ang) and
	
	int i;
	double Lm;
	
	Lm = m->L0;
	for (i=0; i<NMOM; i++) {
		Lm = Lm - m->MA[i]*ang[i];
		dLm_dang[i] = -m->MA[i];
	}
	return Lm;
}


// ===================================================================================
// MuscleDynamics: the implicit muscle dynamics, returns f(a,Lce,Lcedot,Lm) and its derivatives
// ===================================================================================
double MuscleDynamics(muscleprop *m,										// muscle parameters (input)
	double a, double Lce, double Lcedot, double Lm,							// the input variables
	double *df_da, double *df_dLce, double *df_dLcedot, double *df_dLm, 	// the gradients (output)
	double *force, double *dforce_dLm, double *dforce_dLce,					// muscle force (output) and derivatives
    double *forceCE, double *dforceCE_da, double *dforceCE_dLce, double *dforceCE_dLcedot, // muscle force of ce (output) and derivatives
	double *powerCE) {														// power generated by CE (only valid if Lcedot satisfies dynamics)

	double F1, dF1_dLce;
	double f,x,k1;
	double F2, dF2_dLcedot;
	double Fpee, dFpee_dLce;
	double Fsee, dFsee_dLce, dFsee_dLm;
	double Fdamp, dFdamp_dLcedot;
	double Fce;
	
	// F1 is the isometric force-length relationship at max activation, normalized to Fmax
	x = (Lce - 1.0)/m->Width;
	F1 = exp(-x*x);
	dF1_dLce = -2.0*x*F1 / m->Width;
	
	// F2 is the normalized force-velocity relationship
	if (Lcedot < 0) {
		// concentric contraction
		x = (m->Vmax - Lcedot/m->Arel);
		F2 = (m->Vmax + Lcedot)/x;
		dF2_dLcedot = (1.0 + F2/m->Arel)/x;
	}
	else {
		// eccentric contraction
		x = Lcedot + m->c3;
		F2 = (m->Gmax*Lcedot + m->c3) / x;
		dF2_dLcedot = (m->Gmax - F2)/x;
	}
	
	// Fdamp is a small viscous damping of CE (0.001 of Fmax at 1 Lceopt/s) to ensure df/dLcedot is never zero
	Fdamp = 0.001*Lcedot;
	dFdamp_dLcedot = 0.001;
	
	// Fpee is the PEE force-length relationship, normalized to Fmax
	k1 = 0.01 * m->Lceopt;			// stiffness of the linear term is 0.01 Fmax/meter
	x = (Lce - m->PEEslack);		// elongation of PEE, in Lceopt units
	Fpee = k1*x;					// low stiffness linear term
	dFpee_dLce = k1;
	if (x>0) {						// add quadratic term for positive deformation
		Fpee = Fpee + m->kPEE*x*x;
		dFpee_dLce = dFpee_dLce + 2*m->kPEE*x;
	}
	
	//  Fsee is the SEE force-length relationship, normalized to Fmax
	k1 = 0.01;										// stiffness of the linear term is 0.01 Fmax/meter	
	x = Lm - Lce*m->Lceopt - m->SEEslack;			// elongation of SEE, in meters
	Fsee = k1*x;									// low stiffness linear term
	dFsee_dLce = -k1*m->Lceopt;
	dFsee_dLm = k1;
	if (x>0) {										// add quadratic term for positive deformation
		Fsee = Fsee + m->kSEE*x*x;
		dFsee_dLce = dFsee_dLce - 2*m->kSEE*m->Lceopt*x;
		dFsee_dLm = dFsee_dLm + 2*m->kSEE*x;
	}

	// Compute f, the force imbalance in the muscle contraction, and its derivatives
	Fce = a*F1*F2 + Fdamp;
	f = Fsee - Fce - Fpee;
	*df_da = -F1*F2;
	*df_dLce = dFsee_dLce - a*dF1_dLce*F2 - dFpee_dLce;
	*df_dLcedot = -a*F1*dF2_dLcedot - dFdamp_dLcedot;
	*df_dLm = dFsee_dLm;
	
	// Muscle force is the force in SEE
	*force = m->Fmax*Fsee;
	*dforce_dLm = m->Fmax * dFsee_dLm;
	*dforce_dLce = m->Fmax * dFsee_dLce;
    
    // Muscle force is the force in CE
	*forceCE = m->Fmax*Fce;
	*dforceCE_da = m->Fmax * F1*F2;
    *dforceCE_dLce = m->Fmax * a*dF1_dLce*F2;
    *dforceCE_dLcedot = m->Fmax * (a*F1*dF2_dLcedot + dFdamp_dLcedot);
	
	// power (W) generated by CE (positive when shortening)
	// Fce was in Fmax units, and Lce was in Lceopt units, so some scaling needed
	*powerCE = -m->Fmax*Fce*Lcedot*m->Lceopt;
	
	// Return the imbalance
	return f;
}

// =========================================================================
// mexFunction: this is the actual MEX function interface
// =========================================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// working variables
	int i, j, jj, k, kk, nrows, ncols;
	char *command, *parametername;
	double *par, ang, angvel, dmax, dmin, ddmax, ddmin;
	double massfactor, lengthfactor, inertiafactor;
	int derivatives;
	
	
	// muscle variables
	double Lm[NMUS];				// muscle+tendon length, based on skeleton kinematic state
	double dLm_dang[NMUS][NMOM];	// derivatives of Lm with respect to joint angles
	double force[NMUS];				// muscle forces
    double forceCE[NMUS];           // muscle forces of CE
	double CEpower[NMUS];			// power generated by CE of each muscle
	double g[NMUS];					// muscle force imbalance
	double dg_da[NMUS], dg_dLce[NMUS], dg_dLcedot[NMUS], dg_dLm[NMUS]; 	// derivatives of muscle imbalance
	double dforce_dLm[NMUS], dforce_dLce[NMUS];
    double dforceCE_da[NMUS], dforceCE_dLce[NMUS], dforceCE_dLcedot[NMUS];

	// multibody dynamics variables
	double *q, *qd, *qdd;
	double mom[NMOM];
	double zero[NDOF];
	double dz_dq[NDOF][NDOF];
	double dz_dqd[NDOF][NDOF];
	double dz_dqdd[NDOF][NDOF];
	double dz_dmom[NDOF][NMOM];
	double dz_dgrf[NDOF][NGRF];
	double dmom_dang[NMOM][NMOM];
	double dmom_dangvel[NMOM];		// assumed to be a diagonal matrix (no coupling between joints)
	double dmom_dLce[NMOM][NMUS];
	double grf[NGRF];
	double dgrfdx[NGRF][4*MAXCONTACTS];		// dense dgrf/dx for use in dynamics calculations
	double fk[NFK], dfk_dq[NFK][NDOF];
	double fkdot[NFK], dfkdot_dq[NFK][NDOF];
	double Stick[NSTICK][2];

	//CHANGE
	// accelerometer model variables
	double acc[42];
	double dacc_dq[42][NDOF];
	double dacc_dqd[42][NDOF];
	double dacc_dqdd[42][NDOF];
	
	// MEX function pointers to inputs from Matlab
	double *x, *xdot, *u, *M;
	
	// pointers for MEX function outputs
	double *f, *f2;
	
	double *df_dx; 
	mwIndex Ndf_dx, *df_dx_irs, *df_dx_jcs;
	mwSize Nalloc_df_dx;
    
    double *dfCE_dx; 
	mwIndex NdfCE_dx, *dfCE_dx_irs, *dfCE_dx_jcs;
	mwSize Nalloc_dfCE_dx;
	
    double *dfCE_dxdot; 
	mwIndex NdfCE_dxdot, *dfCE_dxdot_irs, *dfCE_dxdot_jcs;
	mwSize Nalloc_dfCE_dxdot;
    
	double *df_dxdot;
	mwIndex Ndf_dxdot, *df_dxdot_irs, *df_dxdot_jcs;
	mwSize Nalloc_df_dxdot;
	
	double *df_du;
	mwIndex Ndf_du, *df_du_irs, *df_du_jcs;
	mwSize Nalloc_df_du;
	
	double *df_dM;
	mwIndex Ndf_dM, *df_dM_irs, *df_dM_jcs;
	mwSize Nalloc_df_dM;

	double *dgrf_dx;
	mwIndex Ndgrf_dx, *dgrf_dx_irs, *dgrf_dx_jcs;
	mwSize Nalloc_dgrf_dx;

	double *dM_dx;
	mwIndex NdM_dx, *dM_dx_irs, *dM_dx_jcs;
	mwSize Nalloc_dM_dx;

	double *dF_dx;
	mwIndex NdF_dx, *dF_dx_irs, *dF_dx_jcs;
	mwSize Nalloc_dF_dx;
	
	//CHANGE
	// for accelerometer model
	double *fa, *dfa_dq, *dfa_dqd, *dfa_dqdd;
	mwIndex Ndfa_dq;
	mwIndex *dfa_dq_irs, *dfa_dq_jcs;
	mwIndex Ndfa_dqd;
	mwIndex *dfa_dqd_irs, *dfa_dqd_jcs;
	mwIndex Ndfa_dqdd;
	mwIndex *dfa_dqdd_irs, *dfa_dqdd_jcs;
	
	// the number of nonzeros of the sparse accelerometer Jacobians is already known:
	#define NNZ_DA_DQ 	38
	#define NNZ_DA_DQD 	58
	#define NNZ_DA_DQDD 86	
    
	//CHANGE
	// command flags
	int CMDdynamics, CMDjointmoments, CMDstick, CMDFKin, CMDmuscleforces, CMDmuscleCEforces, CMDscale;
	int CMDgrf, CMDget, CMDset, CMDmuscleCEpower, CMDaccelerometer;
		
	// get the first argument, see if it is a command string
	if (nrhs < 1)
		mexErrMsgTxt("At least one input argument is needed.");
    if (!mxIsChar(prhs[0])) 
		mexErrMsgTxt("First argument must be a string.");
 	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("First argument must be a row vector.");
	command = mxArrayToString(prhs[0]);
	if(command == NULL) 
	  mexErrMsgTxt("First argument was not a command string.");
	
	//CHANGE
	// decode the command string
	CMDdynamics = 0;
	CMDjointmoments = 0;
	CMDstick = 0;
	CMDFKin = 0;
	CMDmuscleforces = 0;
    CMDmuscleCEforces = 0;
	CMDmuscleCEpower = 0;
	CMDgrf = 0;
	CMDget = 0;
	CMDset = 0;
	CMDaccelerometer = 0;
	CMDscale = 0;
	

	if (strcmp(command, "Dynamics") == 0) {
		CMDdynamics = 1;		// the most common use first
		
		// Check number of inputs
		if ((nrhs < 4) || (nrhs > 5)) {
				mexErrMsgTxt("gait2dc(Dynamics): must have 4 or 5 inputs.");		
		}

		// Check number of outputs, determine whether derivatives are requested
		if (nlhs == 1) {
			derivatives = 0;
		}
		else {
			if (nlhs > 5 || nlhs < 1)
				mexErrMsgTxt("gait2dc(Dynamics): has 1 to 5 outputs.");
			if (nlhs == 5 && nrhs != 5)
				mexErrMsgTxt("gait2dc(Dynamics): 5th output df_dM requires 5th input M");
			derivatives = 1;
		}
	}
	else if (strcmp(command, "Initialize") == 0) {
		if (nrhs != 2)
			mexErrMsgTxt("gait2dc:Initialize: Two input arguments required.");
		nrows = mxGetM(prhs[1]);
		ncols = mxGetN(prhs[1]);
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
			mexErrMsgTxt("gait2dc: Incorrect type for par, must be double.");
		par = mxGetPr(prhs[1]);
		printf("Initializing MEX function gait2dc...\n");
		ExtractParameters(par);
		initialized = 1959;
		nstates = 2*NDOF + 2*NMUS + 4*ncontacts;	
		for (i=0; i<(int)nstates; i++) {
			zeros[i] = 0.0;
		}

		// return a somewhat neutral state to Matlab
		plhs[0] = mxCreateDoubleMatrix(nstates, 1, mxREAL);
		x = mxGetPr(plhs[0]);
		for (i=0; i<(int)nstates; i++) x[i] = 0.0;		// start with zeros everywhere
		x[1] = 1.0;									// vertical position of hip
		for (i=0; i<NMUS; i++) x[2*NDOF+i] = 1.4;	// a neutral LCE value
		for (i=0; i<ncontacts; i++) {
			j = 2*NDOF+2*NMUS+4*i;					// where contact variables are stored for this point
			x[j]   = contacts[i].x;					// put contact point in undeformed state
			x[j+1] = x[1] - param.ThighLen - param.ShankLen + contacts[i].y;
		}
	
		return;
	}
	else if	(strcmp(command, "Jointmoments") == 0) {
		CMDjointmoments = 1;

		// Check number of inputs
		if (nrhs != 2) {
				mexErrMsgTxt("gait2dc(Jointmoments): must have 2 inputs.");		
		}

		// Check number of outputs, determine whether derivatives are requested
		if (nlhs == 1) {
			derivatives = 0;
		}
		else {
			if (nlhs > 2 || nlhs < 1)
				mexErrMsgTxt("gait2dc(Jointmoments): has 1 or 2 outputs.");
			derivatives = 1;
		}
	}
	else if (strcmp(command, "Stick") == 0) 		CMDstick = 1;
	else if (strcmp(command, "Fkin") == 0)  {
		CMDFKin = 1;

		// Check number of inputs
		if (nrhs != 2) {
				mexErrMsgTxt("gait2dc(Fkin): must have 2 inputs.");		
		}

		// Check number of outputs, determine whether derivatives are requested
		if (nlhs == 1) {
			derivatives = 0;
		}
		else {
			if (nlhs > 3 || nlhs < 1)
				mexErrMsgTxt("gait2dc(Fkin): has 1, 2, or 3 outputs.");
			derivatives = 1;
		}
	}		
	else if	(strcmp(command, "Muscleforces") == 0) 	{
		CMDmuscleforces = 1;

		// Check number of inputs
		if (nrhs != 2) {
				mexErrMsgTxt("gait2dc(Muscleforces): must have 2 inputs.");		
		}

		// Check number of outputs, determine whether derivatives are requested
		if (nlhs == 1) {
			derivatives = 0;
		}
		else {
			if (nlhs > 2 || nlhs < 1)
				mexErrMsgTxt("gait2dc(Muscleforces): has 1 or 2 outputs.");
			derivatives = 1;
		}
	}
    else if	(strcmp(command, "MuscleCEforces") == 0) 	{
		CMDmuscleCEforces = 1;

		// Check number of inputs
		if (nrhs != 3) {
				mexErrMsgTxt("gait2dc(MuscleCEforces): must have 3 inputs.");		
		}

		// Check number of outputs, determine whether derivatives are requested
		if (nlhs == 1) {
			derivatives = 0;
		}
		else {
			if (nlhs == 2 || nlhs > 3)
				mexErrMsgTxt("gait2dc(MuscleCEforces): has 1 or 3 outputs.");
			derivatives = 1;
		}
	}
	else if	(strcmp(command, "MuscleCEpower") == 0) CMDmuscleCEpower = 1;
	else if (strcmp(command, "GRF") == 0) {
		CMDgrf = 1;

		// Check number of inputs
		if (nrhs != 2) {
				mexErrMsgTxt("gait2dc(GRF): must have 2 inputs.");		
		}

		// Check number of outputs, determine whether derivatives are requested
		if (nlhs == 1) {
			derivatives = 0;
		}
		else {
			if (nlhs > 2 || nlhs < 1)
				mexErrMsgTxt("gait2dc(GRF): has 1 or 2 outputs.");
			derivatives = 1;
		}
	}
	else if (strcmp(command, "Get") == 0) {
		CMDget = 1;
		if (nrhs != 2) {
			mexErrMsgTxt("gait2dc: Get: exactly two input arguments needed.");	
		}
		if (nlhs != 1) {
			mexErrMsgTxt("gait2dc: Get: exactly one output arguments needed.");	
		}
		if (!mxIsChar(prhs[1])) 
			mexErrMsgTxt("gait2dc: Get: second argument must be a string.");
		if (mxGetM(prhs[1])!=1)
		mexErrMsgTxt("gait2dc: Get second argument must be a string.");
		parametername= mxArrayToString(prhs[1]);
		if (parametername == NULL) 
			mexErrMsgTxt("gait2dc: Get: second argument was not a command string.");
	}
	else if (strcmp(command, "Set") == 0) {
		CMDset = 1;
		if (nrhs < 2 || nrhs > 4) {
			mexErrMsgTxt("gait2dc: Set: between two and four input arguments needed.");	
		}
		if (!mxIsChar(prhs[1])) 
			mexErrMsgTxt("gait2dc: Set: second argument must be a string.");
		if (mxGetM(prhs[1])!=1)
		mexErrMsgTxt("gait2dc: Set: second argument must be a string.");
		parametername= mxArrayToString(prhs[1]);
		if (parametername == NULL) 
			mexErrMsgTxt("gait2dc: Set: second argument was not a command string.");
	}
	//CHANGE
	// If command is "Accelerometer", we need q,qdot,qdotdot then we can finish up
	else if(strcmp(command, "Accelerometer") == 0) {
		CMDaccelerometer = 1;
    }

	//CHANGE
	else if(strcmp(command, "Scale") == 0) {
		CMDscale = 1;
	}
	else mexErrMsgTxt("gait2dc: Command not recognized.");
	
		
	// Not initializing, so make sure that Initialize was done already
	if (initialized != 1959) {
		mexErrMsgTxt("gait2dc: model was not initialized.");
	}
	
	//CHANGE
	if (CMDscale) {
		if (nrhs != 3) {
			mexErrMsgTxt("gait2d: Scale: exactly three input arguments needed.");				
		}
		
		// make sure we have initialized with the original parameters from Winter, and model was not scaled yet:
		if ((fabs(param.TrunkMass - 50.85)>1e-6) || (fabs(param.TrunkCMy - 0.315504)>1e-6)) {
			printf("trunk mass: %f\n",param.TrunkMass);
			printf("trunk CM: %f\n",param.TrunkCMy);
			mexErrMsgTxt("gait2d: Scale: model was already scaled.");				
		}
		
		// check the input arguments
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect type for Body Mass, must be double.");
		}
		nrows = mxGetM(prhs[1]);
		ncols = mxGetN(prhs[1]);
		if ((nrows != 1) || (ncols != 1) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect size for Body Mass, must be a scalar.");
		}
		if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect type for Body Mass, must be double.");
		}
		nrows = mxGetM(prhs[2]);
		ncols = mxGetN(prhs[2]);
		if ((nrows != 1) || (ncols != 1) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect size for Body Mass, must be a scalar.");
		}
		
		// get the body mass and height and calculate scaling factors
		massfactor = *mxGetPr(prhs[1]) / 75.0;		// divide subject mass by Winter standard body mass
		lengthfactor = *mxGetPr(prhs[2]) / 1.8;		// divide subject height by Winter standard body height
		inertiafactor = massfactor * lengthfactor * lengthfactor;

		// apply the scaling factors
		param.TrunkMass *= massfactor;
		param.ThighMass *= massfactor;
		param.ShankMass *= massfactor;
		param.FootMass *= massfactor;

		param.TrunkInertia *= inertiafactor;
		param.ThighInertia *= inertiafactor;
		param.ShankInertia *= inertiafactor;
		param.FootInertia *= inertiafactor;

		param.ThighLen *= lengthfactor;
		param.ShankLen *= lengthfactor;

		param.TrunkCMy *= lengthfactor;
		param.ThighCMy *= lengthfactor;
		param.ShankCMy *= lengthfactor;
		param.FootCMx *= lengthfactor;
		param.FootCMy *= lengthfactor;

		//CHANGE
		//param.ContactY0 *= lengthfactor;
		//param.ContactHeelX *= lengthfactor;
		//param.ContactToeX  *= lengthfactor;
		
		return;
	}
	
	
	// Is it a "Get" command?
	if (CMDget) {			
		// Decode the parameter name
		if (strcmp(parametername, "Lceopt") == 0) {
			plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].Lceopt;
		}
		else if (strcmp(parametername, "Fmax") == 0) {
			plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].Fmax;
		}
        else if (strcmp(parametername, "FT") == 0) {
            plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].FT;
		}
        else if (strcmp(parametername, "Width") == 0) {
            plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].Width;
		}
        else if (strcmp(parametername, "Ahill") == 0) {
            plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].Arel;
		}
        else if (strcmp(parametername, "Gmax") == 0) {
            plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].Gmax;
		}
        else if (strcmp(parametername, "PEEslack") ==0) {
            plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].PEEslack;
		}
        else if (strcmp(parametername, "kPEE") ==0) {
            plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].kPEE;
		}
        
		else if (strcmp(parametername, "Total Mass") == 0) {
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			*f = bodymass;
		}
		else {
			printf("Parameter name: %s\n", parametername);
			mexErrMsgTxt("gait2dc: getparameter: parameter name not recognized.");		
		}
		return;
	}

	// Is it a "Set" command?
	if (CMDset) {
		// Decode the parameter name
		if (strcmp(parametername, "Extra Mass") == 0) {
			if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgTxt("gait2dc: Incorrect type for Extra Mass, must be double.");
			nrows = mxGetM(prhs[2]);
			ncols = mxGetN(prhs[2]);
			if ((nrows != 1) || (ncols != 4) )
				mexErrMsgTxt("gait2dc: Incorrect size for Extra Mass, must be a 1x4 matrix.");
			f = mxGetPr(prhs[2]);
			if (f[0] != 3 && f[2] != 0.0) {
				mexErrMsgTxt("gait2dc: On Trunk, Thigh, Shank, extra mass must be at X=0.");
			}
			if (f[0] == 0) {
				printf("Adding %8.3f kg to Trunk at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.TrunkCMy = (param.TrunkMass * param.TrunkCMy + f[1]*f[3])/(param.TrunkMass + f[1]);
				// new mass is sum
				param.TrunkMass = param.TrunkMass + f[1];
				// new inertia according to parallel axes theorem
				param.TrunkInertia = param.TrunkInertia + f[1]*(f[3]-param.TrunkCMy)*(f[3]-param.TrunkCMy);
			}
			else if (f[0] == 1) {
				printf("Adding %8.3f kg to Thigh segments at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.ThighCMy = (param.ThighMass * param.ThighCMy + f[1]*f[3])/(param.ThighMass + f[1]);
				// new mass is sum
				param.ThighMass = param.ThighMass + f[1];
				// new inertia according to parallel axes theorem
				param.ThighInertia = param.ThighInertia + f[1]*(f[3]-param.ThighCMy)*(f[3]-param.ThighCMy);
			}
			else if (f[0] == 2) {
				printf("Adding %8.3f kg to Shank segments at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.ShankCMy = (param.ShankMass * param.ShankCMy + f[1]*f[3])/(param.ShankMass + f[1]);
				// new mass is sum
				param.ShankMass = param.ShankMass + f[1];
				// new inertia according to parallel axes theorem
				param.ShankInertia = param.ShankInertia + f[1]*(f[3]-param.ShankCMy)*(f[3]-param.ShankCMy);
			}
			else if (f[0] == 3) {
				printf("Adding %8.3f kg to Foot segments at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.FootCMx = (param.FootMass * param.FootCMx + f[1]*f[2])/(param.FootMass + f[1]);
				param.FootCMy = (param.FootMass * param.FootCMy + f[1]*f[3])/(param.FootMass + f[1]);
				// new mass is sum
				param.FootMass = param.FootMass + f[1];
				// new inertia according to parallel axes theorem
				param.FootInertia = param.FootInertia 
					+ f[1]*(f[2]-param.FootCMx)*(f[2]-param.FootCMx)
					+ f[1]*(f[3]-param.FootCMy)*(f[3]-param.FootCMy);
			}
			else {
				mexErrMsgTxt("gait2dc: Incorrect segment number in Set Extra Mass.  Must be 1 or 2 or 3.");
			}
			// update the bodymass
			bodymass = param.TrunkMass + 2*(param.ThighMass + param.ShankMass + param.FootMass);
			param.bodyweight = param.gravity*bodymass;
		}
        else if (strcmp(parametername, "Right BKA") == 0) {
            // Get the ankle stiffness, specified by user as 3rd argument of the MEX function
            if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
                mexErrMsgTxt("gait2d: Incorrect type for ankle stiffness, must be double.");
            nrows = mxGetM(prhs[2]);
            ncols = mxGetN(prhs[2]);
            if ((nrows != 1) || (ncols != 1) )
                mexErrMsgTxt("gait2d: Incorrect size for ankle stiffness, must be a scalar.");

            // Get the ankle damping, specified by user as 4th argument of the MEX function
            if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
                mexErrMsgTxt("gait2d: Incorrect type for ankle damping, must be double.");
            nrows = mxGetM(prhs[3]);
            ncols = mxGetN(prhs[3]);
            if ((nrows != 1) || (ncols != 1) )
                mexErrMsgTxt("gait2d: Incorrect size for ankle damping, must be a scalar.");

            // Remove the ankle muscles in the right leg
            for (i=5; i<8; i++) muscles[i].Fmax = 0.0;

            // Passive elastic ankle on right leg
            param.MinAngle[2] = 0.0;
            param.MaxAngle[2] = 0.0;
            param.JointK2[2] = *mxGetPr(prhs[2]);
            param.JointB[2] = *mxGetPr(prhs[3]);
        }
        else if (strcmp(parametername, "Right opt BKA") == 0) {
            // Prosthesis parameters will be optimized, only remove ankle muscles
            for (i=5; i<8; i++) muscles[i].Fmax = 0.0;
        }
        else if (strcmp(parametername, "Right AKA") == 0) {
            // Get the ankle stiffness, specified by user as 3rd argument of the MEX function
            if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
                mexErrMsgTxt("gait2d: Incorrect type for ankle stiffness, must be double.");
            nrows = mxGetM(prhs[2]);
            ncols = mxGetN(prhs[2]);
            if ((nrows != 1) || (ncols != 1) )
                mexErrMsgTxt("gait2d: Incorrect size for ankle stiffness, must be a scalar.");
            f = mxGetPr(prhs[2]);

            // Get the ankle damping, specified by user as 4th argument of the MEX function
            if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
                mexErrMsgTxt("gait2d: Incorrect type for ankle damping, must be double.");
            nrows = mxGetM(prhs[3]);
            ncols = mxGetN(prhs[3]);
            if ((nrows != 1) || (ncols != 1) )
                mexErrMsgTxt("gait2d: Incorrect size for ankle damping, must be a scalar.");

            // Remove the ankle muscles in the right leg
            for (i=2; i<8; i++) muscles[i].Fmax = 0.0;

            // Passive elastic ankle on right leg
            param.MinAngle[2] = 0.0;
            param.MaxAngle[2] = 0.0;
            param.JointK2[2] = *mxGetPr(prhs[2]);
            param.JointB[2] = *mxGetPr(prhs[3]);
        }
		else if (strcmp(parametername, "Lambda") == 0) {
			if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgTxt("gait2dc: Incorrect type for Lambda, must be double.");
			nrows = mxGetM(prhs[2]);
			ncols = mxGetN(prhs[2]);
			if ((nrows != 1) || (ncols != 2) )
				mexErrMsgTxt("gait2dc: Incorrect size for Lambda, must be a 1x2 matrix.");
			f = mxGetPr(prhs[2]);
			for (i=0; i<ncontacts; i++) {
				contacts[i].LambdaX = f[0];
				contacts[i].LambdaY = f[1];
			}
		}
		else {
			printf("Parameter name: %s\n", parametername);
			mexErrMsgTxt("gait2dc: Set: parameter name not recognized.");		
		}
		return;
	}

    
    if (CMDaccelerometer) {
		// get the inputs from Matlab
		if (nrhs < 4) mexErrMsgTxt("gait2d: Accelerometer command needs q,qd,qdd.");
		nrows = mxGetM(prhs[1]);
		ncols = mxGetN(prhs[1]);
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) mexErrMsgTxt("gait2d: Incorrect type for q, must be double.");
		if ((nrows != NDOF) || (ncols != 1)) mexErrMsgTxt("gait2d: Incorrect size for q, must be a 9 x 1 column vector.");
		q = mxGetPr(prhs[1]);
		nrows = mxGetM(prhs[2]);
		ncols = mxGetN(prhs[2]);
		if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) mexErrMsgTxt("gait2d: Incorrect type for qd, must be double.");
		if ((nrows != NDOF) || (ncols != 1)) mexErrMsgTxt("gait2d: Incorrect size for qd, must be a 9 x 1 column vector.");
		qd = mxGetPr(prhs[2]);
		nrows = mxGetM(prhs[3]);
		ncols = mxGetN(prhs[3]);
		if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) mexErrMsgTxt("gait2d: Incorrect type for qd, must be double.");
		if ((nrows != NDOF) || (ncols != 1)) mexErrMsgTxt("gait2d: Incorrect size for qd, must be a 9 x 1 column vector.");
		qdd = mxGetPr(prhs[3]);
		
		// call the Autolev C code
		acc_al(&param, q, qd, qdd, acc, dacc_dq, dacc_dqd, dacc_dqdd);  

		// Create the return argument fa
		if (nlhs < 1) mexErrMsgTxt("gait2d: Accelerometer command needs at least one output.");
		plhs[0] = mxCreateDoubleMatrix(42, 1, mxREAL);
		fa = mxGetPr(plhs[0]);
		for (i=0; i<42; i++) {
			fa[i] = acc[i];
		}
		
		// Create the sparse Jacobian dfa_dq
		if (nlhs < 2) return;
		plhs[1] = mxCreateSparse(42, NDOF, NNZ_DA_DQ, 0);
		dfa_dq = mxGetPr(plhs[1]);
		dfa_dq_irs = mxGetIr(plhs[1]);
		dfa_dq_jcs = mxGetJc(plhs[1]);
		Ndfa_dq = 0;
		for (i=0; i<NDOF; i++) {			// derivatives of fa w.r.t. q, columns 1..NDOF
			dfa_dq_jcs[i] = Ndfa_dq;		// store element number where this column starts			
			for (j=0; j<42; j++) {			// go through all elements of this row
				dfa_dq_irs[Ndfa_dq] = j;			// store row number of this matrix element
				if (dacc_dq[j][i] != 0){				// if this matrix element is not zero
					dfa_dq[Ndfa_dq++] = dacc_dq[j][i];	// store its value
				}
			}
		}
		dfa_dq_jcs[NDOF] = Ndfa_dq;		// store final element number
		if (Ndfa_dq > NNZ_DA_DQ) {
			printf("Number of nonzeros was %d and should be %d\n",Ndfa_dq,NNZ_DA_DQ);
			mexErrMsgTxt("gait2d: dfa_dq has incorrect number of nonzeros.");
		}
		
		// Create the sparse Jacobian dfa_dqd
		if (nlhs < 3) return;
		plhs[2] = mxCreateSparse(42, NDOF, NNZ_DA_DQD, 0);
		dfa_dqd = mxGetPr(plhs[2]);
		dfa_dqd_irs = mxGetIr(plhs[2]);
		dfa_dqd_jcs = mxGetJc(plhs[2]);
		Ndfa_dqd = 0;
		for (i=0; i<NDOF; i++) {			// derivatives of fa w.r.t. q, columns 1..NDOF
			dfa_dqd_jcs[i] = Ndfa_dqd;		// store element number where this column starts			
			for (j=0; j<42; j++) {			// go through all elements of this row
				dfa_dqd_irs[Ndfa_dqd] = j;			// store row number of this matrix element
				if (dacc_dqd[j][i] != 0){					// if this matrix element is not zero
					dfa_dqd[Ndfa_dqd++] = dacc_dqd[j][i];	// store its value
				}
			}
		}
		dfa_dqd_jcs[NDOF] = Ndfa_dqd;		// store final element number
		if (Ndfa_dqd > NNZ_DA_DQD) {
			printf("Number of nonzeros was %d and should be %d\n",Ndfa_dqd,NNZ_DA_DQD);
			mexErrMsgTxt("gait2d: dfa_dqd has incorrect number of nonzeros.");
		}
		
		// Create the sparse Jacobian dfa_dqdd
		if (nlhs < 4) return;
		plhs[3] = mxCreateSparse(42, NDOF, NNZ_DA_DQDD, 0);
		dfa_dqdd = mxGetPr(plhs[3]);
		dfa_dqdd_irs = mxGetIr(plhs[3]);
		dfa_dqdd_jcs = mxGetJc(plhs[3]);
		Ndfa_dqdd = 0;
		for (i=0; i<NDOF; i++) {			// derivatives of fa w.r.t. q, columns 1..NDOF
			dfa_dqdd_jcs[i] = Ndfa_dqdd;		// store element number where this column starts			
			for (j=0; j<42; j++) {			// go through all elements of this row
				
				dfa_dqdd_irs[Ndfa_dqdd] = j;			// store row number of this matrix element
				if (dacc_dqdd[j][i] != 0){						// if this matrix element is not zero
					dfa_dqdd[Ndfa_dqdd++] = dacc_dqdd[j][i];	// store its value
				}
			}
		}
		dfa_dqdd_jcs[NDOF] = Ndfa_dqdd;		// store final element number
		if (Ndfa_dqdd > NNZ_DA_DQDD) {
			printf("Number of nonzeros was %d and should be %d\n",Ndfa_dqdd,NNZ_DA_DQDD);
			mexErrMsgTxt("gait2d: dfa_dqdd has incorrect number of nonzeros.");
		}
		
		return;
    }
    
	
	// State of model is the second argument, for all other uses of this MEX function
	if (nrhs < 2) {
		mexErrMsgTxt("gait2dc: state vector x is needed as second argument.");
	}
	nrows = mxGetM(prhs[1]);
	ncols = mxGetN(prhs[1]);
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
		mexErrMsgTxt("gait2dc: Incorrect type for state vector x, must be double.");
	}
	if ((nrows != nstates) || (ncols != 1) ) {
		printf("number of states: %d\n", nstates);
		mexErrMsgTxt("gait2dc: Incorrect size for state vector x, must be a state column vector.");
	}
	x = mxGetPr(prhs[1]);
	
	// Use zeros for xdot and u, if they are not important
	xdot = zeros;		// xdot now points to zeros
	u = zeros;
	
	// If command is "Dynamics", "MuscleCEpower", "MuscleCEforces", or "accelerometer" we also need xdot
	if (CMDdynamics || CMDmuscleCEpower || CMDmuscleCEforces|| CMDaccelerometer) {
		if (nrhs < 3) mexErrMsgTxt("gait2dc: 3rd input argument must be xdot.");
		nrows = mxGetM(prhs[2]);
		ncols = mxGetN(prhs[2]);
		if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) mexErrMsgTxt("gait2d: Incorrect type for xdot, must be double.");
		if ((nrows != nstates) || (ncols != 1) ) mexErrMsgTxt("gait2dx: Incorrect size for xdot, must be a state column vector.");
		xdot = mxGetPr(prhs[2]);
	}

    // If command is "fkin", we can directly call the autolev function
    if (CMDFKin) {

        // Call autolev function
        q = &x[0];
        qd = &x[NDOF];	
        qdd = &xdot[NDOF]; // Might be able to remove qdd from the input in future.
        gait2dc_FK_al(&param, q, qd, qdd, fk, dfk_dq, fkdot, dfkdot_dq);


        // Assemble the MEX function outputs for the "Fkin" command
        int i;
        double *fk_matlab;

        plhs[0] = mxCreateDoubleMatrix(NFK, 1, mxREAL);
        fk_matlab = mxGetPr(plhs[0]);
        for (i=0; i<NFK; i++) {
            fk_matlab[i] = fk[i];
        }

        // output the Jacobian dFK/dq, if requested
        // this is a full matrix (maybe we will make it sparse later)
        if (derivatives) {
            int i,j;
            double *fkjac;
            plhs[1] = mxCreateDoubleMatrix(NFK, NDOF, mxREAL);
            fkjac = mxGetPr(plhs[1]);
            for (j=0; j<NDOF; j++) for (i=0; i<NFK; i++) *fkjac++ = dfk_dq[i][j];
        }

        // if a 3rd output was requested, also output dFKdot/dq
        if (nlhs > 2) {
            int i,j;
            double *fkdotjac;
            plhs[2] = mxCreateDoubleMatrix(NFK, NDOF, mxREAL);
            fkdotjac = mxGetPr(plhs[2]);
            for (j=0; j<NDOF; j++) for (i=0; i<NFK; i++) *fkdotjac++ = dfkdot_dq[i][j];
        }

        return;
    }

    // If command is "stick", we can directly call the autolev function
    if (CMDstick) {

        // Call autolev function
        q = &x[0];
        qd = &x[NDOF];	
        qdd = &xdot[NDOF]; // Might be able to remove qdd from the input in future.
        gait2dc_stick_al(&param, q, qd, qdd, fk, dfk_dq, fkdot, dfkdot_dq, Stick);


        // Assemble the MEX function outputs for the "stick" command
        // stick figure output are two matrix with two columns
		// right stick points and left stick points
		mwSize NRstick, NLstick;
		
		// We use 5 and 4 stick points generated by Autolev (remember that ankle is used twice)
		NRstick = 5;
		NLstick = 4;
		
		// find out how many contact points there are on each side
		for (i=0; i<ncontacts; i++) {
			if (contacts[i].foot == 0) NRstick++;
			if (contacts[i].foot == 1) NLstick++;
		}
		
		// generate the right side stick points
		plhs[0] = mxCreateDoubleMatrix(NRstick, 2, mxREAL);
		f = mxGetPr(plhs[0]);
		*f++ = Stick[0][0];				// x coordinate of trunk CM
		*f++ = Stick[1][0];				// x coordinate of hip
		*f++ = Stick[2][0];				// x coordinate of Rknee
		*f++ = Stick[3][0];				// x coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			k = 2*NDOF+2*NMUS+4*j;		// place in state vector x where X coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[3][0];				// x coordinate of Rankle
		
		*f++ = Stick[0][1];				// y coordinate of trunk CM
		*f++ = Stick[1][1];				// y coordinate of hip
		*f++ = Stick[2][1];				// y coordinate of Rknee
		*f++ = Stick[3][1];				// y coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			k = 2*NDOF+2*NMUS+4*j+1;		// place in state vector x where Y coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[3][1];				// y coordinate of Rankle
		
		// generate the left side stick points
		plhs[1] = mxCreateDoubleMatrix(NLstick, 2, mxREAL);
		f = mxGetPr(plhs[1]);
		*f++ = Stick[1][0];				// x coordinate of hip
		*f++ = Stick[4][0];				// x coordinate of Lknee
		*f++ = Stick[5][0];				// x coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			k = 2*NDOF+2*NMUS+4*j;		// place in state vector x where X coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[5][0];				// x coordinate of Lankle

		*f++ = Stick[1][1];				// y coordinate of hip
		*f++ = Stick[4][1];				// y coordinate of Lknee
		*f++ = Stick[5][1];				// y coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			k = 2*NDOF+2*NMUS+4*j+1;		// place in state vector x where Y coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[5][1];				// y coordinate of Lankle

		if (nlhs == 2) return;

		// If 3rd and 4th outputs were requested, output undeformed foot shapes
		
		// Count the points again (ankle, contact points, ankle)
		NRstick = 2;
		NLstick = 2;
		for (i=0; i<ncontacts; i++) {
			if (contacts[i].foot == 0) NRstick++;
			if (contacts[i].foot == 1) NLstick++;
		}
		
		// generate the right side stick points
		plhs[2] = mxCreateDoubleMatrix(NRstick, 2, mxREAL);
		f = mxGetPr(plhs[2]);
		*f++ = Stick[3][0];				// x coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			// calculate x coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[18] + fk[20]*contacts[j].x + fk[21]*contacts[j].y;
		}
		*f++ = Stick[3][0];				// x coordinate of Rankle
		
		*f++ = Stick[3][1];				// y coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			// calculate y coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[19] + fk[22]*contacts[j].x + fk[23]*contacts[j].y;
		}
		*f++ = Stick[3][1];				// y coordinate of Rankle
		
		// generate the left side stick points
		plhs[3] = mxCreateDoubleMatrix(NLstick, 2, mxREAL);
		f = mxGetPr(plhs[3]);
		*f++ = Stick[5][0];				// x coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			// calculate x coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[36] + fk[38]*contacts[j].x + fk[39]*contacts[j].y;
		}
		*f++ = Stick[5][0];				// x coordinate of Lankle

		*f++ = Stick[5][1];				// y coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			// calculate y coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[37] + fk[40]*contacts[j].x + fk[41]*contacts[j].y;
		}
		*f++ = Stick[5][1];				// y coordinate of Lankle
		
		return;

    }

	// If command is "Dynamics", we also need u and optionally M as inputs
	if (CMDdynamics) {
		nrows = mxGetM(prhs[3]);
		ncols = mxGetN(prhs[3]);
		if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
			mexErrMsgTxt("gait2dc: Incorrect type for u, must be double.");
		if ((nrows != NMUS) || (ncols != 1) )
			mexErrMsgTxt("gait2dc(Dynamics): Incorrect size for u, must be a 16 x 1 column vector.");
		u = mxGetPr(prhs[3]);
		
		if (nrhs == 5) {
			// additional joint moments were given as inputs
			nrows = mxGetM(prhs[4]);
			ncols = mxGetN(prhs[4]);
			if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[3]) )
				mexErrMsgTxt("gait2dc(Dynamics): Incorrect type for M, must be double.");
			if ((nrows != NMOM) || (ncols != 1) )
				mexErrMsgTxt("gait2dc(Dynamics): Incorrect size for M, must be a 6 x 1 column vector.");
			M = mxGetPr(prhs[4]);
		}
	}

	// Compute the muscle contraction dynamics, and get muscle forces
	for(i=0; i<NMUS; i++) {
		
		// Calculate muscle length Lm and derivatives dLm/dq from generalized coordinates in x
		Lm[i] = MusclePath(&muscles[i], &x[3], dLm_dang[i]);
	
		// Calculate muscle force imbalance, normalized to Fmax
		g[i] = MuscleDynamics(&muscles[i],
			x[2*NDOF+NMUS+i],		// active state of muscle i
			x[2*NDOF+i],			// Lce of muscle i
			xdot[2*NDOF+i],			// Lcedot
			Lm[i],					// muscle length
			&dg_da[i], &dg_dLce[i], &dg_dLcedot[i], &dg_dLm[i],
			&force[i], &dforce_dLm[i], &dforce_dLce[i],
            &forceCE[i], &dforceCE_da[i], &dforceCE_dLce[i], &dforceCE_dLcedot[i],
			&CEpower[i]);				
	}

	// Compute the joint moments
	for (i=0; i<NMOM; i++) {
		// initialize derivatives to zero
		if (derivatives) {
			for (j=0; j<NMOM; j++) dmom_dang[i][j] = 0.0;
			for (j=0; j<NMUS; j++) dmom_dLce[i][j] = 0.0;
			dmom_dangvel[i] = 0.0;							// this is a diagonal matrix so we only store diagonal
		}
		
		// start with passive joint moment
		ang = x[i+3];											// joint angle is one of the state variables
		angvel = x[NDOF+i+3];									// the corresponding angular velocity
		// calculate deformations (positive when outside ROM)
		dmax = Softpos(ang - param.MaxAngle[i], param.JointD[i], &ddmax);
		dmin = Softpos(param.MinAngle[i] - ang, param.JointD[i], &ddmin);

		// calculate moments due to lienar term & both endpoints and add them, also add damping
		mom[i] = -param.JointK1[i]*(ang-param.JointPhi0[i]) - param.JointK2[i] * (dmax - dmin) - param.JointB[i] * angvel;

		if (derivatives) {
			dmom_dang[i][i] = -param.JointK1[i] - param.JointK2[i] * (ddmax + ddmin);
			dmom_dangvel[i] = -param.JointB[i];
		}
		
		// add the muscle moments
		for (j=0; j<NMUS; j++) {
			mom[i] = mom[i] - dLm_dang[j][i]*force[j];		// moment arm is -dLm/dang
			if (derivatives) {
				for (k=0; k<NMOM; k++) {
					dmom_dang[i][k] = dmom_dang[i][k] - dLm_dang[j][i]*dforce_dLm[j]*dLm_dang[j][k];
				}
				dmom_dLce[i][j] = dmom_dLce[i][j] - dLm_dang[j][i]*dforce_dLce[j];
			}
		}
		
		// add the extra moments given as inputs
		if (CMDdynamics && nrhs == 5) {
			mom[i] = mom[i] + M[i];
		}
	}
	
	// Assemble the MEX function outputs for the "Jointmoments" command
	if (CMDjointmoments) {
		plhs[0] = mxCreateDoubleMatrix(NMOM, 1, mxREAL);
		f = mxGetPr(plhs[0]);
		// fill the column vector
		for (i=0; i<NMOM; i++) f[i] = mom[i];
		
		// if requested, return derivatives dM/dx, copy these from the dense matrices that we already have
		if (derivatives) {							
			Nalloc_dM_dx = 2*nonzeromomentarms + NMOM;			// is a bit too generous, some angle-angle combinations arise from more than one muscle
			plhs[1] = mxCreateSparse(nstates, NMOM, Nalloc_dM_dx, mxREAL);
			dM_dx = mxGetPr(plhs[1]);
			dM_dx_irs = mxGetIr(plhs[1]);
			dM_dx_jcs = mxGetJc(plhs[1]);
			NdM_dx = 0;
			for (i=0; i<NMOM; i++) {
				dM_dx_jcs[i] = NdM_dx;					// store element number where this column starts			

				// rows from dM/dang
				for (j=0; j<NMOM; j++) {				
					if (dmom_dang[i][j] != 0.0) {
						dM_dx_irs[NdM_dx] = j+3;					// store row number of this matrix element
						dM_dx[NdM_dx++] = dmom_dang[i][j];			// store the value of this matrix element
					}
				}
					
				// one row, on diagonal, from dM/dangvel
				dM_dx_irs[NdM_dx] = NDOF+i+3;				// store row number of this matrix element
				dM_dx[NdM_dx++] = dmom_dangvel[i];			// store the value of this matrix element
				
				// rows from dM/dLce
				for (j=0; j<NMUS; j++) {				
					if (dmom_dLce[i][j] != 0.0) {
						dM_dx_irs[NdM_dx] = 2*NDOF+j;		// store row number of this matrix element
						dM_dx[NdM_dx++] = dmom_dLce[i][j];			// store the value of this matrix element
					}
				}			
				
			}
			dM_dx_jcs[NMOM] = NdM_dx;						// store final element number

			if (Nalloc_dM_dx < NdM_dx) {
				printf("dM_dx: allocated %d, actual %d\n", Nalloc_dM_dx, NdM_dx);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for dM_dx.");			
			}
		}
		return;			

	}
	
	// Assemble the MEX function outputs for the "Muscleforces" command
	if (CMDmuscleforces) {
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		f = mxGetPr(plhs[0]);	
		// fill the column vector with muscle forces
		for (j=0; j<NMUS; j++) *f++ = force[j];
        
        // if requested, also return derivatives dF/dx, copy these from the dense matrices that we already have
		// muscle force (SEE force) depends only on joint angles and on Lce
		if (derivatives) {							
			Nalloc_dF_dx = 2*nonzeromomentarms + NMUS;	// is a bit too generous, but small enough
			plhs[1] = mxCreateSparse(nstates, NMUS, Nalloc_dF_dx, mxREAL);
			dF_dx = mxGetPr(plhs[1]);
			dF_dx_irs = mxGetIr(plhs[1]);
			dF_dx_jcs = mxGetJc(plhs[1]);
			NdF_dx = 0;
			for (i=0; i<NMUS; i++) {
				dF_dx_jcs[i] = NdF_dx;					// store element number where this column starts			

				// rows from dF/dang (for all NMOM joint angles)
				for (j=0; j<NMOM; j++) {
					double dF_dang = dforce_dLm[i]*dLm_dang[i][j];	// chain rule
					if (dF_dang != 0.0) {
						dF_dx_irs[NdF_dx] = j+3;			// store row number of this matrix element
						dF_dx[NdF_dx++] = dF_dang;			// store the value of this matrix element
					}
				}
					
				// one element for each muscle, on diagonal, from dF/dLce
				dF_dx_irs[NdF_dx] = 2*NDOF+i;				// store row number of this matrix element
				dF_dx[NdF_dx++] = dforce_dLce[i];			// store the value of this matrix element
								
			}
			dF_dx_jcs[NMUS] = NdF_dx;						// store final element number

			if (Nalloc_dF_dx < NdF_dx) {
				printf("dF_dx: allocated %d, actual %d\n", Nalloc_dF_dx, NdF_dx);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for dF_dx.");			
			}
        }
		return;
	}
    
    // Assemble the MEX function outputs for the "MuscleCEforces" command
	if (CMDmuscleCEforces) {
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		f = mxGetPr(plhs[0]);	
		// fill the column vector with muscle forces
		for (j=0; j<NMUS; j++) *f++ = forceCE[j];
        
        // if requested, also return derivatives dFCE/dx and dFCE/dxdot, copy these from the dense matrices that we already have
        if (derivatives) {
            Nalloc_dfCE_dx = 2*nonzeromomentarms + NMUS;	// is a bit too generous, but small enough
            Nalloc_dfCE_dxdot = 2*nonzeromomentarms + NMUS;	// is a bit too generous, but small enough
			plhs[1] = mxCreateSparse(nstates, NMUS, Nalloc_dfCE_dx, mxREAL);
            plhs[2] = mxCreateSparse(nstates, NMUS, Nalloc_dfCE_dxdot, mxREAL);
			dfCE_dx = mxGetPr(plhs[1]);
			dfCE_dx_irs = mxGetIr(plhs[1]);
			dfCE_dx_jcs = mxGetJc(plhs[1]);
            dfCE_dxdot = mxGetPr(plhs[2]);
			dfCE_dxdot_irs = mxGetIr(plhs[2]);
			dfCE_dxdot_jcs = mxGetJc(plhs[2]);
			NdfCE_dx = 0;
            NdfCE_dxdot = 0;
			for (i=0; i<NMUS; i++) {
				dfCE_dx_jcs[i] = NdfCE_dx;					     // store element number where this column starts			
                dfCE_dxdot_jcs[i] = NdfCE_dxdot;		         // store element number where this column starts			

                // one element for each muscle, on diagonal, from dF/dLce
				dfCE_dx_irs[NdfCE_dx] = 2*NDOF+i;				 // store row number of this matrix element
				dfCE_dx[NdfCE_dx++] = dforceCE_dLce[i];			 // store the value of this matrix element
                
                // one element for each muscle, on diagonal, from dF/dLcedot
				dfCE_dxdot_irs[NdfCE_dxdot] = 2*NDOF+i;			 // store row number of this matrix element
				dfCE_dxdot[NdfCE_dxdot++] = dforceCE_dLcedot[i]; // store the value of this matrix element
                
                // one element for each muscle, on diagonal, from dF/da
				dfCE_dx_irs[NdfCE_dx] = 2*NDOF+NMUS+i;	         // store row number of this matrix element
				dfCE_dx[NdfCE_dx++] = dforceCE_da[i];			 // store the value of this matrix element
                								
			}
			dfCE_dx_jcs[NMUS] = NdfCE_dx;						 // store final element number
            dfCE_dxdot_jcs[NMUS] = NdfCE_dxdot;					 // store final element number

			if (Nalloc_dfCE_dx < NdfCE_dx) {
				printf("dfCE_dx: allocated %d, actual %d\n", Nalloc_dfCE_dx, NdfCE_dx);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for dfCE_dx.");			
			}      
            if (Nalloc_dfCE_dxdot < NdfCE_dxdot) {
				printf("dfCE_dxdot: allocated %d, actual %d\n", Nalloc_dfCE_dxdot, NdfCE_dxdot);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for dFCE_dxdot.");			
			}   
        }
        return;
	}
	
	// Assemble the MEX function outputs for the "MuscleCEpower" command
	if (CMDmuscleCEpower) {
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		f = mxGetPr(plhs[0]);	
		f2 = mxGetPr(plhs[1]);	
		// fill the column vectors with CE power and contraction imbalance
		for (j=0; j<NMUS; j++) {
			*f++ = CEpower[j];
			*f2++ = g[j];
		}
		return;
	}
	
	// For those commands that need it, combine the point contact forces into total force and moment on each foot (grf)
	if (CMDgrf || CMDdynamics) {
		for (i=0; i<NGRF; i++) {
			grf[i] = 0.0;	// start with zero, then add up all contact points
			for (j=0; j<4*ncontacts; j++) {
				dgrfdx[i][j]   = 0.0;
			}
		}
		for (i=0; i<ncontacts; i++) {
			j = 2*NDOF+2*NMUS+4*i;			// index within x for the contact variables x,y,Fx,Fy for this point
			if (contacts[i].foot == 0) {	// right foot
				k = 0;						// index within grf where this point should be accumulated
			}
			else {
				k = 3;
			}
			grf[k]   += x[j+2];						// total Fx on foot
				dgrfdx[k][4*i]   += 0;
				dgrfdx[k][4*i+1] += 0;
				dgrfdx[k][4*i+2] += 1;
				dgrfdx[k][4*i+3] += 0;
			grf[k+1] += x[j+3];						// total Fy on foot
				dgrfdx[k+1][4*i]   += 0;
				dgrfdx[k+1][4*i+1] += 0;
				dgrfdx[k+1][4*i+2] += 0;
				dgrfdx[k+1][4*i+3] += 1;
			grf[k+2] += x[j]*x[j+3] - x[j+1]*x[j+2];	// total Mz, add x*Fy-y*Fz
				dgrfdx[k+2][4*i]   += x[j+3];
				dgrfdx[k+2][4*i+1] += -x[j+2];
				dgrfdx[k+2][4*i+2] += -x[j+1];
				dgrfdx[k+2][4*i+3] += x[j];		
		}
	}

	// Assemble the MEX function outputs for the "GRF" command
	if (CMDgrf) {	
		plhs[0] = mxCreateDoubleMatrix(NGRF, 1, mxREAL);
		f = mxGetPr(plhs[0]);
		for (i=0; i<NGRF; i++) f[i] = grf[i];		// copy from the grf array we already have
		
		// if requested, return derivatives dGRF/dx, copy these from the dense matrix dgrfdx that we already have
		if (derivatives) {
			Nalloc_dgrf_dx = 6*ncontacts;			// 6 derivatives for each contact point
			plhs[1] = mxCreateSparse(nstates, NGRF, Nalloc_dgrf_dx, mxREAL);
			dgrf_dx = mxGetPr(plhs[1]);
			dgrf_dx_irs = mxGetIr(plhs[1]);
			dgrf_dx_jcs = mxGetJc(plhs[1]);
			Ndgrf_dx = 0;
			for (i=0; i<NGRF; i++) {
				dgrf_dx_jcs[i] = Ndgrf_dx;						// store element number where this column starts			
				for (j=0; j<4*ncontacts; j++) {					// GRF only depends on contact states in x, rows 2*NDOF+2*NMUS+1..nstates
					if (dgrfdx[i][j] != 0.0) {
						dgrf_dx_irs[Ndgrf_dx] = 2*NDOF+2*NMUS+j;	// store row number of this matrix element
						dgrf_dx[Ndgrf_dx++] = dgrfdx[i][j];			// store the value of this matrix element
					}
				}
			}
			dgrf_dx_jcs[NGRF] = Ndgrf_dx;		// store final element number

			if (Nalloc_dgrf_dx < Ndgrf_dx) {
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for dgrf_dx.");			
			}
		}
				
		return;	
	}	
	
	// Call the C function that was generated by Autolev
	q = &x[0];
	qd = &x[NDOF];	
	qdd = &xdot[NDOF];
	gait2dc_dynamics_al(&param, q, qd, qdd, mom, grf, zero, dz_dq, dz_dqd, dz_dqdd, dz_dmom, dz_dgrf, fk, dfk_dq, fkdot, dfkdot_dq);
	
	// Assemble the MEX function outputs for the "Dynamics" command
	if (CMDdynamics) {
		plhs[0] = mxCreateDoubleMatrix(nstates, 1, mxREAL);
		f = mxGetPr(plhs[0]);

		// create MEX outputs for the sparse Jacobians
		if (derivatives) {
			// The sparse Jacobians have to be filled in column order, using Matlab sparse data structure
			// For this reason we will return the transpose of the Jacobians: (df/dx)' etc.

			// --------Jacobian df/dx
			Nalloc_df_dx = NDOF*(1 + 2*NDOF) - 52 + 3*NMUS + 2*nonzeromomentarms + 54*ncontacts;
			plhs[1] = mxCreateSparse(nstates, nstates, Nalloc_df_dx, mxREAL);
			df_dx = mxGetPr(plhs[1]);
			df_dx_irs = mxGetIr(plhs[1]);
			df_dx_jcs = mxGetJc(plhs[1]);
			Ndf_dx = 0;

			// --------Jacobian df/dxdot
			Nalloc_df_dxdot = NDOF*(NDOF+1) + 2*NMUS + 5*ncontacts;
			plhs[2] = mxCreateSparse(nstates, nstates, Nalloc_df_dxdot, mxREAL);
			df_dxdot = mxGetPr(plhs[2]);
			df_dxdot_irs = mxGetIr(plhs[2]);
			df_dxdot_jcs = mxGetJc(plhs[2]);
			Ndf_dxdot = 0;

			// --------Jacobian df/du
			Nalloc_df_du = NMUS;
			plhs[3] = mxCreateSparse(NMUS, nstates, Nalloc_df_du, mxREAL);
			df_du = mxGetPr(plhs[3]);
			df_du_irs = mxGetIr(plhs[3]);
			df_du_jcs = mxGetJc(plhs[3]);
			Ndf_du = 0;
			
			//----------Jacobian df/dM
			if (nlhs > 4) {
				Nalloc_df_dM = NDOF*NMOM;		// NOTE: actual number of nonzeros will be a lot smaller!
				plhs[4] = mxCreateSparse(NMOM, nstates, Nalloc_df_dM, mxREAL);  
				df_dM = mxGetPr(plhs[4]);
				df_dM_irs = mxGetIr(plhs[4]);
				df_dM_jcs = mxGetJc(plhs[4]);
				Ndf_dM = 0;
			}
		}
					
		// the first NDOF elements of the implicit differential equation are: qdot-dq/dt = 0
		for (i=0; i<NDOF; i++) {
			f[i] = x[NDOF+i] - xdot[i];
			if (derivatives) {
				df_dx_jcs[i]    = Ndf_dx;				// index for start of column i of (df/dx)'
				df_dxdot_jcs[i] = Ndf_dxdot;			// index for start of column i of (df/dxdot)'
				df_du_jcs[i]    = Ndf_du;				// index for start of column i of (df/du)'
				if (nlhs > 4) df_dM_jcs[i] = Ndf_dM;	// index for start of column i of (df/dM)'

				// ----------- df/dx
				df_dx_irs[Ndf_dx] = NDOF+i;			// store row number of this matrix element
				df_dx[Ndf_dx++] = 1.0;				// store the value of this matrix element

				// ----------- df/dxdot
				df_dxdot_irs[Ndf_dxdot] = i;		// store row number of this matrix element
				df_dxdot[Ndf_dxdot++] = -1.0;		// store the value of this matrix element
			}
		}
			
		// the next NDOF elements of the IDE are the equations of motion from Autolev (the ZERO expressions)
		for (i=0; i<NDOF; i++) {
			j = NDOF+i;
			f[j] = zero[i];
			if (derivatives) {
				df_dx_jcs[j]    = Ndf_dx;				// index for start of column j of (df/dx)'
				df_dxdot_jcs[j] = Ndf_dxdot;			// index for start of column j of (df/dxdot)'
				df_du_jcs[j]    = Ndf_du;				// index for start of column j of (df/du)'
				if (nlhs > 4) df_dM_jcs[j] = Ndf_dM;	// index for start of column j of (df/dM)'

				// ----------- df/dx
				// derivatives of ZERO[i] with respect to q go in rows 1..NDOF
				for (k=0; k<NDOF; k++) {
					double tmp;
					tmp = dz_dq[i][k];				// start with what came from Autolev
					// add the contributions dz/dmom * dmom/dq, these are zero unless q is an angle (k>2)
					if (k>2) for (kk=0; kk<NMOM; kk++) {
						tmp = tmp + dz_dmom[i][kk]*dmom_dang[kk][k-3];
					}
					if (fabs(tmp) > ALMOSTZERO) {
						df_dx_irs[Ndf_dx] = k;		// store row number of this matrix element
						df_dx[Ndf_dx++] = tmp;		// store the value of this matrix element
					}
				}
				
				// derivatives of ZERO[i] with respect to qd go in rows NDOF+1 to 2*NDOF
				for (k=0; k<NDOF; k++) {
					double tmp;		
					tmp = dz_dqd[i][k];			// start with what came from Autolev
					// add the contributions dz/dmom * dmom/dqdot, these are nonzero when qdot is an ang. vel (k>2)
					if (k>2) {
						tmp = tmp + dz_dmom[i][k-3]*dmom_dangvel[k-3];
					}
					if (fabs(tmp) > ALMOSTZERO) {
						df_dx_irs[Ndf_dx] = NDOF+k;		// store row number of this matrix element
						df_dx[Ndf_dx++] = tmp;			// store the value of this matrix element
					}
				}	
				
				// derivatives of ZERO[i] with respect to Lce should go into rows 2*NDOF+(1..NMUS)
				for (k=0; k<NMUS; k++) {
					double tmp;
					tmp = 0.0;
					for (kk=0; kk<NMOM; kk++) {
						// next line is a bit tricky, neg sign for moment arm, neg sign for dF/dLce cancel
						tmp = tmp + dz_dmom[i][kk]*dmom_dLce[kk][k];
					}
					if (fabs(tmp) > ALMOSTZERO) {
						df_dx_irs[Ndf_dx] = 2*NDOF+k;		// store row number of this matrix element
						df_dx[Ndf_dx++] = tmp;				// store the value of this matrix element
					}
				}
				
				// derivatives of ZERO[i] with respect to ground contact variables should go into rows 2*NDOF+2*MUS+1 to end
				for (k=0; k<4*ncontacts; k++) {
					double tmp;
					tmp = 0.0;
					// add the contributions dz/dgrf * dgrf/dx, these are zero unless x is a contact force variable
					for (kk=0; kk<NGRF; kk++) {
						tmp = tmp + dz_dgrf[i][kk]*dgrfdx[kk][k];
					}
					if (fabs(tmp) > ALMOSTZERO) {
						df_dx_irs[Ndf_dx] = 2*NDOF+2*NMUS+k;	// store row number of this matrix element
						df_dx[Ndf_dx++] = tmp;					// store the value of this matrix element
					}
				}				
				
				// ----------- df/dxdot
				// derivatives of ZERO[i] with respect to qdd should go into rows NDOF+1 to 2*NDOF of df/dxdot
				for (k=0; k<NDOF; k++) {
					if (fabs(dz_dqdd[i][k]) > ALMOSTZERO) {
						df_dxdot_irs[Ndf_dxdot] = NDOF+k;		// store row number of this matrix element
						df_dxdot[Ndf_dxdot++] = dz_dqdd[i][k];	// store the value of this matrix element
					}
				}
				
				// ------------ df/dM
				if (nlhs > 4) {
					for (k=0; k<NMOM; k++) {			
						// derivatives of Zero with respect to each moment, rows 1..NMOM
						if (fabs(dz_dmom[i][k]) > ALMOSTZERO) {
							df_dM_irs[Ndf_dM] = k;				// store row number of this matrix element
							df_dM[Ndf_dM++] = dz_dmom[i][k];		// store the value of this matrix element		
						}
					}
				}
			}
		}
		
		// the next NMUS elements of the IDE are the muscle contraction dynamics
		for (i=0; i<NMUS; i++) {
			j = 2*NDOF+i;
			f[j] = g[i];
			if (derivatives) {
				df_dx_jcs[j]    = Ndf_dx;				// index for start of column j of (df/dx)'
				df_dxdot_jcs[j] = Ndf_dxdot;			// index for start of column j of (df/dxdot)'
				df_du_jcs[j]    = Ndf_du;				// index for start of column j of (df/du)'
				if (nlhs > 4) df_dM_jcs[j] = Ndf_dM;	// index for start of column j of (df/dM)'

				// ----------- df/dx
				// derivatives of muscle imbalance with respect to q are in rows 1..NDOF
				for (k=0; k<NDOF; k++) {
					// element only exists if k is an angle and muscle i spans DOF k
					if (k>2) if (fabs(dLm_dang[i][k-3]) > ALMOSTZERO) {
						df_dx_irs[Ndf_dx] = k;							// store row number of this matrix element
						df_dx[Ndf_dx++] = dg_dLm[i]*dLm_dang[i][k-3];	// store the value of this matrix element
					}
				}

				// derivatives of muscle imbalance with respect to Lce are diagonal, rows 2*NDOF+1 to 2*NDOF+NMUS
				df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
				df_dx[Ndf_dx++] = dg_dLce[i];		// store the value of this matrix element

				// derivatives of muscle imbalance with respect to Act are diagonal, rows 2*NDOF+NMUS+(1..NMUS)
				df_dx_irs[Ndf_dx] = 2*NDOF+NMUS+i;	// store row number of this matrix element
				df_dx[Ndf_dx++] = dg_da[i];			// store the value of this matrix element
				
				// derivatives of muscle imbalance with respect to Lcedot are diagonal, rows 2*NDOF to 2*NDOF+NMUS
				df_dxdot_irs[Ndf_dxdot] = 2*NDOF+i;			// store row number of this matrix element
				df_dxdot[Ndf_dxdot++] = dg_dLcedot[i];		// store the value of this matrix element
			}			
		}
		
		// the next NMUS elements of the IDE are the muscle activation dynamics: da/dt - (u-a)(c1*u + c2) = 0
		for (i=0; i<NMUS; i++) {
			double rate = muscles[i].c1*u[i] + muscles[i].c2;
			j = 2*NDOF+NMUS+i;
			f[j] = xdot[j] - (u[i] - x[j])*rate;
			if (derivatives) {
				df_dx_jcs[j]    = Ndf_dx;				// index for start of column j of (df/dx)'
				df_dxdot_jcs[j] = Ndf_dxdot;			// index for start of column j of (df/dxdot)'
				df_du_jcs[j]    = Ndf_du;				// index for start of column j of (df/du)'
				if (nlhs > 4) df_dM_jcs[j] = Ndf_dM;	// index for start of column j of (df/dM)'

				// ----------- df/dx
				// derivatives of activation dynamics with respect to Act are diagonal, rows 2*NDOF+NMUS+(1..NMUS)
				df_dx_irs[Ndf_dx] = 2*NDOF+NMUS+i;		// store row number of this matrix element
				df_dx[Ndf_dx++] = rate;					// store the value of this matrix element

				// ----------- df/dxdot
				// derivatives of activation dynamics with respect to Actdot are diagonal, rows 2*NDOF+NMUS+(1..NMUS)
				df_dxdot_irs[Ndf_dxdot] = 2*NDOF+NMUS+i;		// store row number of this matrix element
				df_dxdot[Ndf_dxdot++] = 1.0;				// store the value of this matrix element

				// ----------- df/du
				// derivatives of activation dynamics with respect to u are diagonal, rows 1..NMUS
				df_du_irs[Ndf_du] = i;									// store row number of this matrix element
				df_du[Ndf_du] =  -rate - muscles[i].c1*(u[i] - x[j]);	// store the value of this matrix element
				Ndf_du++;				
			}
		}

		// then there are 4 equations for each contact point
		for (i=0; i<ncontacts; i++) {
			int ifkptr;							// pointer to forward kinematics results for foot
			double df_dfk[4][6];
			double df_dfkdot[4][6];
			double df_dxc[4][4];
			double df_dxcdot[4][2];
			int ixc;
			
			// index to contact variables x,y,Fx,Fy of this contact point
			ixc = 2*NDOF + 2*NMUS + 4*i;			

			// get pointer to forward kinematic results for this foot
			if (contacts[i].foot == 0) {
				ifkptr = 18;					//  pointer to forward kinematics variables fk[] for right foot
			} 
			else {
				ifkptr = 36;					//  pointer to forward kinematics variables fk[] for left foot
			}	

			// use the C code generated by contact.al
			contact_al(&contacts[i], &fk[ifkptr], &fkdot[ifkptr], &x[ixc], &xdot[ixc],
				&f[ixc], df_dfk, df_dfkdot, df_dxc, df_dxcdot);
					
			if (derivatives) {
				double tmp;
				for (j=0; j<4; j++) {
					df_dx_jcs[ixc+j]    = Ndf_dx;				// index for start of column j+2 of (df/dx)'
					df_dxdot_jcs[ixc+j] = Ndf_dxdot;			// index for start of column j+2 of (df/dxdot)'
					df_du_jcs[ixc+j]    = Ndf_du;				// index for start of column j+2 of (df/du)'
					if (nlhs > 4) df_dM_jcs[ixc+j] = Ndf_dM;	// index for start of column j+2 of (df/dM)'

					// compute derivatives of f[ixc+j] with respect to q
					for (k=0; k<NDOF; k++) {
						tmp = 0.0;
						for (kk=0; kk<6; kk++) {
							tmp = tmp + df_dfk[j][kk]*dfk_dq[ifkptr+kk][k] + df_dfkdot[j][kk]*dfkdot_dq[ifkptr+kk][k];
						}
						if (fabs(tmp) > ALMOSTZERO) {
							df_dx_irs[Ndf_dx] = k;				// row number of matrix element
							df_dx[Ndf_dx++] = tmp;   			// store the value of this matrix element
						}
					}

					// compute derivatives of f[ixc+j] with respect to qdot
					for (k=0; k<NDOF; k++) {
						tmp = 0.0;
						for (kk=0; kk<6; kk++) {
							tmp = tmp + df_dfkdot[j][kk]*dfk_dq[ifkptr+kk][k];		// dfkd_dqdot = dfk_dq!
						}
						if (fabs(tmp) > ALMOSTZERO) {
							df_dx_irs[Ndf_dx] = NDOF+k;			// row number of matrix element
							df_dx[Ndf_dx++] = tmp;   			// store the value of this matrix element
						}
					}
						
					// compute derivatives of f[ixc+j] with respect to contact variables in x
					for (k=0; k<4; k++) {
						tmp = df_dxc[j][k];
						if (fabs(tmp) > ALMOSTZERO) {
							df_dx_irs[Ndf_dx] = ixc+k;			// row number of matrix element
							df_dx[Ndf_dx++] = tmp;				// store the value of this matrix element
						}
					}
					// compute derivatives of f[ixc+j] with respect to contact point velocities in xdot
					for (k=0; k<2; k++) {
						tmp = df_dxcdot[j][k];
						if (fabs(tmp) > ALMOSTZERO) {
							df_dxdot_irs[Ndf_dxdot] = ixc+k;	// row number of matrix element
							df_dxdot[Ndf_dxdot++] = tmp;		// store the value of this matrix element
						}
					}
				}
			}

		}
		
		// store final element numbers, and check whether memory allocation was correct
		if (derivatives) {
			df_dx_jcs[nstates]    			= Ndf_dx;			
			df_dxdot_jcs[nstates] 			= Ndf_dxdot;		
			df_du_jcs[nstates]       		= Ndf_du;				
			if (nlhs > 4) df_dM_jcs[nstates]= Ndf_dM;
	
			// printf("dfdx: allocated %d, actual %d\n", Nalloc_df_dx, Ndf_dx);
			// printf("dfdxdot: allocated %d, actual %d\n", Nalloc_df_dxdot, Ndf_dxdot);
			// printf("dfdu: allocated %d, actual %d\n", Nalloc_df_du, Ndf_du);
			// if (nlhs > 4) printf("dfdM: allocated %d, actual %d\n", Nalloc_df_dM, Ndf_dM);
			
			if (Nalloc_df_dx  < Ndf_dx) {
				printf("dfdx: allocated %d, actual %d\n", Nalloc_df_dx, Ndf_dx);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for df_dx.");			
			}
			if (Nalloc_df_dxdot < Ndf_dxdot) {
				printf("dfdxdot: allocated %d, actual %d\n", Nalloc_df_dxdot, Ndf_dxdot);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for df_dxdot.");	
			}
			if (Nalloc_df_du < Ndf_du) {
				printf("dfdu: allocated %d, actual %d\n", Nalloc_df_du, Ndf_du);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for df_du.");			
			}
			if (nlhs > 4) if (Nalloc_df_dM < Ndf_dM) {
				printf("dfdM: allocated %d, actual %d\n", Nalloc_df_dM, Ndf_dM);
				mexErrMsgTxt("gait2dc: Incorrect memory allocation for df_dM.");			
			}
		}
		return;
	}

	// Assemble the MEX function outputs for the "Stick" command
	if (CMDstick) {
		// stick figure output are two matrix with two columns
		// right stick points and left stick points
		mwSize NRstick, NLstick;
		
		// We use 5 and 4 stick points generated by Autolev (remember that ankle is used twice)
		NRstick = 5;
		NLstick = 4;
		
		// find out how many contact points there are on each side
		for (i=0; i<ncontacts; i++) {
			if (contacts[i].foot == 0) NRstick++;
			if (contacts[i].foot == 1) NLstick++;
		}
		
		// generate the right side stick points
		plhs[0] = mxCreateDoubleMatrix(NRstick, 2, mxREAL);
		f = mxGetPr(plhs[0]);
		*f++ = Stick[0][0];				// x coordinate of trunk CM
		*f++ = Stick[1][0];				// x coordinate of hip
		*f++ = Stick[2][0];				// x coordinate of Rknee
		*f++ = Stick[3][0];				// x coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			k = 2*NDOF+2*NMUS+4*j;		// place in state vector x where X coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[3][0];				// x coordinate of Rankle
		
		*f++ = Stick[0][1];				// y coordinate of trunk CM
		*f++ = Stick[1][1];				// y coordinate of hip
		*f++ = Stick[2][1];				// y coordinate of Rknee
		*f++ = Stick[3][1];				// y coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			k = 2*NDOF+2*NMUS+4*j+1;		// place in state vector x where Y coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[3][1];				// y coordinate of Rankle
		
		// generate the left side stick points
		plhs[1] = mxCreateDoubleMatrix(NLstick, 2, mxREAL);
		f = mxGetPr(plhs[1]);
		*f++ = Stick[1][0];				// x coordinate of hip
		*f++ = Stick[4][0];				// x coordinate of Lknee
		*f++ = Stick[5][0];				// x coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			k = 2*NDOF+2*NMUS+4*j;		// place in state vector x where X coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[5][0];				// x coordinate of Lankle

		*f++ = Stick[1][1];				// y coordinate of hip
		*f++ = Stick[4][1];				// y coordinate of Lknee
		*f++ = Stick[5][1];				// y coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			k = 2*NDOF+2*NMUS+4*j+1;		// place in state vector x where Y coords of contact points are stored
			*f++ = x[k];				
		}
		*f++ = Stick[5][1];				// y coordinate of Lankle

		if (nlhs == 2) return;

		// If 3rd and 4th outputs were requested, output undeformed foot shapes
		
		// Count the points again (ankle, contact points, ankle)
		NRstick = 2;
		NLstick = 2;
		for (i=0; i<ncontacts; i++) {
			if (contacts[i].foot == 0) NRstick++;
			if (contacts[i].foot == 1) NLstick++;
		}
		
		// generate the right side stick points
		plhs[2] = mxCreateDoubleMatrix(NRstick, 2, mxREAL);
		f = mxGetPr(plhs[2]);
		*f++ = Stick[3][0];				// x coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			// calculate x coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[18] + fk[20]*contacts[j].x + fk[21]*contacts[j].y;
		}
		*f++ = Stick[3][0];				// x coordinate of Rankle
		
		*f++ = Stick[3][1];				// y coordinate of Rankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 0) {
			// calculate y coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[19] + fk[22]*contacts[j].x + fk[23]*contacts[j].y;
		}
		*f++ = Stick[3][1];				// y coordinate of Rankle
		
		// generate the left side stick points
		plhs[3] = mxCreateDoubleMatrix(NLstick, 2, mxREAL);
		f = mxGetPr(plhs[3]);
		*f++ = Stick[5][0];				// x coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			// calculate x coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[36] + fk[38]*contacts[j].x + fk[39]*contacts[j].y;
		}
		*f++ = Stick[5][0];				// x coordinate of Lankle

		*f++ = Stick[5][1];				// y coordinate of Lankle
		for (j=0; j<ncontacts; j++) if (contacts[j].foot == 1) {
			// calculate y coordinate of undeformed contact point using forward kinematics result
			*f++ = fk[37] + fk[40]*contacts[j].x + fk[41]*contacts[j].y;
		}
		*f++ = Stick[5][1];				// y coordinate of Lankle
		
		return;
	}
    
    // Assemble the MEX function outputs for the "Fkin" command
		if (CMDFKin) {
			int i;
			double *fk_matlab;
			
			plhs[0] = mxCreateDoubleMatrix(NFK, 1, mxREAL);
			fk_matlab = mxGetPr(plhs[0]);
			for (i=0; i<NFK; i++) {
				fk_matlab[i] = fk[i];
			}
			
			// output the Jacobian dFK/dq, if requested
			// this is a full matrix (maybe we will make it sparse later)
			if (derivatives) {
				int i,j;
				double *fkjac;
				plhs[1] = mxCreateDoubleMatrix(NFK, NDOF, mxREAL);
				fkjac = mxGetPr(plhs[1]);
				for (j=0; j<NDOF; j++) for (i=0; i<NFK; i++) *fkjac++ = dfk_dq[i][j];
			}
			
			// if a 3rd output was requested, also output dFKdot/dq
			if (nlhs > 2) {
				int i,j;
				double *fkdotjac;
				plhs[2] = mxCreateDoubleMatrix(NFK, NDOF, mxREAL);
				fkdotjac = mxGetPr(plhs[2]);
				for (j=0; j<NDOF; j++) for (i=0; i<NFK; i++) *fkdotjac++ = dfkdot_dq[i][j];
			}
			
			return;
		}
			
    mexErrMsgTxt("gait2dc: Programming error.");

}
