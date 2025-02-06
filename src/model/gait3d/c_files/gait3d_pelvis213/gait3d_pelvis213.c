/*=================================================================
 *
 * gait3d_pelvis213.c
 *
 * Implicit differential equation for 3D musculoskeletal model : f(x,dx/dt,u) = 0

 * This is the source code for the MEX function gait3d_pelvis213.mexw32
 * The musculoskeletal model is documented in the file gait3d_reference.docx
 * The function API documentation is in gait3d.m
 *
 * Copyright 2013-2017 Orchard Kinetics LLC
 *
 *=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mex.h"
#include "gait3d_pelvis213.h"

// set this to true when we're determining size of the Jacobians
#define JACOBIANSIZING 0

// maximum size of the model (some other size constants are in gait3d_pelvis213.h)
#define MAXCONTACTS 200				// maximum number of contact points in the model
#define MAXMUSCLES 200				// maximum number of muscles in the model
#define MAXSTATES (2*NDOF+2*MAXMUSCLES+NCVAR*MAXCONTACTS)		// max number of system state variables, including contact variables

// macro to extract all mass properties of a segment and store them as scalars in the parameters struct
// S: segment name, see gait3d_pelvis213.h.  NUM: index for the segments array that came out of Opensim
#define EXTRACTINERTIAL(S,NUM) {										\
	parameters.S##_M 	= extract2(P,"segments",(NUM),"mass",1);		\
	parameters.S##_CMx 	= extract2(P,"segments",(NUM),"mass_center",1);	\
	parameters.S##_CMy 	= extract2(P,"segments",(NUM),"mass_center",2);	\
	parameters.S##_CMz 	= extract2(P,"segments",(NUM),"mass_center",3);	\
	parameters.S##_Ixx 	= extract2(P,"segments",(NUM),"inertia",1);		\
	parameters.S##_Ixy 	= extract2(P,"segments",(NUM),"inertia",2);		\
	parameters.S##_Ixz 	= extract2(P,"segments",(NUM),"inertia",3);		\
	parameters.S##_Iyy 	= extract2(P,"segments",(NUM),"inertia",5);		\
	parameters.S##_Iyz 	= extract2(P,"segments",(NUM),"inertia",6);		\
	parameters.S##_Izz 	= extract2(P,"segments",(NUM),"inertia",9);		\
}

// static global variables are used to preserve model parameters betweem MEX function calls
static int initialized = 0;					// this will be set to 1 when initialization is complete
static param_struct parameters;				// contains model parameters that must be made available to Autolev multibody model
static muscleprop muscles[MAXMUSCLES];		// contains muscle parameters
static contactprop contacts[MAXCONTACTS];	// contains contact element parameters
static jointprop joints[NDOF];				// joint moment parameters
static int nmuscles, ncontacts, nstates, nf, nstrainterms;		// model size parameters.  nf is number of elements in the model expressions f(x,xdot,u,M)
static double zeros[MAXSTATES];				// an array of zeros
static double straincoeff[MAXPOLTERMS];			// coefficients of strain energy polynomial
static int strainexpon[MAXPOLTERMS][NDOF];		// exponents of strain energy polynomial

//===================================================================================
// MusclePaths:
// returns the lengths Lm, gradients dLm/dq, and Hessians d2Lm/dq2 for all muscles
//===================================================================================
void MusclePaths(
	double q[NDOF], 									// input: joint angles
	double Lm[MAXMUSCLES], 								// output: muscle-tendon lengths Lm
	double dLm_dq[MAXMUSCLES][NDOF], 					// output: first derivatives of Lm w.r.t. q
	double dLm_dqdq[MAXMUSCLES][NDOF][NDOF],			// output: second derivatives of Lm w.r.t. q
	int derivatives) {									// input: 1 if second derivatives should be calculated, 0 otherwise

	// local variables
	#define MAXPWR 4
	double qp[NDOF][MAXPWR+1], qp1[NDOF][MAXPWR+1], qp2[NDOF][MAXPWR+1];
	int i,j,k,k1,k2,p;
	muscleprop *m;
	
	// compute the powers of the generalized coordinates q, and 1st and 2nd derivatives
	// do a linear extrapolation outside of qmin_muscleMoment and qmax_muscleMoment
	for (i=6; i<NDOF; i++) {
		double qmin_muscleMoment = joints[i].qmin_muscleMoment;
		double qmax_muscleMoment = joints[i].qmax_muscleMoment;
		
		// q^0 and q^1 
		qp[i][0] = 1;
		qp1[i][0] = 0;
		qp2[i][0] = 0;

		qp[i][1] = q[i];
		qp1[i][1] = 1;
		qp2[i][1] = 0;

		// q^2, q^3 and q^4
		if (q[i] < qmin_muscleMoment) {
			// linear extrapolation below qmin_muscleMoment
			qp[i][2]  = qmin_muscleMoment * (2*q[i] - qmin_muscleMoment);
			qp1[i][2] = 2*qmin_muscleMoment;
			qp2[i][2] = 0;
			
			qp[i][3] = qmin_muscleMoment*qmin_muscleMoment * (3*q[i] - 2*qmin_muscleMoment);
			qp1[i][3] = 3*qmin_muscleMoment*qmin_muscleMoment;
			qp2[i][3] = 0;

			qp[i][4] = qmin_muscleMoment*qmin_muscleMoment*qmin_muscleMoment * (4*q[i] - 3*qmin_muscleMoment);;
			qp1[i][4] = 4*qmin_muscleMoment*qmin_muscleMoment*qmin_muscleMoment;
			qp2[i][4] = 0;
		}
		else if (q[i] < qmax_muscleMoment) {
			// q is between qmin_muscleMoment and qmax_muscleMoment
			qp[i][2] = q[i]*q[i];
			qp1[i][2] = 2*q[i];
			qp2[i][2] = 2;

			qp[i][3] = q[i]*q[i]*q[i];
			qp1[i][3] = 3*q[i]*q[i];
			qp2[i][3] = 6*q[i];

			qp[i][4] = q[i]*q[i]*q[i]*q[i];
			qp1[i][4] = 4*q[i]*q[i]*q[i];
			qp2[i][4] = 12*q[i]*q[i];
		}
		else {
			// linear extrapolation above qmax_muscleMoment
			qp[i][2]  = qmax_muscleMoment * (2*q[i] - qmax_muscleMoment);
			qp1[i][2] = 2*qmax_muscleMoment;
			qp2[i][2] = 0;

			qp[i][3] = qmax_muscleMoment*qmax_muscleMoment * (3*q[i] - 2*qmax_muscleMoment);
			qp1[i][3] = 3*qmax_muscleMoment*qmax_muscleMoment;
			qp2[i][3] = 0;

			qp[i][4] = qmax_muscleMoment*qmax_muscleMoment*qmax_muscleMoment * (4*q[i] - 3*qmax_muscleMoment);;
			qp1[i][4] = 4*qmax_muscleMoment*qmax_muscleMoment*qmax_muscleMoment;
			qp2[i][4] = 0;
		}

	}
	
	// loop through all muscles
	for (i=0; i<nmuscles; i++) {
		// initialize length and all derivatives to 0 
		Lm[i] = 0.0;
		for (j=0; j<NDOF; j++) {
			dLm_dq[i][j] = 0.0;    
			if (derivatives) {
				for (k=0; k<NDOF; k++) {	
					dLm_dqdq[i][j][k] = 0.0;
				}
			}
		}

		// add contributions from each polynomial term
		m = &muscles[i];
		for (j=0; j < m->npolterms; j++) {

			// add this term's contribution to Lm
			double term = m->polcoeff[j];
			for (k=0; k < m->nmusdof; k++) {
				p = m->expon[j][k];
				if (p > MAXPWR) {
					mexErrMsgTxt("gait3d_pelvis213: Polynomial power too high.");				
				}
				term = term * qp[m->musdof[k]][p];
			}
			Lm[i] = Lm[i] + term;	

			// add this term's contributions to the elements of dLm/dq
			for (k1=0; k1 < m->nmusdof; k1++) {
				int dof1 = m->musdof[k1];
				term = m->polcoeff[j] * qp1[dof1][m->expon[j][k1]];
				if (term != 0.0) {
					for (k=0; k < m->nmusdof; k++) {
						if (k != k1) term = term * qp[m->musdof[k]][m->expon[j][k]];
					}
					dLm_dq[i][dof1] = dLm_dq[i][dof1] + term;
				}
			}
		
			// add this term's contributions to the elements of the Hessian
			if (derivatives) {						// Hessian is only needed if Dynamics derivatives were requested
				for (k1=0; k1 < m->nmusdof; k1++) {
					int dof1 = m->musdof[k1];
					// diagonal elements (second derivative with respect to q(k1) )
					term = m->polcoeff[j] * qp2[dof1][m->expon[j][k1]];
					if (term != 0.0) {
						for (k=0; k < m->nmusdof; k++) {
							if (k != k1) term = term * qp[m->musdof[k]][m->expon[j][k]];
						}
						dLm_dqdq[i][dof1][dof1] = dLm_dqdq[i][dof1][dof1] + term;
						// printf("muscle %d term %d hessian element %d %d: %f\n", i, j, k1,k1, term);
					}
				
					// below-diagonal elements (row k1, column k2)
					for(k2=0; k2 < k1; k2++) {
						int dof2 = m->musdof[k2];
						term = m->polcoeff[j] * qp1[dof1][m->expon[j][k1]] * qp1[dof2][m->expon[j][k2]];
						if (term != 0.0) {
							for (k=0; k < m->nmusdof; k++) {
								if (k != k1) if (k != k2) term = term * qp[m->musdof[k]][m->expon[j][k]];
							}
							dLm_dqdq[i][dof1][dof2] = dLm_dqdq[i][dof1][dof2] + term;
							// printf("muscle %d term %d hessian element %d %d: %f\n", i, j, k1,k2, term);
							// above-diagonal elements can just be copied because Hessian is symmetric	
							dLm_dqdq[i][dof2][dof1] = dLm_dqdq[i][dof1][dof2];
						}
					}
				}
			}
		}
    }
}

//===================================================================================================
// MuscleDynamics: the implicit muscle contraction dynamics f()=0, returns f(a,s,sdot,Lm) and its derivatives
// For pennated muscle model, the state variable s is Lce*cos(p)
// state variable s is dimensionless (normalized to Lceopt)
// imbalance f is dimensionless (normalized to Fmax)
//===================================================================================================
double MuscleDynamics(muscleprop *m,									// muscle parameters (input)
	double a, double s, double sdot, double Lm,							// the input variables
	double *df_da, double *df_ds, double *df_dsdot, double *df_dLm, 	// the gradients (output)													
	double *force, double *dforce_dLm, double *dforce_ds,				// muscle force (N) and derivatives (output)
    double *forceCE, double *dforceCE_da, double *dforceCE_ds, double *dforceCE_dsdot, // muscle force of CE (N) and derivatives (output)
	double *powerCE) {													// CE power (W) (output)
	
	double Lce, cosp, dLce_ds, dcosp_ds, b;
	double Lcedot, dLcedot_dsdot, dLcedot_ds;
	double F1, dF1_dLce, dF1_ds;
	double f,x,k1;
	double F2, dF2_dLcedot, dF2_dsdot, dF2_ds, dF2_da;
	double F3, dF3_dLce, dF3_ds;
	double F4, dF4_ds, dF4_dLm;
	double F5, dF5_dsdot;
	double lambda, dlambda_da;		// factor for activation-dependent maximal shortening velocity
	double c;						// continuity parameter for eccentric force-velocity relationship
	
	// If there is pennation, compute Lce and cos(p) from state s using the constant volume constraint: Lce*sin(p) = Lceopt*sin(popt)
	// Lce is dimensionless (normalized to Lceopt)
	// If pennation is zero, we can't do this because volume is zero, and Lce is equal to s
	if (m->pennopt < 0.01) {
		cosp = 1;
		Lce = s;
		dLce_ds = 1;
		dcosp_ds = 0;
		}
	else {
		b = sin(m->pennopt);			
		Lce = sqrt(s*s + b*b);
		cosp = s/Lce;
		dLce_ds = cosp;
		dcosp_ds = b*b/Lce/Lce/Lce;
		}
	// Compute Lcedot and its derivatives with respect to sdot and s
	Lcedot = sdot*cosp;
	dLcedot_dsdot = cosp;
	dLcedot_ds = sdot*dcosp_ds;
		
	// F1 is the normalized isometric force-length relationship at max activation
	x = (Lce - 1.0)/m->width;
	F1 = exp(-x*x);
	dF1_dLce = -2.0*x*F1 / m->width;
	dF1_ds = dF1_dLce * dLce_ds;
		
	// F2 is the normalized force-velocity relationship
	lambda = 0.5025 + 0.5341*a;				// lambda is the factor for activation dependence of Vmax
	dlambda_da = 0.5341;
	// lambda = 1.0;
	// dlambda_da = 0.0;
	if (Lcedot < 0) {
		// concentric contraction
		x = lambda * m->Vmax - Lcedot/m->HillA;
		F2 = (lambda * m->Vmax + Lcedot)/x;
		dF2_dLcedot = (1.0 + F2/m->HillA)/x;
		dF2_da = -dlambda_da * m->Vmax * Lcedot * (1.0 + 1.0/m->HillA) / x / x;
	}
	else {
		// eccentric contraction
		c = lambda * m->c3;
		x = Lcedot + c;
		F2 = (m->gmax*Lcedot + c) / x;
		dF2_dLcedot = (m->gmax - F2)/x;
		dF2_da = dlambda_da * m->c3 * Lcedot * (1.0 - m->gmax)/x/x;
	}
	dF2_dsdot =  dF2_dLcedot * dLcedot_dsdot;
	dF2_ds = dF2_dLcedot * dLcedot_ds;
	
	// F3 is the PEE force-length relationship
	k1 = 1.0/m->Fmax*m->Lceopt;	// stiffness of the linear term is 1 N/m, convert to Fmax/Lceopt units	
	x = (Lce - m->PEEslack);		// elongation of PEE, relative to Lceopt
	F3 = k1*x;						// low stiffness linear term
	dF3_dLce = k1;
	if (x>0) {						// add quadratic term for positive elongation						
		F3 = F3 + m->kPEE*x*x;
		dF3_dLce = dF3_dLce + 2*m->kPEE*x;
	}
	dF3_ds = dF3_dLce * dLce_ds;
	
	//  F4 is the SEE force-length relationship
	k1 = 1.0/m->Fmax;			// stiffness of the linear term is 1 N/m, convert to Fmax/meter	
	x = Lm - s * m->Lceopt - m->SEEslack;			// elongation of SEE, in meters
	F4 = k1*x;										// low stiffness linear term
	dF4_ds = -k1*m->Lceopt;
	dF4_dLm = k1;
	if (x>0) {										// add quadratic term for positive deformation
		F4 = F4 + m->kSEE*x*x;
		dF4_ds = dF4_ds - 2 * m->kSEE * m->Lceopt * x;
		dF4_dLm = dF4_dLm + 2 * m->kSEE * x;
	}

	// F5 is viscous damping in the projected CE (0.001 of Fmax at 1 Lceopt/s) to ensure df/dLcedot is never zero
	// this is only really needed if we want to solve Lcedot explicitly for all possible muscle states (including a=0)
	F5 = .001*sdot;
	dF5_dsdot = .001;
		
	// Compute f, the force imbalance in the muscle contraction, and its derivatives
	f = F4 - (a*F1*F2 + F3)*cosp - F5;
	*df_da = -(F1*F2 + a*F1*dF2_da)*cosp;
	*df_ds = dF4_ds - (a*(dF1_ds*F2 + F1*dF2_ds) + dF3_ds)*cosp - (a*F1*F2 + F3)*dcosp_ds;
	*df_dsdot = -a*F1*dF2_dsdot*cosp - dF5_dsdot;
	*df_dLm = dF4_dLm;
	
	// Muscle force is the force in SEE
	*force = m->Fmax*F4;
	*dforce_dLm = m->Fmax * dF4_dLm;
	*dforce_ds = m->Fmax * dF4_ds;
    
    // Muscle force is the force in CE
	*forceCE = m->Fmax * (a*F1*F2*cosp + F5);
	*dforceCE_da = m->Fmax * (F1*F2*cosp+a*F1*dF2_da*cosp);
    *dforceCE_ds = m-> Fmax * a * ((dF1_ds*F2 + F1*dF2_ds)*cosp +F1*F2*dcosp_ds);
    *dforceCE_dsdot = m->Fmax * (a*F1*dF2_dsdot*cosp + dF5_dsdot);
    
	// power (W) generated by CE (positive when shortening)
	// Fce was in Fmax units, and Lce was in Lceopt units, so some scaling needed to get power in Watts
	*powerCE = -m->Fmax * a*F1*F2 * Lcedot * m->Lceopt;
	
	// printf("\n%s\n", m->name);
	// printf("Lm %f\n", Lm);
	// printf("s %f\n", s);
	// printf("Fsee %f\n", F4);
	// printf("Fce %f\n", a*F1*F2);
	// printf("Fpee %f\n", F3);

	// Return the imbalance
	return f;
}

//===================================================================================================
// normalize: normalize a vector to length 1, warn if length was not close enough to 1 to begin with
//===================================================================================================
void normalize(double *x, double *y, double *z) {
	double length;
	length = sqrt((*x)*(*x) + (*y)*(*y) + (*z)*(*z));
	if (fabs(length - 1.0) > 1e-5) {
		printf("gait3d_pelvis213: warning: vector was not normalized: %f %f %f -> length %f\n", *x, *y, *z, length);		
	}
	*x = *x / length;
	*y = *y / length;
	*z = *z / length;
}

//=========================================================================
// extract: returns value of P.fieldname1(i1) (double)
//=========================================================================
double extract(const mxArray *P, char *fieldname1, unsigned int i1) {
	
	mxArray *field;

	if ( ((field = mxGetField(P,0,fieldname1))==NULL) || 
	    (mxGetNumberOfDimensions(field)!=2) ||
		!mxIsDouble(field)||mxIsComplex(field) ) {
			printf("field name: %s\n", fieldname1);
			mexErrMsgTxt("gait3d_pelvis213: Valid field with this name not found in parameters structure.");
	}
	return mxGetPr(field)[i1-1];
}

//=========================================================================
// extract2: returns value of P.fieldname1{i1}.fieldname2(i2) (double) 
//=========================================================================
double extract2(const mxArray *P, char *fieldname1, unsigned int i1, char *fieldname2, unsigned int i2) {
	
	mxArray *field;

	if ( ((field = mxGetField(P,0,fieldname1))==NULL) || 
	    (mxGetNumberOfDimensions(field)!=2) || 
		!mxIsCell(field)||(i1==0)) {
			printf("field name: %s\n", fieldname1);
			mexErrMsgTxt("gait3d_pelvis213: Initialize: Valid field with this name not found in parameters structure.");
	}
    
    field = mxGetCell(field, i1-1);
    
	return extract(field, fieldname2, i2);
}

// =========================================================================
// mexFunction: this is the actual MEX function interface
// =========================================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	char *command;				// command string given as first input to the MEX function

	// get the first argument, make sure that it is a command string
	if (nrhs < 1)
		mexErrMsgTxt("gait3d_pelvis213: At least one input argument is needed.");
    if (!mxIsChar(prhs[0])) 
		mexErrMsgTxt("gait3d_pelvis213: First argument must be a string.");
 	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("gait3d_pelvis213: First argument must be a row vector.");
	command = mxArrayToString(prhs[0]);
	if(command == NULL) 
	  mexErrMsgTxt("gait3d_pelvis213: First argument was not a command string.");
	  
	// initialize?
	if (strcmp(command, "Initialize") == 0) {
		const mxArray *P;									// pointer to the parameter struct
		printf("***************************************************\n");
		printf("*          MEX Function gait3d_pelvis213          *\n");  
		printf("*       (c) 2013-2015 Orchard Kinetics LLC        *\n");
		printf("***************************************************\n");
		printf("Initializing...\n");
		if (nrhs != 2)
			mexErrMsgTxt("gait3d_pelvis213:Initialize: Two input arguments required.");
			
		// get a pointer P to the second argument and check that it is a 1x1 struct
		P = prhs[1];
		if ((mxGetM(P) != 1) || (mxGetN(P) != 1)) {
			mexErrMsgTxt("gait3d_pelvis213:Dynamics: Incorrect size for parameters, must be 1 x 1.");
		}
		if (!mxIsStruct(P)) {
			mexErrMsgTxt("gait3d_pelvis213:Initialize: Incorrect type for parameters, must be struct.");	
		}
		
		// from P, extract the parameters needed by the Autolev generated code, and store them in the C struct "parameters"
		// see gait3d_pelvis213.h for information about the fields in the parameters.struct
		parameters.gravity_x = extract(P, "gravity", 1);
		parameters.gravity_y = extract(P, "gravity", 2);
		parameters.gravity_z = extract(P, "gravity", 3);
		
		// parameters of the air drag model
		parameters.airdrag 	= extract(P, "drag_coefficient", 1);
		parameters.wind 	= extract(P, "wind_speed", 1);
		
		// Check that the nDofs in the model struct is the same as the constant NDOF we use internally here
		if (NDOF - extract(P, "nDofs", 1) != 0) {
			mexErrMsgTxt("gait3d_pelvis213:Initialize: nDofs is not consistent with this version of the MEX function.");		
		}
		
		// Check that the Nsegments in the model struct is one plus the NFKSEG we use internally here (Opensim has ground as a segment)
		if (NFKSEG + 1 - extract(P, "nSegments", 1) != 0) {
			mexErrMsgTxt("gait3d_pelvis213:Initialize: nSegments is not consistent with this version of the MEX function.");		
		}
		
		// from P, extract the joint position parameters
		{
		// joints are referred to by their index in the joints{} cell array
		// TO DO: better use joint names, in case we insert joints or change their order in Opensim
		parameters.back_x 		= extract2(P,"joints",12,"location",1);		
		parameters.back_y 		= extract2(P,"joints",12,"location",2);		

		parameters.Rhip_x 		= extract2(P,"joints",2,"location",1);		
		parameters.Rhip_y 		= extract2(P,"joints",2,"location",2);		
		parameters.Rhip_z 		= extract2(P,"joints",2,"location",3);		
		parameters.Rknee_x1 	= extract2(P,"joints",3,"t1_coefs",1);		
		parameters.Rknee_x2 	= extract2(P,"joints",3,"t1_coefs",2);		
		parameters.Rknee_x3 	= extract2(P,"joints",3,"t1_coefs",3);		
		parameters.Rknee_x4 	= extract2(P,"joints",3,"t1_coefs",4);		
		parameters.Rknee_x5 	= extract2(P,"joints",3,"t1_coefs",5) + extract2(P,"joints",3,"location",1);		
		parameters.Rknee_y1 	= extract2(P,"joints",3,"t2_coefs",1);		
		parameters.Rknee_y2 	= extract2(P,"joints",3,"t2_coefs",2);		
		parameters.Rknee_y3 	= extract2(P,"joints",3,"t2_coefs",3);		
		parameters.Rknee_y4 	= extract2(P,"joints",3,"t2_coefs",4);		
		parameters.Rknee_y5 	= extract2(P,"joints",3,"t2_coefs",5) + extract2(P,"joints",3,"location",2);	
		parameters.Rankle_y		= extract2(P,"joints",4,"location",2);		
		parameters.Rsubtalar_x	= extract2(P,"joints",5,"location",1);		
		parameters.Rsubtalar_y	= extract2(P,"joints",5,"location",2);		
		parameters.Rsubtalar_z	= extract2(P,"joints",5,"location",3);		
		parameters.Rmtp_x		= extract2(P,"joints",6,"location",1);		
		parameters.Rmtp_y		= extract2(P,"joints",6,"location",2);		
		parameters.Rmtp_z		= extract2(P,"joints",6,"location",3);		
		parameters.Racromial_x	= extract2(P,"joints",13,"location",1);		
		parameters.Racromial_y	= extract2(P,"joints",13,"location",2);		
		parameters.Racromial_z	= extract2(P,"joints",13,"location",3);		
		parameters.Relbow_x		= extract2(P,"joints",14,"location",1);		
		parameters.Relbow_y		= extract2(P,"joints",14,"location",2);		
		parameters.Relbow_z		= extract2(P,"joints",14,"location",3);		
		parameters.Rraduln_x	= extract2(P,"joints",15,"location",1);		
		parameters.Rraduln_y	= extract2(P,"joints",15,"location",2);		
		parameters.Rraduln_z	= extract2(P,"joints",15,"location",3);		
		parameters.Rwrist_x		= extract2(P,"joints",16,"location",1);		
		parameters.Rwrist_y		= extract2(P,"joints",16,"location",2);		
		parameters.Rwrist_z		= extract2(P,"joints",16,"location",3);		

		parameters.Lhip_x 		= extract2(P,"joints",7,"location",1);		
		parameters.Lhip_y 		= extract2(P,"joints",7,"location",2);		
		parameters.Lhip_z 		= extract2(P,"joints",7,"location",3);		
		parameters.Lknee_x1 	= extract2(P,"joints",8,"t1_coefs",1);		
		parameters.Lknee_x2 	= extract2(P,"joints",8,"t1_coefs",2);		
		parameters.Lknee_x3 	= extract2(P,"joints",8,"t1_coefs",3);		
		parameters.Lknee_x4 	= extract2(P,"joints",8,"t1_coefs",4);		
		parameters.Lknee_x5 	= extract2(P,"joints",8,"t1_coefs",5) + extract2(P,"joints",8,"location",1);		
		parameters.Lknee_y1 	= extract2(P,"joints",8,"t2_coefs",1);		
		parameters.Lknee_y2 	= extract2(P,"joints",8,"t2_coefs",2);		
		parameters.Lknee_y3 	= extract2(P,"joints",8,"t2_coefs",3);		
		parameters.Lknee_y4 	= extract2(P,"joints",8,"t2_coefs",4);		
		parameters.Lknee_y5 	= extract2(P,"joints",8,"t2_coefs",5) + extract2(P,"joints",8,"location",2);		
		parameters.Lankle_y		= extract2(P,"joints",9,"location",2);		
		parameters.Lsubtalar_x	= extract2(P,"joints",10,"location",1);		
		parameters.Lsubtalar_y	= extract2(P,"joints",10,"location",2);		
		parameters.Lsubtalar_z	= extract2(P,"joints",10,"location",3);		
		parameters.Lmtp_x		= extract2(P,"joints",11,"location",1);		
		parameters.Lmtp_y		= extract2(P,"joints",11,"location",2);		
		parameters.Lmtp_z		= extract2(P,"joints",11,"location",3);		
		parameters.Lacromial_x	= extract2(P,"joints",17,"location",1);		
		parameters.Lacromial_y	= extract2(P,"joints",17,"location",2);		
		parameters.Lacromial_z	= extract2(P,"joints",17,"location",3);		
		parameters.Lelbow_x		= extract2(P,"joints",18,"location",1);		
		parameters.Lelbow_y		= extract2(P,"joints",18,"location",2);		
		parameters.Lelbow_z		= extract2(P,"joints",18,"location",3);		
		parameters.Lraduln_x	= extract2(P,"joints",19,"location",1);		
		parameters.Lraduln_y	= extract2(P,"joints",19,"location",2);		
		parameters.Lraduln_z	= extract2(P,"joints",19,"location",3);		
		parameters.Lwrist_x		= extract2(P,"joints",20,"location",1);		
		parameters.Lwrist_y		= extract2(P,"joints",20,"location",2);		
		parameters.Lwrist_z		= extract2(P,"joints",20,"location",3);	
		}
		
		// From P, extract joint axis orientations.
		{
		double zero = 0.0;
		parameters.Rankle1 		= extract2(P,"joints",4,"r1_axis",1);		
		parameters.Rankle2 		= extract2(P,"joints",4,"r1_axis",2);		
		parameters.Rankle3 		= extract2(P,"joints",4,"r1_axis",3);		
		parameters.Rsubtalar1 	= extract2(P,"joints",5,"r1_axis",1);		
		parameters.Rsubtalar2 	= extract2(P,"joints",5,"r1_axis",2);		
		parameters.Rsubtalar3 	= extract2(P,"joints",5,"r1_axis",3);		
		parameters.Rmtp1 		= extract2(P,"joints",6,"r1_axis",1);		
		parameters.Rmtp3 		= extract2(P,"joints",6,"r1_axis",3);		
		parameters.Relbow1 		= extract2(P,"joints",14,"r1_axis",1);		
		parameters.Relbow2 		= extract2(P,"joints",14,"r1_axis",2);		
		parameters.Relbow3 		= extract2(P,"joints",14,"r1_axis",3);		
		parameters.Rraduln1 	= extract2(P,"joints",15,"r1_axis",1);		
		parameters.Rraduln2 	= extract2(P,"joints",15,"r1_axis",2);		
		parameters.Rraduln3 	= extract2(P,"joints",15,"r1_axis",3);	
        parameters.Lankle1 		= extract2(P,"joints",9,"r1_axis",1);		
		parameters.Lankle2 		= extract2(P,"joints",9,"r1_axis",2);		
		parameters.Lankle3 		= extract2(P,"joints",9,"r1_axis",3);		
		parameters.Lsubtalar1 	= extract2(P,"joints",10,"r1_axis",1);		
		parameters.Lsubtalar2 	= extract2(P,"joints",10,"r1_axis",2);		
		parameters.Lsubtalar3 	= extract2(P,"joints",10,"r1_axis",3);		
		parameters.Lmtp1 		= extract2(P,"joints",11,"r1_axis",1);		
		parameters.Lmtp3 		= extract2(P,"joints",11,"r1_axis",3);		
		parameters.Lelbow1 		= extract2(P,"joints",18,"r1_axis",1);		
		parameters.Lelbow2 		= extract2(P,"joints",18,"r1_axis",2);		
		parameters.Lelbow3 		= extract2(P,"joints",18,"r1_axis",3);		
		parameters.Lraduln1 	= extract2(P,"joints",19,"r1_axis",1);		
		parameters.Lraduln2 	= extract2(P,"joints",19,"r1_axis",2);		
		parameters.Lraduln3 	= extract2(P,"joints",19,"r1_axis",3);	
		normalize(&parameters.Rankle1,    &parameters.Rankle2,    &parameters.Rankle3);
		normalize(&parameters.Rsubtalar1, &parameters.Rsubtalar2, &parameters.Rsubtalar3);
		normalize(&parameters.Rmtp1,      &zero,                  &parameters.Rmtp3);
		normalize(&parameters.Relbow1,    &parameters.Relbow2,    &parameters.Relbow3);
		normalize(&parameters.Rraduln1,   &parameters.Rraduln2,   &parameters.Rraduln3);
        normalize(&parameters.Lankle1,    &parameters.Lankle2,    &parameters.Lankle3);
		normalize(&parameters.Lsubtalar1, &parameters.Lsubtalar2, &parameters.Lsubtalar3);
		normalize(&parameters.Lmtp1,      &zero,                  &parameters.Lmtp3);
		normalize(&parameters.Lelbow1,    &parameters.Lelbow2,    &parameters.Lelbow3);
		normalize(&parameters.Lraduln1,   &parameters.Lraduln2,   &parameters.Lraduln3);
		}
		
		// extract and store all mass properties
		{
			// segments 1, 5, 10 in the opensim model have no mass properties (ground, Rtalus, Ltalus)
			EXTRACTINERTIAL(pelvis, 2);
			EXTRACTINERTIAL(Rfemur, 3);
			EXTRACTINERTIAL(Rtibia, 4);
			EXTRACTINERTIAL(Rcalcaneus, 6);
			EXTRACTINERTIAL(Rtoes, 7);
			EXTRACTINERTIAL(Lfemur, 8);
			EXTRACTINERTIAL(Ltibia, 9);
			EXTRACTINERTIAL(Lcalcaneus, 11);
			EXTRACTINERTIAL(Ltoes, 12);
			EXTRACTINERTIAL(torso, 13);
			EXTRACTINERTIAL(Rhumerus, 14);
			EXTRACTINERTIAL(Rulna, 15);
			EXTRACTINERTIAL(Rradius, 16);
			EXTRACTINERTIAL(Rhand, 17);	
			EXTRACTINERTIAL(Lhumerus, 18);
			EXTRACTINERTIAL(Lulna, 19);
			EXTRACTINERTIAL(Lradius, 20);
			EXTRACTINERTIAL(Lhand, 21);	
		}
		
		// determine the body mass and weight
		{
			double g = sqrt(	parameters.gravity_x*parameters.gravity_x + 
								parameters.gravity_y*parameters.gravity_y + 
								parameters.gravity_z*parameters.gravity_z); 
			double bodymass = 	parameters.pelvis_M +
								parameters.torso_M +
								parameters.Rfemur_M +
								parameters.Rtibia_M +
								parameters.Rcalcaneus_M +
								parameters.Rtoes_M +
								parameters.Rhumerus_M +
								parameters.Rulna_M +
								parameters.Rradius_M +
								parameters.Rhand_M +
								parameters.Lfemur_M +
								parameters.Ltibia_M +
								parameters.Lcalcaneus_M +
								parameters.Ltoes_M +
								parameters.Lhumerus_M +
								parameters.Lulna_M +
								parameters.Lradius_M +
								parameters.Lhand_M;
			parameters.bodyweight = g * bodymass;// todo add body mass of missing segments (talus) and make mass 0 in osim file
		}
		
		// from P, extract the muscles and their properties
		// Hamner's model used the Thelen 2003 muscle model so we translate those to our muscle model
		{
			mxArray *field;
			int i,j,k;
			nmuscles = (int) extract(P, "nMus", 1);
			if (nmuscles > MAXMUSCLES) {
				mexErrMsgTxt("gait3d_pelvis213:Initialize: too many muscles.");
			}
			if ((field = mxGetField(P,0,"muscles") ) == NULL) {
				mexErrMsgTxt("gait3d_pelvis213:Initialize: no muscles field found in model.");
			}
			if (!mxIsCell(field)) {
					mexErrMsgTxt("gait3d_pelvis213: Initialize: muscles field is not a cell array.");
			}
			for (i=0; i<nmuscles; i++) {
				mxArray *cell, *fieldincell;
				// printf("Extracting fields from muscle %d\n", i);
				cell = mxGetCell(field, i);
				
				// extract muscle name
				if ((fieldincell = mxGetField(cell,0,"name") ) == NULL) {
					mexErrMsgTxt("gait3d_pelvis213:Initialize: no name field found in muscle.");
				}
				mxGetString(fieldincell, muscles[i].name, NAMELENGTH);
				
				// extract muscle properties
				muscles[i].Lceopt 		= extract(cell, "lceopt" , 1);
				muscles[i].pennopt 		= extract(cell, "pennatopt" , 1);
				muscles[i].width 		= sqrt(extract(cell, "kactive" , 1));		// see Thelen 2003, our width is square root of his gamma parameter
				muscles[i].Fmax 		= extract(cell, "fmax" , 1);
				muscles[i].Vmax 		= extract(cell, "vmax" , 1);
				muscles[i].Tact 		= extract(cell, "tact" , 1);
				muscles[i].Tdeact 		= extract(cell, "tdeact" , 1);
				muscles[i].gmax 		= extract(cell, "flen" , 1);
				muscles[i].SEEslack 	= extract(cell, "lslack" , 1);
				muscles[i].PEEslack 	= 1.0;					// resting length equal to Lceopt, same as Thelen
				muscles[i].umax		 	= 0.04;					// Thelen does the SEE differently, but this is not critical
				muscles[i].krel			= 1.0;					// Thelen does the PEE differently, but this is not critical
				muscles[i].HillA		= 0.25;					// Thelen (2003) also used this value, and we used it before

				// some other properties are derived:
				muscles[i].c3 = muscles[i].Vmax * muscles[i].HillA * (muscles[i].gmax - 1.0) / (muscles[i].HillA + 1.0);
				muscles[i].kSEE = 1.0/(muscles[i].umax*muscles[i].umax*muscles[i].SEEslack*muscles[i].SEEslack);
				muscles[i].kPEE = muscles[i].krel / (muscles[i].width * muscles[i].width);

				// initialize a matrix that can tell us quickly if the muscle spans a DOF
				for (j=0; j<NDOF; j++) muscles[i].spans[j] = 0;
				
				// polynomials
				muscles[i].nmusdof		= (int) extract(cell, "dof_count", 1);
				if (muscles[i].nmusdof > MAXMUSDOF) {
					printf("Muscle number: %d\n", i+1);
					mexErrMsgTxt("gait3d_pelvis213:Initialize: too many degrees of freedom in muscle path.");
				}
				for (j=0; j<muscles[i].nmusdof; j++) {
					muscles[i].musdof[j] = (int) extract(cell, "dof_indexes", j+1) - 1;	// our dof index starts at 0 in this C code
				}
				muscles[i].npolterms	= (int) extract(cell, "lparam_count", 1);
				if (muscles[i].npolterms > MAXPOLTERMS) {
					printf("Muscle number: %d\n", i+1);
					mexErrMsgTxt("gait3d_pelvis213:Initialize: too many polynomial terms in muscle path.");
				}
						
				for (j=0; j<muscles[i].npolterms; j++) {
					muscles[i].polcoeff[j] = extract(cell, "lcoefs", j+1);
					for (k=0; k<muscles[i].nmusdof; k++) {
						muscles[i].expon[j][k] = (int) extract(cell, "lparams", muscles[i].npolterms*k+j+1);
						
						// if there is a nonzero exponent in any of the terms, there is a nonzero moment arm
						if (muscles[i].expon[j][k] > 0) {
							muscles[i].spans[muscles[i].musdof[k]] = 1;
						}
					}
				}
			}
		}
		
		// from P, extract the contact elements and their properties
		{
			mxArray *field;
			int i;
			ncontacts = (int) extract(P, "nCPs", 1);
			if (ncontacts > MAXCONTACTS)
				mexErrMsgTxt("gait3d_pelvis213:Initialize: too many contact elements.");
			if ((field = mxGetField(P,0,"CPs") ) == NULL) {
				mexErrMsgTxt("gait3d_pelvis213:Initialize: no field named CPs found in model.");
			}
			if (!mxIsCell(field)) {
				mexErrMsgTxt("gait3d_pelvis213: Initialize: contacts field is not a cell array.");
			}
			for (i=0; i<ncontacts; i++) {
				mxArray *cell;
				// printf("Extracting fields from contact element %d\n", i);
				cell = mxGetCell(field, i);
				// subtract 2 because our FK array not include ground, AND index should start at zero
				contacts[i].segment = (int) extract(cell, "segmentindex", 1) - 2;		
				contacts[i].xp 		= extract(cell, "position", 1);
				contacts[i].yp 		= extract(cell, "position", 2);
				contacts[i].zp 		= extract(cell, "position", 3);
				
				// determine which of the GRF segments this contact element is on
				if 		(contacts[i].segment == 4) {
					contacts[i].grfseg = 0;		// calcn_r
				}
				else if (contacts[i].segment == 5) {
					contacts[i].grfseg = 1;		// toes_r
				}
				else if (contacts[i].segment == 9) {
					contacts[i].grfseg = 2;		// calcn_l
				}
				else if (contacts[i].segment == 10) {
					contacts[i].grfseg = 3;		// toes_l
				}
				else {
					mexErrMsgTxt("gait3d_pelvis213: Initialize: contact point attached to a non-foot segment.");
				}
				
			}

		}
		
		// from P, extract the passive joint properties
		{
			int i;
            
            // translation and orientation relative to the ground
            // => We need only the neutral position to return xneutral during initialization
            for (i=0; i<6   ; i++) {
                joints[i].qneutral  = extract2(P,"dofs",i+1,"neutral_position",1);
            }
            
            // joint angles
			for (i=6; i<NDOF; i++) {
				joints[i].qneutral  = extract2(P,"dofs",i+1,"neutral_position",1);
                joints[i].K1 		= extract2(P,"dofs",i+1,"stiffness_K1",1);
				joints[i].K2 		= extract2(P,"dofs",i+1,"stiffness_K2",1);
				joints[i].B 		= extract2(P,"dofs",i+1,"damping_B",1);
				joints[i].qmin_muscleMoment	 = extract2(P,"dofs",i+1,"range_muscleMoment",1);
				joints[i].qmax_muscleMoment  = extract2(P,"dofs",i+1,"range_muscleMoment",2);
                joints[i].qmin_passiveMoment = extract2(P,"dofs",i+1,"range_passiveMoment",1);
				joints[i].qmax_passiveMoment = extract2(P,"dofs",i+1,"range_passiveMoment",2);
			}
		}
		
		// extract the strain energy polynomial, if there is one
		nstrainterms = 0;
				
		// set the initialized flag and compute some constants
		initialized = 1;
		nstates = 2*NDOF + 2*nmuscles + NCVAR*ncontacts;
		if (nstates > MAXSTATES)
			mexErrMsgTxt("gait3d_pelvis213:Initialize: too many states.");
		nf = 2*NDOF + 2*nmuscles + NCF*ncontacts;					// number of elements in f(x,xdot,u,M)
		
		// fill an array with zeros
		{
			int i;
			for (i=0; i<nstates; i++) zeros[i] = 0.0;
		}	
		
		// create the init struct for output
		{
			const char *fieldnames[] = { "Nx", "Nf", "Bodyweight", "fmin", "fmax"};
			mxArray *field_value;
			plhs[0] = mxCreateStructMatrix(1, 1, 5, fieldnames);
			
			field_value = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(field_value) = nstates;
			mxSetField(plhs[0], 0, "Nx", field_value);
			
			field_value = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(field_value) = nf;
			mxSetField(plhs[0], 0, "Nf", field_value);
			
			field_value = mxCreateDoubleMatrix(1,1,mxREAL); 
			*mxGetPr(field_value) = parameters.bodyweight;
			mxSetField(plhs[0], 0, "Bodyweight", field_value);
		}
		
		// create lower and upper bounds for residuals f		
		{
			double *fmin, *fmax;
			int i,j,k;
			mxArray *field_val_1;
			mxArray *field_val_2;
			field_val_1 = mxCreateDoubleMatrix(nf, 1, mxREAL);
			fmin = mxGetPr(field_val_1);
			field_val_2 = mxCreateDoubleMatrix(nf, 1, mxREAL);
			fmax = mxGetPr(field_val_2);
			j = 0;
			
			// for the multibody dynamics and muscle dynamics, lb=ub=0 (equality constraints f=0)
			for (i=0; i<2*NDOF+2*nmuscles; i++) {
				fmin[j]   = 0.0;
				fmax[j++] = 0.0;
			}

			// for each contact element, 6 equality constraints (f=0)
			for (i=0; i<ncontacts; i++) {
				for (k=0; k<NCF; k++) {
					fmin[j]   = 0.0;
					fmax[j++] = 0.0;
				}
			}
			mxSetField(plhs[0], 0, "fmin", field_val_1);
			mxSetField(plhs[0], 0, "fmax", field_val_2);
		}
		
		return;
	}
	else {
		// commands other than Initialize are done here
		
		// MEX function pointers to inputs from Matlab
		double *x, *xdot, *u, *Mextra, *Gextra;
		
		// Input and output variables for Autolev
		double *q, *qd, *qdd;
		double zero[NDOF];
		double dz_dq[NDOF][NDOF];
		double dz_dqd[NDOF][NDOF];
		double dz_dqdd[NDOF][NDOF];
		double dz_dG[NDOF][NGRF];
		double fk[NFK];
		double dfk_dq[NFK][NDOF];
		double fkdot[NFK];
		double dfkdot_dq[NFK][NDOF];		
		
		// Muscle variables
		double Lm[MAXMUSCLES];							// muscle+tendon length, based on skeleton kinematic state
		double dLm_dq[MAXMUSCLES][NDOF];				// derivatives of Lm with respect to joint angles
		double dLm_dqq[MAXMUSCLES][NDOF][NDOF]; 		// second derivative of Lm w.r.t. joint angles
		double g[MAXMUSCLES];							// muscle force imbalance
		double dg_da[MAXMUSCLES], dg_ds[MAXMUSCLES], dg_dsdot[MAXMUSCLES], dg_dLm[MAXMUSCLES]; 	// derivatives of muscle imbalance
		double h[MAXMUSCLES];							// activation dynamics imbalance
		double dh_da[MAXMUSCLES], dh_dadot[MAXMUSCLES], dh_du[MAXMUSCLES];
		double force[MAXMUSCLES];						// muscle forces
		double dforce_dLm[MAXMUSCLES], dforce_ds[MAXMUSCLES];	// and its derivatives
        double forceCE[MAXMUSCLES];						// muscle forces of CE
		double dforceCE_da[MAXMUSCLES], dforceCE_ds[MAXMUSCLES], dforceCE_dsdot[MAXMUSCLES];	// and its derivatives
		double CEpower[MAXMUSCLES];
		
		// Joint variables
		double M[NDOF];							// total actuation to be supplied to Autolev generated code
		double dM_dq[NDOF][NDOF];
		double dM_dqd[NDOF];					// is a diagonal matrix (no coupling between joints), so no need for full storage
		double dM_ds[NDOF][MAXMUSCLES];
		double Mpas [NDOF];						// passive joint moments
		double dMpas_dq[NDOF][NDOF];
		double dMpas_dqd[NDOF];

		
		// GRF variables
		double G[NGRF];								// 3D ground reaction force/moment for the segments that have contact elements
		double dG_dxc[NGRF][MAXCONTACTS*NCVAR];		// derivatives of G with respect to contact variables in x

		// some logical variables related to the user's request
		int derivatives 	= (nlhs > 1);						// will be true if user requested derivatives;
		int cmdFkin 		= (strcmp(command, "Fkin") == 0);
		int cmdDynamics 	= (strcmp(command, "Dynamics") == 0);
		int cmdGRF 			= (strcmp(command, "GRF") == 0);
		int cmdJointmoments	= (strcmp(command, "Jointmoments") == 0);
		int cmdPassiveJointmoments	= (strcmp(command, "PassiveJointmoments") == 0);
		int cmdMuscleforces	= (strcmp(command, "Muscleforces") == 0);
        int cmdMuscleCEforces	= (strcmp(command, "MuscleCEforces") == 0);
		int cmdMuscleCEpower= (strcmp(command, "MuscleCEpower") == 0);

		// Give error message if model was not initialized
		if (!initialized) {
			printf("gait3d_pelvis213: command given: %s\n", command);
			mexErrMsgTxt("gait3d_pelvis213: model was not initialized.");
		}

		// get all the inputs that are needed, depending on what the command is
		if (cmdDynamics) {
			if (nrhs < 4 || nrhs > 6) {
				mexErrMsgTxt("gait3d_pelvis213:Dynamics: Must have 4,5 or 6 inputs.");			
			}
		
			// get x, xdot, u
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:Dynamics: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);
			
			if ((!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != nstates) || (mxGetN(prhs[2]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:Dynamics: Incorrect input xdot.");
			}
			xdot = mxGetPr(prhs[2]);
			
			if ((!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetM(prhs[3]) != nmuscles) || (mxGetN(prhs[3]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:Dynamics: Incorrect input for u.");
			}
			u = mxGetPr(prhs[3]);	
			
			// get optional input M (additional actuation)
			if (nrhs > 4) {
				if ((!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetM(prhs[4]) != NDOF) || (mxGetN(prhs[4]) != 1)) {
					mexErrMsgTxt("gait3d_pelvis213:Dynamics Incorrect input for M.");
				}
				Mextra = mxGetPr(prhs[4]);
			}
            
            // get G as input 
			if (nrhs > 5) {
				if ((!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetM(prhs[5]) != NGRF) || (mxGetN(prhs[5]) != 1)) {
					mexErrMsgTxt("gait3d_pelvis213:Dynamics Incorrect input for G.");
				}
				Gextra = mxGetPr(prhs[5]);
			}
            
			
			// generate pointers to q, qd, and qdd
			q = &x[0];
			qd = &x[NDOF];	
			qdd = &xdot[NDOF];
			
		}
		if (cmdMuscleforces) {
			if (nrhs != 2) {
				mexErrMsgTxt("gait3d_pelvis213:Muscleforce: Must have 2 inputs.");			
			}
		
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:Dynamics: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);
			xdot = zeros;
		}
        if (cmdMuscleCEforces) {
			if (nrhs != 3) {
				mexErrMsgTxt("gait3d_pelvis213:MuscleCEforce: Must have 3 inputs.");			
			}
		
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:MuscleCEforce: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);
			xdot = mxGetPr(prhs[2]);
		}
        if (cmdMuscleCEpower) {
			if (nrhs != 3) {
				mexErrMsgTxt("gait3d_pelvis213:MuscleCEpower: Must have 3 inputs.");			
			}
		
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:MuscleCEpower: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);
			xdot = mxGetPr(prhs[2]);
		}

		if (cmdGRF) {
			if (nrhs != 2) {
				mexErrMsgTxt("gait3d_pelvis213:GRF: Must have exactly 2 inputs.");			
			}
		
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:GRF: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);	
		}

		if (cmdJointmoments) {
			if (nrhs < 2 || nrhs > 3) {
				mexErrMsgTxt("gait3d_pelvis213:Jointmoments: Must have 2 or 3 inputs.");			
			}
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:Jointmoments: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);	
			xdot = zeros;				// pointer to zeros, to avoid crash when xdot is not initialized
			
			// generate pointers to q, qd
			q = &x[0];
			qd = &x[NDOF];	
            
            if (nrhs > 2) {
                if ((!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != NDOF) || (mxGetN(prhs[2]) != 1)) {
                    mexErrMsgTxt("gait3d_pelvis213:Jointmoments: Incorrect input for M.");
                }
                Mextra = mxGetPr(prhs[2]);
            }
            
		}
		
		if (cmdPassiveJointmoments) {
			if (nrhs < 2 || nrhs > 2) {
				mexErrMsgTxt("gait3d_pelvis213:PassiveJointmoments: Must have 2 inputs.");			
			}
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:PassiveJointmoments: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);	
			xdot = zeros;				// pointer to zeros, to avoid crash when xdot is not initialized
			
			// generate pointers to q, qd
			q = &x[0];
			qd = &x[NDOF];	
		}

		if (cmdFkin) {
			int i;
			
			if (nrhs < 2) {
				mexErrMsgTxt("gait3d_pelvis213:Fkin: must have at least one input.");			
			}
			
			// get q
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != NDOF) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213:Fkin Incorrect input for q.");
			}
			q = mxGetPr(prhs[1]);
			
			// case 1: 1 or 2 outputs requested
			if (nlhs < 3) {
				if (nrhs != 2) {
					mexErrMsgTxt("gait3d_pelvis213:Fkin: Must have exactly 2 inputs when 2 outputs are requested.");			
				}
				qd = zeros;		// velocities are not needed
			}
			else if (nlhs < 4) {
				if (nrhs != 3) {
					mexErrMsgTxt("gait3d_pelvis213:Fkin: Must have exactly 3 inputs when 3 outputs are requested.");			
				}
				// get qd (dq/dt)
				if ((!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != NDOF) || (mxGetN(prhs[2]) != 1)) {
					mexErrMsgTxt("gait3d_pelvis213:Fkin Incorrect input for qd.");
				}
				qd = mxGetPr(prhs[2]);
			}
			else {
				mexErrMsgTxt("gait3d_pelvis213:Fkin: Too many outputs requested.");			
			}

			// set some other things to zero before we call Autolev generated code
			qdd = zeros;
			for (i=0; i<NDOF; i++) M[i] = 0.0;
			for (i=0; i<NGRF; i++) G[i] = 0.0;
		}

		// Compute the muscle dynamics, and get muscle forces
		if (cmdDynamics || cmdMuscleforces || cmdMuscleCEforces || cmdMuscleCEpower || cmdJointmoments) {
			int i;
			// Calculate muscle-tendon lengths Lm and first and second derivatives with respect to q
			MusclePaths(x, Lm, dLm_dq, dLm_dqq, derivatives);
			for(i=0; i<nmuscles; i++) {	
				int is = 2*NDOF + i;													// index to s state variable
				int iact = 2*NDOF + nmuscles + i;										// index to activation state variable

				// Calculate imbalance of activation dynamics equation: adot = rate * (u-a)
				if (cmdDynamics) {
					double rate = u[i]/muscles[i].Tact + (1-u[i])/muscles[i].Tdeact;		// activation rate depends on u (He et al 1991)
					h[i] = xdot[iact] - rate*(u[i] - x[iact]);								// imbalance
					if (derivatives) {
						dh_da[i] = rate;
						dh_dadot[i] = 1.0;
						dh_du[i] = -rate - (1/muscles[i].Tact - 1/muscles[i].Tdeact)*(u[i] - x[iact]);
					}
				}
				
				// Calculate muscle force imbalance, normalized to Fmax
				g[i] = MuscleDynamics(&muscles[i],
					x[iact],				// active state of muscle i
					x[is],					// s (projected Lce) of muscle i
					xdot[is],				// sdot
					Lm[i],					// muscle length
					&dg_da[i], &dg_ds[i], &dg_dsdot[i], &dg_dLm[i],
					&force[i], &dforce_dLm[i], &dforce_ds[i],
                    &forceCE[i], &dforceCE_da[i], &dforceCE_ds[i], &dforceCE_dsdot[i],
					&CEpower[i]);
	
			}
		}

		// Evaluate the strain energy polynomial, with gradient and hessian
	#if 0
		if ((nstrainterms > 0) && (cmdDynamics || cmdJointmoments)) {
			int i,j,k;
			double strainenergy;
			double polgrad[NDOF];
			double polhess[NDOF][NDOF];
			strainenergy = 0.0;
			for (i=0; i<NDOF; i++) {
				polgrad[i] = 0.0;
				for (j=0; j<NDOF; j++) {
					polhess[i][j] = 0.0;
				}
			}
			for (i=0; i < nstrainterms; i++) {
				double term, dterm;
			
				// calculate this term's contribution to the strain energy
				term = straincoeff[i];
				if (term != 0.0) {
					for (j=0; j < NDOF; j++) {
						for (k=0; k < strainexpon[i][j]; k++) {
							term = term * x[3+j];			// this creates coeff[i] * product of exponentiated angles
						}
					}
					strainenergy = strainenergy + term;

					// contribution of this term to the gradient and hessian of E with respect to joint angles
					for (k=0; k < NDOF; k++) {
						int kk;
						// derivative of this term with respect to angle k is nonzero when exponent of angle k is not zero
						// this goes wrong when angle is exactly zero, but that never really happens
						if (strainexpon[i][k] > 0) {
							if (x[3+k] == 0.0) {
								mexErrMsgTxt("Elastic strain effects cannot be calculated when a joint angle is exactly zero.");		
							}
							dterm = strainexpon[i][k] * term / x[3+k];
							polgrad[k] = polgrad[k] + dterm;
						
							// off-diagonal hessian elements
							for (kk=0; kk < k; kk++) {					// no need for kk > k because the hessian is symmetrical
								if (strainexpon[i][kk] > 0) {
									if (x[3+kk] == 0.0) {
										mexErrMsgTxt("Elastic strain effects cannot be calculated when a joint angle is exactly zero.");		
									}
									polhess[k][kk] = polhess[k][kk] + strainexpon[i][kk] * dterm / x[3+kk];
									polhess[kk][k] = polhess[k][kk];		// symmetry
								}
							}
							
							// diagonal hessian elements
							if (strainexpon[i][k] > 1) {
								polhess[k][k] = polhess[k][k] + (strainexpon[i][k]-1) * dterm / x[3+k];
							}
						}
					}
				}
			}
		}
	#endif	

        // Compute the passive joint moments in Nm
        if (cmdPassiveJointmoments)  {
            int i;
			for (i=0; i<NDOF; i++) {
				int j;
			
				// initialize joint moments and derivatives to zero
				Mpas[i] = 0.0;
				if (derivatives) {
					for (j=0; j<NDOF; j++) dMpas_dq[i][j] = 0.0;
					dMpas_dqd[i] = 0.0;
				}	
				// internal moments due to passive limits and muscles,
				// these only exist for DOF 7 and higher
				if (i >= 6) {
					double d,q,qdot;
					// passive joint moment
					q = x[i];									// joint angle i is state variable i
					qdot = x[NDOF+i];							// and the corresponding angular velocity
					Mpas[i] = Mpas[i] - joints[i].K1*(q-joints[i].qneutral) - joints[i].B*qdot;		// a small linear elastic moment at all angles, plus damping
					if (derivatives) {
						dMpas_dq[i][i] = -joints[i].K1;				// and its derivatives
						dMpas_dqd[i]   = -joints[i].B;
					}

					d = q - joints[i].qmin_passiveMoment;	    // are we below min angle?
					if (d < 0) {								// yes, add quadratic term
						Mpas[i] = Mpas[i] + joints[i].K2*d*d;	
						if (derivatives) {						// and its derivative
							dMpas_dq[i][i] = dMpas_dq[i][i] + 2*joints[i].K2*d;
						}
					}
					d = q - joints[i].qmax_passiveMoment;
					if (d > 0) {								// are we above max angle?
						Mpas[i] = Mpas[i] - joints[i].K2*d*d;			// yes, add quadratic term
						if (derivatives) {
							dMpas_dq[i][i] = dMpas_dq[i][i] - 2*joints[i].K2*d;
						}
					}
                }
			}
        }
		// Assemble the MEX function outputs for the "PassiveJointmoments" command
		if (cmdPassiveJointmoments) {
			double *f;
			int i,j;
			
			double *dMpas_dx;
			mwIndex NdMpas_dx, *dMpas_dx_irs, *dMpas_dx_jcs;
			mwSize Nalloc_dMpas_dx;
			
			plhs[0] = mxCreateDoubleMatrix(NDOF, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			// fill the column vector
			for (i=0; i<NDOF; i++) f[i] = Mpas[i];
			
			// if requested, return derivatives dM/dx, copy these from the dense matrices that we already have
			if (derivatives) {	
				#if JACOBIANSIZING
					Nalloc_dMpas_dx = 10000;
				#else
					Nalloc_dMpas_dx = 352;		
				#endif
				plhs[1] = mxCreateSparse(nstates, NDOF, Nalloc_dMpas_dx, mxREAL);  // this is Jacobian transposed!
				dMpas_dx = mxGetPr(plhs[1]);
				dMpas_dx_irs = mxGetIr(plhs[1]);
				dMpas_dx_jcs = mxGetJc(plhs[1]);
				NdMpas_dx = 0;
				for (i=0; i<NDOF; i++) {
					dMpas_dx_jcs[i] = NdMpas_dx;					// store element number where this column starts			

					// rows from dM/dq
					for (j=0; j<NDOF; j++) {				
						if (dMpas_dq[i][j] != 0.0) {
							dMpas_dx_irs[NdMpas_dx] = j;				// store row number of this matrix element
							dMpas_dx[NdMpas_dx++] = dMpas_dq[i][j];		// store the value of this matrix element
						}
					}
						
					// one row, on diagonal, from dM/dqdot
					if (dMpas_dqd[i] != 0.0) {
						dMpas_dx_irs[NdMpas_dx] = NDOF+i;					// store row number of this matrix element
						dMpas_dx[NdMpas_dx++] = dMpas_dqd[i];				// store the value of this matrix element
					}
				}
				dMpas_dx_jcs[NDOF] = NdMpas_dx;						// store final element number

				if (Nalloc_dMpas_dx < NdMpas_dx) {
					mexErrMsgTxt("gait3d_pelvis213c: Insufficient memory allocation for dMpas_dx.");	
				}
				#if JACOBIANSIZING
					printf("dMpas_dx: allocated %d, actual %d\n", Nalloc_dMpas_dx, NdMpas_dx);
				#endif

			}
			return;			

		}

		// Compute the joint moments in Nm
		if (cmdDynamics || cmdJointmoments) {
			int i;
			for (i=0; i<NDOF; i++) {
				int j;
			
				// initialize joint moments and derivatives to zero
				M[i] = 0.0;
				if (derivatives) {
					for (j=0; j<NDOF; j++) dM_dq[i][j] = 0.0;
					for (j=0; j<nmuscles; j++) dM_ds[i][j] = 0.0;
					dM_dqd[i] = 0.0;
				}	
				// internal moments due to passive limits and muscles,
				// these only exist for DOF 7 and higher
				if (i >= 6) {
					double d,q,qdot;
					// passive joint moment
					q = x[i];									// joint angle i is state variable i
					qdot = x[NDOF+i];							// and the corresponding angular velocity
					M[i] = M[i] - joints[i].K1*(q-joints[i].qneutral) - joints[i].B*qdot;		// a small linear elastic moment at all angles, plus damping
					if (derivatives) {
						dM_dq[i][i] = -joints[i].K1;				// and its derivatives
						dM_dqd[i]   = -joints[i].B;
					}

					d = q - joints[i].qmin_passiveMoment;	    // are we below min angle?
					if (d < 0) {								// yes, add quadratic term
						M[i] = M[i] + joints[i].K2*d*d;	
						if (derivatives) {						// and its derivative
							dM_dq[i][i] = dM_dq[i][i] + 2*joints[i].K2*d;
						}
					}
					d = q - joints[i].qmax_passiveMoment;
					if (d > 0) {								// are we above max angle?
						M[i] = M[i] - joints[i].K2*d*d;			// yes, add quadratic term
						if (derivatives) {
							dM_dq[i][i] = dM_dq[i][i] - 2*joints[i].K2*d;
						}
					}
					
					// add the muscle moments
					for (j=0; j<nmuscles; j++) {
						int k;
						if ( muscles[j].spans[i] ) {
							M[i] = M[i] - dLm_dq[j][i]*force[j];		// moment arm is -dLm/dq
							if (derivatives) {
								for (k=0; k<NDOF; k++) {
									if ( muscles[j].spans[k] ) {
										// printf("dLm_dqq[%d][%d][%d] = %f\n",j,i,k,dLm_dqq[j][i][k]);
										dM_dq[i][k] = dM_dq[i][k] - dLm_dq[j][i]*dforce_dLm[j]*dLm_dq[j][k] - dLm_dqq[j][i][k]*force[j];
									}
								}
								dM_ds[i][j] = dM_ds[i][j] - dLm_dq[j][i]*dforce_ds[j];
							}
						}
					}
				}
			
				// add the extra actuation
				if ((cmdDynamics && nrhs > 4) || (cmdJointmoments && nrhs > 2)){
					M[i] = M[i] + Mextra[i];
                }              
                
			}
		}
		// Assemble the MEX function outputs for the "Jointmoments" command
		if (cmdJointmoments) {
			double *f;
			int i,j;
			
			double *dM_dx; 
			mwIndex NdM_dx, *dM_dx_irs, *dM_dx_jcs;
			mwSize Nalloc_dM_dx;
			
			plhs[0] = mxCreateDoubleMatrix(NDOF, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			// fill the column vector
			for (i=0; i<NDOF; i++) f[i] = M[i];
			
			// if requested, return derivatives dM/dx, copy these from the dense matrices that we already have
			if (derivatives) {	
				#if JACOBIANSIZING
					Nalloc_dM_dx = 10000;
				#else
					Nalloc_dM_dx = 352;		
				#endif
				plhs[1] = mxCreateSparse(nstates, NDOF, Nalloc_dM_dx, mxREAL);  // this is Jacobian transposed!
				dM_dx = mxGetPr(plhs[1]);
				dM_dx_irs = mxGetIr(plhs[1]);
				dM_dx_jcs = mxGetJc(plhs[1]);
				NdM_dx = 0;
				for (i=0; i<NDOF; i++) {
					dM_dx_jcs[i] = NdM_dx;					// store element number where this column starts			

					// rows from dM/dq
					for (j=0; j<NDOF; j++) {				
						if (dM_dq[i][j] != 0.0) {
							dM_dx_irs[NdM_dx] = j;				// store row number of this matrix element
							dM_dx[NdM_dx++] = dM_dq[i][j];		// store the value of this matrix element
						}
					}
						
					// one row, on diagonal, from dM/dqdot
					if (dM_dqd[i] != 0.0) {
						dM_dx_irs[NdM_dx] = NDOF+i;					// store row number of this matrix element
						dM_dx[NdM_dx++] = dM_dqd[i];				// store the value of this matrix element
					}
					
					// rows from dM/ds
					for (j=0; j<nmuscles; j++) {				
						if ( muscles[j].spans[i] ) {
							dM_dx_irs[NdM_dx] = 2*NDOF+j;	// store row number of this matrix element
							dM_dx[NdM_dx++] = dM_ds[i][j];	// store the value of this matrix element
						}
					}			
					
				}
				dM_dx_jcs[NDOF] = NdM_dx;						// store final element number

				if (Nalloc_dM_dx < NdM_dx) {
					mexErrMsgTxt("gait3d_pelvis213c: Insufficient memory allocation for dM_dx.");	
				}
				#if JACOBIANSIZING
					printf("dM_dx: allocated %d, actual %d\n", Nalloc_dM_dx, NdM_dx);
				#endif

			}
			return;			

		}
	
		// Assemble the MEX function outputs for the "Muscleforces" command
		if (cmdMuscleforces) {
			double *f;
			int i,j;
			
			// first output: muscle forces
			plhs[0] = mxCreateDoubleMatrix(nmuscles, 1, mxREAL);
			f = mxGetPr(plhs[0]);	
			for (j=0; j<nmuscles; j++) {
					*f++ = force[j];
			}
			
			// second output: muscle-tendon lengths
			plhs[1] = mxCreateDoubleMatrix(nmuscles, 1, mxREAL);
			f = mxGetPr(plhs[1]);	
			for (j=0; j<nmuscles; j++) {
					*f++ = Lm[j];
			}
			
			// third output: moment arm matrix Dij = moment arm of muscle i at DOF j
			plhs[2] = mxCreateDoubleMatrix(nmuscles, NDOF, mxREAL);
			f = mxGetPr(plhs[2]);	
			for (i=0; i<NDOF; i++) {
				for (j=0; j<nmuscles; j++) {
					*f++ = -dLm_dq[j][i];
				}
			}
			return;
		}
        
        // Assemble the MEX function outputs for the "MuscleCEforces" command
		if (cmdMuscleCEforces) {
			double *f;
			int i,j;
            
            double *dfCE_dx; 
			mwIndex NdfCE_dx, *dfCE_dx_irs, *dfCE_dx_jcs;
			mwSize Nalloc_dfCE_dx;
			double *dfCE_dxdot; 
			mwIndex NdfCE_dxdot, *dfCE_dxdot_irs, *dfCE_dxdot_jcs;
			mwSize Nalloc_dfCE_dxdot;
            
			// first output: muscle forces of CE
			plhs[0] = mxCreateDoubleMatrix(nmuscles, 1, mxREAL);
			f = mxGetPr(plhs[0]);	
			for (j=0; j<nmuscles; j++) {
					*f++ = forceCE[j];
			}
			
            // if requested, return derivatives dfCE/dx and dfCE/dxdot, copy these from the dense matrices that we already have
			if (derivatives) {	
				#if JACOBIANSIZING
					Nalloc_dfCE_dx = 10000;
                    Nalloc_dfCE_dxdot = 10000;
				#else
					Nalloc_dfCE_dx = 352;	// Used the same number as for cmdJointmoments	
                    Nalloc_dfCE_dxdot = 352;	// Used the same number as for cmdJointmoments	
				#endif
				plhs[1] = mxCreateSparse(nstates, nmuscles, Nalloc_dfCE_dx, mxREAL);  // this is Jacobian transposed!
				dfCE_dx = mxGetPr(plhs[1]);
				dfCE_dx_irs = mxGetIr(plhs[1]);
				dfCE_dx_jcs = mxGetJc(plhs[1]);
				NdfCE_dx = 0;
                plhs[2] = mxCreateSparse(nstates, nmuscles, Nalloc_dfCE_dxdot, mxREAL);  // this is Jacobian transposed!
				dfCE_dxdot = mxGetPr(plhs[2]);
				dfCE_dxdot_irs = mxGetIr(plhs[2]);
				dfCE_dxdot_jcs = mxGetJc(plhs[2]);
				NdfCE_dxdot = 0;
				for (i=0; i<nmuscles; i++) {
					dfCE_dx_jcs[i] = NdfCE_dx;					// store element number where this column starts			
                    dfCE_dxdot_jcs[i] = NdfCE_dxdot;					// store element number where this column starts			

                    // one element for each muscle, on diagonal, from dF/ds
                    dfCE_dx_irs[NdfCE_dx] = 2*NDOF+i;				 // store row number of this matrix element
                    dfCE_dx[NdfCE_dx++] = dforceCE_ds[i];			 // store the value of this matrix element

                    // one element for each muscle, on diagonal, from dF/dsdot
                    dfCE_dxdot_irs[NdfCE_dxdot] = 2*NDOF+i;			 // store row number of this matrix element
                    dfCE_dxdot[NdfCE_dxdot++] = dforceCE_dsdot[i]; // store the value of this matrix element

                    // one element for each muscle, on diagonal, from dF/da
                    dfCE_dx_irs[NdfCE_dx] = 2*NDOF+nmuscles+i;	         // store row number of this matrix element
                    dfCE_dx[NdfCE_dx++] = dforceCE_da[i];			 // store the value of this matrix element			
					
				}
				dfCE_dx_jcs[nmuscles] = NdfCE_dx;						 // store final element number
                dfCE_dxdot_jcs[nmuscles] = NdfCE_dxdot;					 // store final element number

                if (Nalloc_dfCE_dx < NdfCE_dx) {
                    mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for dfCE_dx.");			
                }      
                if (Nalloc_dfCE_dxdot < NdfCE_dxdot) {
                    mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for dFCE_dxdot.");			
                } 

				#if JACOBIANSIZING
					printf("dfce_dx: allocated %d, actual %d\n", Nalloc_dfCE_dx, NdfCE_dx);
                    printf("dfce_dxdot: allocated %d, actual %d\n", Nalloc_dfCE_dxdot, NdfCE_dxdot);
				#endif

			}
            
			return;
		}
		
		// Assemble the MEX function outputs for the "MuscleCEpower" command
		if (cmdMuscleCEpower) {
			double *f, *f2;
			int j;
			plhs[0] = mxCreateDoubleMatrix(nmuscles, 1, mxREAL);
			plhs[1] = mxCreateDoubleMatrix(nmuscles, 1, mxREAL);
			f = mxGetPr(plhs[0]);	
			f2 = mxGetPr(plhs[1]);	
			// fill the column vectors with CE power and contraction imbalance
			for (j=0; j<nmuscles; j++) {
				*f++ = CEpower[j];
				*f2++ = g[j];
			}
			return;
		}

		// combine the point contact forces from x (in BW) into total force and moment G on each segment,
		// in BW.  G will be input for the Autolev multibody dynamics function, in which it is
		// converted into N and Nm.
		if (cmdDynamics || cmdGRF) {
			int i;
			// first initialize everything to zero
			for (i=0; i<NGRF; i++) {
				int j;
				G[i] = 0.0;
                if (nrhs == 6) {
					G[i] += Gextra[i];
				}
				if (derivatives) {
					for (j=0; j<NCVAR*ncontacts; j++) {
						dG_dxc[i][j]   = 0.0;
					}
				}
			}
            
			// then add up the contact forces/moments on each segment
			for (i=0; i<ncontacts; i++) {

				int k = 6*contacts[i].grfseg;		// index within G where forces from this contact element should be added
				int jc = NCVAR*i;					// index within the contact variables in x for this contact element: : Fx,Fy,Fz,xc,yc,zc
				int j = 2*NDOF + 2*nmuscles + jc;	// index within x for variables of this contact element: : Fx,Fy,Fz,xc,yc,zc
				
				// add the XYZ forces that act on this segment (these are in global reference frame) 
				G[k]   += x[j];						
				G[k+1] += x[j+1];						
				G[k+2] += x[j+2];	

				// add the XYZ moments that act on this segment (these are in global reference frame)
				G[k+3] += x[j+4]*x[j+2] - x[j+5]*x[j+1];		// to the total Mx, add y*Fz-z*Fy
				G[k+4] += x[j+5]*x[j  ] - x[j+3]*x[j+2];		// to the total My, add z*Fx-x*Fz
				G[k+5] += x[j+3]*x[j+1] - x[j+4]*x[j  ];		// to the total Mz, add x*Fy-y*Fx
				
				if (derivatives) {
					dG_dxc[k][jc]     += 1;
					dG_dxc[k+1][jc+1] += 1;
					dG_dxc[k+2][jc+2] += 1;
					dG_dxc[k+3][jc+4] += x[j+2];
					dG_dxc[k+3][jc+2] += x[j+4];
					dG_dxc[k+3][jc+5] -= x[j+1];
					dG_dxc[k+3][jc+1] -= x[j+5];
					dG_dxc[k+4][jc+5] += x[j  ];
					dG_dxc[k+4][jc  ] += x[j+5];
					dG_dxc[k+4][jc+3] -= x[j+2];
					dG_dxc[k+4][jc+2] -= x[j+3];
					dG_dxc[k+5][jc+3] += x[j+1];
					dG_dxc[k+5][jc+1] += x[j+3];
					dG_dxc[k+5][jc+4] -= x[j  ];
					dG_dxc[k+5][jc  ] -= x[j+4];
				}
			}
             
		}
		
		// Assemble the MEX function outputs for the "GRF" command, outputs 3D force/moment on each foot
		if (cmdGRF) {	
			// To DO: this code assumes that each foot has just two segments, we should make this more general
		
			int igrf, ifoot;
			double *grf;
			double *dgrf_dx; 
			mwIndex Ndgrf_dx, *dgrf_dx_irs, *dgrf_dx_jcs;
			mwSize Nalloc_dgrf_dx;	
			
			// create the first output
			plhs[0] = mxCreateDoubleMatrix(12, 1, mxREAL);
			grf = mxGetPr(plhs[0]);

			if (derivatives) {
				// 15 derivatives for each contact point (15 nonzeros in dG_dxc, see above)
				Nalloc_dgrf_dx = 15*ncontacts;		
				plhs[1] = mxCreateSparse(nstates, 12, Nalloc_dgrf_dx, mxREAL);
				dgrf_dx = mxGetPr(plhs[1]);
				dgrf_dx_irs = mxGetIr(plhs[1]);
				dgrf_dx_jcs = mxGetJc(plhs[1]);
				Ndgrf_dx = 0;
			}
			
			igrf = 0;
			for (ifoot=0; ifoot<2; ifoot++) {				// 0: right foot, 1: left foot
				int j;
				int iG1 = 12*ifoot;						// pointer to G variables for first segment on this foot
				int iG2 = 12*ifoot + 6;					// pointer to G variables for second segment on this foot
				for (j=0; j<6; j++) {
					grf[igrf] = G[iG1+j] + G[iG2+j];	// simply add the 3D force and moment from the two segments in the foot
					if (derivatives) {
						int k;
						dgrf_dx_jcs[igrf] = Ndgrf_dx;		// store element number where this column starts			
						for (k=0; k<NCVAR*ncontacts; k++) {
							double tmp = dG_dxc[iG1+j][k] + dG_dxc[iG2+j][k];
							if (tmp != 0.0) {
								dgrf_dx_irs[Ndgrf_dx] = 2*NDOF + 2*nmuscles + k;	// store row number of this matrix element
								dgrf_dx[Ndgrf_dx++] = tmp;						// store the value of this matrix element					
							}
						}
					}
					igrf++;
				}
			}
			if (derivatives) {
				dgrf_dx_jcs[12] = Ndgrf_dx;		// store final element number
				if (Nalloc_dgrf_dx != Ndgrf_dx) {
					printf("dgrfdx: allocated %d, actual %d\n", Nalloc_dgrf_dx, Ndgrf_dx);
					if (Nalloc_dgrf_dx < Ndgrf_dx) {
						mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for dgrf_dx.");	
					}
				}
			}
					
			return;	
		}	

		// For dynamics : call the C function that was generated by Autolev
		// TODO: we can create separate functions for results with and without Jacobians
		if (cmdDynamics) {
			if (derivatives) {
				gait3d_pelvis213_al(&parameters, q, qd, qdd, G, zero, dz_dq, dz_dqd, dz_dqdd, dz_dG,
					fk, dfk_dq, fkdot, dfkdot_dq);
            } else {
				gait3d_pelvis213_NoDer_al(&parameters, q, qd, qdd, G, zero,
					fk, dfk_dq, fkdot, dfkdot_dq); // derivatives of fk could also be removed in future
			}
		}

		// For forward kinematics: call the C function that was generated by Autolev
		// TODO: we can create separate functions for results with and without Jacobians
		if (cmdFkin) {
			gait3d_pelvis213_FK_al(&parameters, q, qd, qdd, fk, dfk_dq, fkdot, dfkdot_dq);
		}
		
		// Assemble the MEX function outputs for the "Fkin" command
		if (cmdFkin) {
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
		
		// Assemble the MEX function outputs for the "Dynamics" command
		if (cmdDynamics) {
			double *f;
			int i,j,k;
			
			double *df_dx; 
			mwIndex Ndf_dx, *df_dx_irs, *df_dx_jcs;
			mwSize Nalloc_df_dx;
			
			double *df_dxdot;
			mwIndex Ndf_dxdot, *df_dxdot_irs, *df_dxdot_jcs;
			mwSize Nalloc_df_dxdot;
			
			double *df_du;
			mwIndex Ndf_du, *df_du_irs, *df_du_jcs;
			mwSize Nalloc_df_du;
			
			double *df_dM;
			mwIndex Ndf_dM, *df_dM_irs, *df_dM_jcs;
			mwSize Nalloc_df_dM;
            
            double *df_dG;
			mwIndex Ndf_dG, *df_dG_irs, *df_dG_jcs;
			mwSize Nalloc_df_dG;

			plhs[0] = mxCreateDoubleMatrix(nf, 1, mxREAL);
			f = mxGetPr(plhs[0]);

			// create MEX outputs for the sparse Jacobians
			if (derivatives) {
				// The sparse Jacobians have to be filled in column order, using Matlab sparse data structure
				// For this reason it is more convenient to return the transpose of the Jacobians: (df/dx)' etc.

				// the number of nonzeros was determined experimentally, start with 10000 and then
				// reduce based on the warning messages

				#if JACOBIANSIZING
					Nalloc_df_dx = 10000;		
					Nalloc_df_dxdot = 10000;	
				#else
					Nalloc_df_dx = 4758;		
					Nalloc_df_dxdot = 1050;	
				#endif

				// --------Jacobian (df/dx)'
				plhs[1] = mxCreateSparse(nstates, nf, Nalloc_df_dx, mxREAL);
				df_dx = mxGetPr(plhs[1]);
				df_dx_irs = mxGetIr(plhs[1]);
				df_dx_jcs = mxGetJc(plhs[1]);
				Ndf_dx = 0;

				// --------Jacobian (df/dxdot)'
				plhs[2] = mxCreateSparse(nstates, nf, Nalloc_df_dxdot, mxREAL);
				df_dxdot = mxGetPr(plhs[2]);
				df_dxdot_irs = mxGetIr(plhs[2]);
				df_dxdot_jcs = mxGetJc(plhs[2]);
				Ndf_dxdot = 0;

				// --------Jacobian (df/du)'
				Nalloc_df_du = nmuscles;				// we know there are only this many elements
				plhs[3] = mxCreateSparse(nmuscles, nf, Nalloc_df_du, mxREAL);
				df_du = mxGetPr(plhs[3]);
				df_du_irs = mxGetIr(plhs[3]);
				df_du_jcs = mxGetJc(plhs[3]);
				Ndf_du = 0;
				
				//----------Jacobian (df/dM)'
				if (nlhs > 4) {
					Nalloc_df_dM = NDOF;		// this is the exact number of nonzeros, because M are the generalized forces
					plhs[4] = mxCreateSparse(NDOF, nf, Nalloc_df_dM, mxREAL);  
					df_dM = mxGetPr(plhs[4]);
					df_dM_irs = mxGetIr(plhs[4]);
					df_dM_jcs = mxGetJc(plhs[4]);
					Ndf_dM = 0;
				}
                
                if (nlhs > 5) {
					Nalloc_df_dG = NGRF*NDOF;		// ?
					plhs[5] = mxCreateSparse(NGRF, nf, Nalloc_df_dG, mxREAL);  
					df_dG = mxGetPr(plhs[5]);
					df_dG_irs = mxGetIr(plhs[5]);
					df_dG_jcs = mxGetJc(plhs[5]);
					Ndf_dG = 0;
				}
			}
			
			// the first NDOF elements of f are: qdot-dq/dt = 0
			for (i=0; i<NDOF; i++) {
				f[i] = x[NDOF+i] - xdot[i];
				if (derivatives) {
					df_dx_jcs[i]    = Ndf_dx;				// index for start of column i of (df/dx)'
					df_dxdot_jcs[i] = Ndf_dxdot;			// index for start of column i of (df/dxdot)'
					df_du_jcs[i]    = Ndf_du;				// index for start of column i of (df/du)'
					if (nlhs > 4) df_dM_jcs[i] = Ndf_dM;	// index for start of column i of (df/dM)'
                    if (nlhs > 5) df_dG_jcs[i] = Ndf_dG;	// index for start of column i of (df/dG)'

					// ----------- df/dx
					df_dx_irs[Ndf_dx] = NDOF+i;			// store row number of this matrix element
					df_dx[Ndf_dx++] = 1.0;				// store the value of this matrix element

					// ----------- df/dxdot
					df_dxdot_irs[Ndf_dxdot] = i;		// store row number of this matrix element
					df_dxdot[Ndf_dxdot++] = -1.0;		// store the value of this matrix element
				}
			}
				
			// the next NDOF elements of f are the equations of motion from Autolev (the ZERO expressions from Kane's equations)
			for (i=0; i<NDOF; i++) {
				j = NDOF+i;
				f[j] = zero[i] + M[i] / parameters.bodyweight;	// see NOTES 7/29/2013
				if (derivatives) {
					df_dx_jcs[j]    = Ndf_dx;				// index for start of column j of (df/dx)'
					df_dxdot_jcs[j] = Ndf_dxdot;			// index for start of column j of (df/dxdot)'
					df_du_jcs[j]    = Ndf_du;				// index for start of column j of (df/du)'
					if (nlhs > 4) df_dM_jcs[j] = Ndf_dM;	// index for start of column j of (df/dM)'
                    if (nlhs > 5) df_dG_jcs[j] = Ndf_dG;	// index for start of column j of (df/dG)'

					// ----------- (df/dx)'
					// derivatives of ZERO[i] with respect to q go in rows 1..NDOF
					for (k=0; k<NDOF; k++) {
						double tmp = dz_dq[i][k] + dM_dq[i][k] / parameters.bodyweight;			// because f = z + M
						// printf("dM_dq[%d][%d]=%f\n",i,k,dM_dq[i][k]);
						if (fabs(tmp) > 1e-10) {
							df_dx_irs[Ndf_dx] = k;					// store row number of this matrix element
							df_dx[Ndf_dx++] = tmp;					// store the value of this matrix element
						}
					}
					
					// derivatives of ZERO[i] with respect to qd go in rows NDOF+1 to 2*NDOF
					for (k=0; k<NDOF; k++) {
						double tmp = dz_dqd[i][k];					// start with what came from Autolev
						// add the contribution dM/dqdot, but that's only on the diagonal
						if (k == i) tmp = tmp + dM_dqd[k] / parameters.bodyweight;
						if (tmp != 0.0) {
							df_dx_irs[Ndf_dx] = NDOF+k;				// store row number of this matrix element
							df_dx[Ndf_dx++] = tmp;					// store the value of this matrix element
						}
					}	
					
					// derivatives of ZERO[i] with respect to s should go into rows 2*NDOF+(1..nmuscles)
					for (k=0; k<nmuscles; k++) {
						if ( muscles[k].spans[i] ) {
							df_dx_irs[Ndf_dx] = 2*NDOF+k;		// store row number of this matrix element
							df_dx[Ndf_dx++] = dM_ds[i][k] / parameters.bodyweight;				// store the value of this matrix element
						}
					}
					
					// derivatives of ZERO[i] with respect to ground contact forces should go into rows 2*NDOF+2*NMUS+1 to end
					for (k=0; k<NCVAR*ncontacts; k++) {
						int kk;
						double tmp = 0.0;
						// add the contributions dz/dG * dG/dx
						for (kk=0; kk<NGRF; kk++) {
							tmp = tmp + dz_dG[i][kk]*dG_dxc[kk][k];
						}
						if (fabs(tmp) > 1e-10) {
							df_dx_irs[Ndf_dx] = 2*NDOF + 2*nmuscles + k;		// store row number of this matrix element
							df_dx[Ndf_dx++] = tmp;								// store the value of this matrix element
						}
					}				

					// ----------- (df/dxdot)'
					// derivatives of ZERO[i] with respect to qdd should go into rows NDOF+1 to 2*NDOF of df/dxdot
					for (k=0; k<NDOF; k++) {
						if (dz_dqdd[i][k] != 0.0) {
							df_dxdot_irs[Ndf_dxdot] = NDOF+k;		// store row number of this matrix element
							df_dxdot[Ndf_dxdot++] = dz_dqdd[i][k];	// store the value of this matrix element
						}
					}
					
					// ------------ df/dM comes from dzero/dM which is, by definition, the identity matrix (see NOTES 7/29/2013)
					if (nlhs > 4) {
						df_dM_irs[Ndf_dM] = i;		// store row number of this matrix element of (df/dM)'
						df_dM[Ndf_dM++] = 1.0 / parameters.bodyweight;		// store the value of this matrix element		
					}
                    
                    //------------ df/dG derivatives of ZERO[i] with respect to G 
                    if (nlhs > 5) {
                        for (k=0; k<NGRF; k++) {
                            if  (dz_dG[i][k]!= 0.0) {
                                df_dG_irs[Ndf_dG] = k;		// store row number of this matrix element of (df/dG)'
                                df_dG[Ndf_dG++] = dz_dG[i][k];
                            }
						}
					}
				}
			}
			
			// the next nmuscles elements of f are the muscle contraction dynamics
			for (i=0; i<nmuscles; i++) {
				j = 2*NDOF+i;
				f[j] = g[i];

				if (derivatives) {
					int k;
					
					df_dx_jcs[j]    = Ndf_dx;				// index for start of column j of (df/dx)'
					df_dxdot_jcs[j] = Ndf_dxdot;			// index for start of column j of (df/dxdot)'
					df_du_jcs[j]    = Ndf_du;				// index for start of column j of (df/du)'
					if (nlhs > 4) df_dM_jcs[j] = Ndf_dM;	// index for start of column j of (df/dM)'
                    if (nlhs > 5) df_dG_jcs[j] = Ndf_dG;	// index for start of column j of (df/dG)'

					// ----------- df/dx
					// derivatives of muscle imbalance with respect to q are in rows 1..NDOF
					for (k=0; k<NDOF; k++) {
						// element only exists if muscle i spans DOF k
						if ( muscles[i].spans[k] ) {
							df_dx_irs[Ndf_dx] = k;							// store row number of this matrix element
							df_dx[Ndf_dx++] = dg_dLm[i]*dLm_dq[i][k];		// store the value of this matrix element
						}
					}

					// derivatives of muscle imbalance with respect to s are diagonal, rows 2*NDOF+1 to 2*NDOF+nmuscles
					df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
					df_dx[Ndf_dx++] = dg_ds[i];		// store the value of this matrix element

					// derivatives of muscle imbalance with respect to Act are diagonal, rows 2*NDOF+nmuscles+(1..nmuscles)
					df_dx_irs[Ndf_dx] = 2*NDOF+nmuscles+i;	// store row number of this matrix element
					df_dx[Ndf_dx++] = dg_da[i];				// store the value of this matrix element
					
					// derivatives of muscle imbalance with respect to sdot are diagonal, rows 2*NDOF to 2*NDOF+nmuscles
					df_dxdot_irs[Ndf_dxdot] = 2*NDOF+i;			// store row number of this matrix element
					df_dxdot[Ndf_dxdot++] = dg_dsdot[i];		// store the value of this matrix element
				}			

			}
			
			// the next nmuscles elements of f are the muscle activation dynamics: h(a, da/dt, u) = 0
			for (i=0; i<nmuscles; i++) {
				j = 2*NDOF+nmuscles+i;
				f[j] = h[i];

				if (derivatives) {
					df_dx_jcs[j]    = Ndf_dx;				// index for start of column j of (df/dx)'
					df_dxdot_jcs[j] = Ndf_dxdot;			// index for start of column j of (df/dxdot)'
					df_du_jcs[j]    = Ndf_du;				// index for start of column j of (df/du)'
					if (nlhs > 4) df_dM_jcs[j] = Ndf_dM;	// index for start of column j of (df/dM)'
                    if (nlhs > 5) df_dG_jcs[j] = Ndf_dG;	// index for start of column j of (df/dG)'

					// ----------- df/dx
					// derivatives of activation dynamics with respect to Act are diagonal, rows 2*NDOF+nmuscles+(1..nmuscles)
					df_dx_irs[Ndf_dx] = 2*NDOF+nmuscles+i;		// store row number of this matrix element
					df_dx[Ndf_dx++] = dh_da[i];					// store the value of this matrix element

					// ----------- df/dxdot
					// derivatives of activation dynamics with respect to Actdot are diagonal, rows 2*NDOF+nmuscles+(1..nmuscles)
					df_dxdot_irs[Ndf_dxdot] = 2*NDOF+nmuscles+i;		// store row number of this matrix element
					df_dxdot[Ndf_dxdot++] = dh_dadot[i];				// store the value of this matrix element

					// ----------- df/du
					// derivatives of activation dynamics with respect to u are diagonal, rows 1..nmuscles
					df_du_irs[Ndf_du] = i;								// store row number of this matrix element
					df_du[Ndf_du] =  dh_du[i];							// store the value of this matrix element
					Ndf_du++;				

				}

			}

			// then there are NCF elements of f for each contact point
			for (i=0; i<ncontacts; i++) {
				double df_dfk[NCF][12];
				double df_dfkdot[NCF][12];
				double df_dxc[NCF][NCVAR];
				double df_dxcdot[NCF][NCVAR];
				
				// index to the NCVAR contact variables of this contact point
				int ixc = 2*NDOF + 2*nmuscles + NCVAR*i;			

				// index to the NCF contact residuals of this contact point
				int ifc = 2*NDOF + 2*nmuscles + NCF*i;			

				// get pointer to forward kinematic results for this segment
				int ifkptr = 12*contacts[i].segment;			// there are 12 values in the fk array for each segment

				// use the C code generated by contact.al
				contact_al(&contacts[i], &fk[ifkptr], &fkdot[ifkptr], &x[ixc], &xdot[ixc],
					&f[ifc], df_dfk, df_dfkdot, df_dxc, df_dxcdot);
					
				if (derivatives) {
					for (j=0; j<NCF; j++) {
						df_dx_jcs[ifc+j]    = Ndf_dx;				// index for start of column ifc+j of (df/dx)'
						df_dxdot_jcs[ifc+j] = Ndf_dxdot;			// index for start of column ifc+j of (df/dxdot)'
						df_du_jcs[ifc+j]    = Ndf_du;				// index for start of column ifc+j of (df/du)'
						if (nlhs > 4) df_dM_jcs[ifc+j] = Ndf_dM;	// index for start of column ifc+j of (df/dM)'
                        if (nlhs > 5) df_dG_jcs[ifc+j] = Ndf_dG;	// index for start of column ifc+j of (df/dG)'
						
						// compute derivatives of f[ifc+j] with respect to q
						for (k=0; k<NDOF; k++) {
							int kk;
							double tmp = 0.0;
							for (kk=0; kk<12; kk++) {
								tmp = tmp + df_dfk[j][kk]*dfk_dq[ifkptr+kk][k] + df_dfkdot[j][kk]*dfkdot_dq[ifkptr+kk][k];
							}
							if (tmp != 0.0) {
								df_dx_irs[Ndf_dx] = k;				// row number of matrix element
								df_dx[Ndf_dx++] = tmp;   			// store the value of this matrix element
							}
						}

						// compute derivatives of f[ifc+j] with respect to qdot
						for (k=0; k<NDOF; k++) {
							int kk;
							double tmp = 0.0;
							for (kk=0; kk<12; kk++) {
								tmp = tmp + df_dfkdot[j][kk]*dfk_dq[ifkptr+kk][k];		// dfkd_dqdot = dfk_dq!
							}
							if (tmp != 0.0) {
								df_dx_irs[Ndf_dx] = NDOF+k;			// row number of matrix element
								df_dx[Ndf_dx++] = tmp;   			// store the value of this matrix element
							}
						}
							
						// compute derivatives of f[ifc+j] with respect to contact variables in x
						for (k=0; k<NCVAR; k++) {
							double tmp = df_dxc[j][k];
							if (tmp != 0.0) {
								df_dx_irs[Ndf_dx] = ixc+k;			// row number of matrix element
								df_dx[Ndf_dx++] = tmp;				// store the value of this matrix element
							}
						}
						// compute derivatives of f[ifc+j] with respect to contact point velocities in xdot
						for (k=0; k<NCVAR; k++) {
							double tmp = df_dxcdot[j][k];
							if (tmp != 0.0) {
								df_dxdot_irs[Ndf_dxdot] = ixc+k;	// row number of matrix element
								df_dxdot[Ndf_dxdot++] = tmp;		// store the value of this matrix element
							}
						}
					}
				}
			}
			
			// store final column pointers for the sparse Jacobians, and check whether memory allocation was correct
			if (derivatives) {
				df_dx_jcs[nf]    			= Ndf_dx;			
				df_dxdot_jcs[nf] 			= Ndf_dxdot;		
				df_du_jcs[nf]       		= Ndf_du;				
				if (nlhs > 4) df_dM_jcs[nf] = Ndf_dM;
                if (nlhs > 5) df_dG_jcs[nf] = Ndf_dG;
		
				// Try to catch insufficient memory allocation.
				// Unfortunately, that usuallly causes a crash before we even get here.
				if (Nalloc_df_dx < Ndf_dx) {
					mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for df_dx.");	
				}
				if (Nalloc_df_dxdot < Ndf_dxdot) {
					mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for df_dxdot.");	
				}
				if (Nalloc_df_du < Ndf_du) {
					mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for df_du.");		
				}
				if (nlhs > 4) if (Nalloc_df_dM < Ndf_dM) {
						mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for df_dM.");	
				}
                if (nlhs > 5) if (Nalloc_df_dG < Ndf_dG) {
						mexErrMsgTxt("gait3d_pelvis213: Incorrect memory allocation for df_dG.");	
				}

				#if JACOBIANSIZING
					printf("dfdx: allocated %d, actual %d\n", Nalloc_df_dx, Ndf_dx);
					printf("dfdxdot: allocated %d, actual %d\n", Nalloc_df_dxdot, Ndf_dxdot);
					printf("dfdu: allocated %d, actual %d\n", Nalloc_df_du, Ndf_du);
					if (nlhs > 4) {
						printf("dfdM: allocated %d, actual %d\n", Nalloc_df_dM, Ndf_dM);
					}
                    if (nlhs > 5) {
						printf("dfdG: allocated %d, actual %d\n", Nalloc_df_dG, Ndf_dG);
					}
				#endif
				
			}
			
			return;
			
		}
			
		// if we reach this point, the command was not recognized
		mexErrMsgTxt("gait3d_pelvis213: Command not recognized.");
	}	
		

}
