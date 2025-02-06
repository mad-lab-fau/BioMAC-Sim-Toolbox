/*=================================================================
 *
 * gait3d_pelvis213_torque.c
 *
 * Implicit differential equation for 3D torque-driven model : f(x,dx/dt,u) = 0

 * This is the source code for the MEX function gait3d_pelvis213_torque.mexw32
 * The realted musculoskeletal model is documented in the file gait3d_reference.docx
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
#include "gait3d_pelvis213_torque.h"

// set this to true when we're determining size of the Jacobians
#define JACOBIANSIZING 0

// maximum size of the model (some other size constants are in gait3d_pelvis213_torque.h)
#define MAXCONTACTS 200				// maximum number of contact points in the model
#define MAXSTATES (2*NDOF+NCVAR*MAXCONTACTS)		// max number of system state variables, including contact variables

// macro to extract all mass properties of a segment and store them as scalars in the parameters struct
// S: segment name, see gait3d_pelvis213_torque.h.  NUM: index for the segments array that came out of Opensim
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
static contactprop contacts[MAXCONTACTS];	// contains contact element parameters
static jointprop joints[NDOF];				// joint moment parameters
static int ncontacts, nstates, nf;		// model size parameters.  nf is number of elements in the model expressions f(x,xdot,u,M)
static double zeros[MAXSTATES];				// an array of zeros

//===================================================================================================
// normalize: normalize a vector to length 1, warn if length was not close enough to 1 to begin with
//===================================================================================================
void normalize(double *x, double *y, double *z) {
	double length;
	length = sqrt((*x)*(*x) + (*y)*(*y) + (*z)*(*z));
	if (fabs(length - 1.0) > 1e-5) {
		printf("gait3d_pelvis213_torque: warning: vector was not normalized: %f %f %f -> length %f\n", *x, *y, *z, length);		
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
			mexErrMsgTxt("gait3d_pelvis213_torque: Valid field with this name not found in parameters structure.");
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
			mexErrMsgTxt("gait3d_pelvis213_torque: Initialize: Valid field with this name not found in parameters structure.");
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
		mexErrMsgTxt("gait3d_pelvis213_torque: At least one input argument is needed.");
    if (!mxIsChar(prhs[0])) 
		mexErrMsgTxt("gait3d_pelvis213_torque: First argument must be a string.");
 	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("gait3d_pelvis213_torque: First argument must be a row vector.");
	command = mxArrayToString(prhs[0]);
	if(command == NULL) 
	  mexErrMsgTxt("gait3d_pelvis213_torque: First argument was not a command string.");
	  
	// initialize?
	if (strcmp(command, "Initialize") == 0) {
		const mxArray *P;									// pointer to the parameter struct
		printf("***************************************************\n");
		printf("*      MEX Function gait3d_pelvis213_torque       *\n");  
		printf("*       (c) 2013-2015 Orchard Kinetics LLC        *\n");
		printf("***************************************************\n");
		printf("Initializing...\n");
		if (nrhs != 2)
			mexErrMsgTxt("gait3d_pelvis213_torque:Initialize: Two input arguments required.");
			
		// get a pointer P to the second argument and check that it is a 1x1 struct
		P = prhs[1];
		if ((mxGetM(P) != 1) || (mxGetN(P) != 1)) {
			mexErrMsgTxt("gait3d_pelvis213_torque:Dynamics: Incorrect size for parameters, must be 1 x 1.");
		}
		if (!mxIsStruct(P)) {
			mexErrMsgTxt("gait3d_pelvis213_torque:Initialize: Incorrect type for parameters, must be struct.");	
		}
		
		// from P, extract the parameters needed by the Autolev generated code, and store them in the C struct "parameters"
		// see gait3d_pelvis213_torque.h for information about the fields in the parameters.struct
		parameters.gravity_x = extract(P, "gravity", 1);
		parameters.gravity_y = extract(P, "gravity", 2);
		parameters.gravity_z = extract(P, "gravity", 3);
		
		// parameters of the air drag model
		parameters.airdrag 	= extract(P, "drag_coefficient", 1);
		parameters.wind 	= extract(P, "wind_speed", 1);
		
		// Check that the nDofs in the model struct is the same as the constant NDOF we use internally here
		if (NDOF - extract(P, "nDofs", 1) != 0) {
			mexErrMsgTxt("gait3d_pelvis213_torque:Initialize: nDofs is not consistent with this version of the MEX function.");		
		}
		
		// Check that the Nsegments in the model struct is one plus the NFKSEG we use internally here (Opensim has ground as a segment)
		if (NFKSEG + 1 - extract(P, "nSegments", 1) != 0) {
			mexErrMsgTxt("gait3d_pelvis213_torque:Initialize: nSegments is not consistent with this version of the MEX function.");		
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
		parameters.Rknee_x5 	= extract2(P,"joints",3,"t1_coefs",5);		
		parameters.Rknee_y1 	= extract2(P,"joints",3,"t2_coefs",1);		
		parameters.Rknee_y2 	= extract2(P,"joints",3,"t2_coefs",2);		
		parameters.Rknee_y3 	= extract2(P,"joints",3,"t2_coefs",3);		
		parameters.Rknee_y4 	= extract2(P,"joints",3,"t2_coefs",4);		
		parameters.Rknee_y5 	= extract2(P,"joints",3,"t2_coefs",5);		
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
		parameters.Lknee_x5 	= extract2(P,"joints",8,"t1_coefs",5);		
		parameters.Lknee_y1 	= extract2(P,"joints",8,"t2_coefs",1);		
		parameters.Lknee_y2 	= extract2(P,"joints",8,"t2_coefs",2);		
		parameters.Lknee_y3 	= extract2(P,"joints",8,"t2_coefs",3);		
		parameters.Lknee_y4 	= extract2(P,"joints",8,"t2_coefs",4);		
		parameters.Lknee_y5 	= extract2(P,"joints",8,"t2_coefs",5);		
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
		
		// from P, extract the contact elements and their properties
		{
			mxArray *field;
			int i;
			ncontacts = (int) extract(P, "nCPs", 1);
			if (ncontacts > MAXCONTACTS)
				mexErrMsgTxt("gait3d_pelvis213_torque:Initialize: too many contact elements.");
			if ((field = mxGetField(P,0,"CPs") ) == NULL) {
				mexErrMsgTxt("gait3d_pelvis213_torque:Initialize: no field named CPs found in model.");
			}
			if (!mxIsCell(field)) {
				mexErrMsgTxt("gait3d_pelvis213_torque: Initialize: contacts field is not a cell array.");
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
					mexErrMsgTxt("gait3d_pelvis213_torque: Initialize: contact point attached to a non-foot segment.");
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
                joints[i].qmin_passiveMoment = extract2(P,"dofs",i+1,"range_passiveMoment",1);
				joints[i].qmax_passiveMoment = extract2(P,"dofs",i+1,"range_passiveMoment",2);
			}
		}
				
		// set the initialized flag and compute some constants
		initialized = 1;
		nstates = 2*NDOF + NCVAR*ncontacts;
		if (nstates > MAXSTATES)
			mexErrMsgTxt("gait3d_pelvis213_torque:Initialize: too many states.");
		nf = 2*NDOF + NCF*ncontacts;					// number of elements in f(x,xdot,u,M)
		
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
			
			// for the multibody dynamics, lb=ub=0 (equality constraints f=0)
			for (i=0; i<2*NDOF; i++) {
				fmin[j]   = 0.0;
				fmax[j++] = 0.0;
			}

			// for each contact element, 9 equality constraints (f=0) and 4 inequality constraints (f>0)
			for (i=0; i<ncontacts; i++) {
				for (k=0; k<NCEQ; k++) {
					fmin[j]   = 0.0;
					fmax[j++] = 0.0;
				}
				for (k=0; k<NCINEQ; k++) {
					fmin[j]   = 0.0;
					// fmax[j++] = mxGetInf();
					fmax[j++] = 10;
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
		
		// Joint variables
		double M[NDOF];							// total actuation to be supplied to Autolev generated code
		double dM_dq[NDOF][NDOF];
		double dM_dqd[NDOF];					// is a diagonal matrix (no coupling between joints), so no need for full storage
		
		// GRF variables
		double G[NGRF];								// 3D ground reaction force/moment for the segments that have contact elements
		double dG_dxc[NGRF][MAXCONTACTS*NCVAR];		// derivatives of G with respect to contact variables in x

		// some logical variables related to the user's request
		int derivatives 	= (nlhs > 1);						// will be true if user requested derivatives;
		int cmdFkin 		= (strcmp(command, "Fkin") == 0);
		int cmdDynamics 	= (strcmp(command, "Dynamics") == 0);
		int cmdGRF 			= (strcmp(command, "GRF") == 0);
		int cmdJointmoments	= (strcmp(command, "Jointmoments") == 0);
        int cmdMuscleforces	= (strcmp(command, "Muscleforces") == 0);
        int cmdMuscleCEforces = (strcmp(command, "MuscleCEforces") == 0);
		int cmdMuscleCEpower= (strcmp(command, "MuscleCEpower") == 0);

		// Give error message if model was not initialized
		if (!initialized) {
			printf("gait3d_pelvis213_torque: command given: %s\n", command);
			mexErrMsgTxt("gait3d_pelvis213_torque: model was not initialized.");
		}

		// get all the inputs that are needed, depending on what the command is
		if (cmdDynamics) {
			if (nrhs < 4 || nrhs > 6) {
				mexErrMsgTxt("gait3d_pelvis213_torque:Dynamics: Must have 4,5 or 6 inputs.");			
			}
		
			// get x, xdot, u
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213_torque:Dynamics: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);
			
			if ((!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != nstates) || (mxGetN(prhs[2]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213_torque:Dynamics: Incorrect input xdot.");
			}
			xdot = mxGetPr(prhs[2]);
			
			u = mxGetPr(prhs[3]); // We actually do not use u.
			
			// get optional input M (additional actuation)
			if (nrhs > 4) {
				if ((!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetM(prhs[4]) != NDOF) || (mxGetN(prhs[4]) != 1)) {
					mexErrMsgTxt("gait3d_pelvis213_torque:Dynamics Incorrect input for M.");
				}
				Mextra = mxGetPr(prhs[4]);
			}
            
            // get G as input 
			if (nrhs > 5) {
				if ((!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetM(prhs[5]) != NGRF) || (mxGetN(prhs[5]) != 1)) {
					mexErrMsgTxt("gait3d_pelvis213_torque:Dynamics Incorrect input for G.");
				}
				Gextra = mxGetPr(prhs[5]);
			}
            
			
			// generate pointers to q, qd, and qdd
			q = &x[0];
			qd = &x[NDOF];	
			qdd = &xdot[NDOF];
			
		}
        
        if (cmdMuscleforces || cmdMuscleCEforces || cmdMuscleCEpower) {
			mexErrMsgTxt("gait3d_pelvis213: The model does not contain muscles. Hence, this function is not supported.");
		}

		if (cmdGRF) {
			if (nrhs != 2) {
				mexErrMsgTxt("gait3d_pelvis213_torque:GRF: Must have exactly 2 inputs.");			
			}
		
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213_torque:GRF: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);	
		}

		if (cmdJointmoments) {
			if (nrhs < 2 || nrhs > 3) {
				mexErrMsgTxt("gait3d_pelvis213_torque:jointmoments: Must have 2 or 3 inputs.");			
			}
			// get x
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != nstates) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213_torque:Dynamics: Incorrect input for x.");
			}
			x = mxGetPr(prhs[1]);	
			xdot = zeros;				// pointer to zeros, to avoid crash when xdot is not initialized
			
			// generate pointers to q, qd
			q = &x[0];
			qd = &x[NDOF];	
            
            if (nrhs > 2) {
                if ((!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != NDOF) || (mxGetN(prhs[2]) != 1)) {
                    mexErrMsgTxt("gait3d_pelvis213_torque:Jointmoments: Incorrect input for M.");
                }
                Mextra = mxGetPr(prhs[2]);
            }
            
		}
		
		if (cmdFkin) {
			int i;
			
			if (nrhs < 2) {
				mexErrMsgTxt("gait3d_pelvis213_torque:Fkin: must have at least one input.");			
			}
			
			// get q
			if ((!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != NDOF) || (mxGetN(prhs[1]) != 1)) {
				mexErrMsgTxt("gait3d_pelvis213_torque:Fkin Incorrect input for q.");
			}
			q = mxGetPr(prhs[1]);
			
			// case 1: 1 or 2 outputs requested
			if (nlhs < 3) {
				if (nrhs != 2) {
					mexErrMsgTxt("gait3d_pelvis213_torque:Fkin: Must have exactly 2 inputs when 2 outputs are requested.");			
				}
				qd = zeros;		// velocities are not needed
			}
			else if (nlhs < 4) {
				if (nrhs != 3) {
					mexErrMsgTxt("gait3d_pelvis213_torque:Fkin: Must have exactly 3 inputs when 3 outputs are requested.");			
				}
				// get qd (dq/dt)
				if ((!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != NDOF) || (mxGetN(prhs[2]) != 1)) {
					mexErrMsgTxt("gait3d_pelvis213_torque:Fkin Incorrect input for qd.");
				}
				qd = mxGetPr(prhs[2]);
			}
			else {
				mexErrMsgTxt("gait3d_pelvis213_torque:Fkin: Too many outputs requested.");			
			}

			// set some other things to zero before we call Autolev generated code
			qdd = zeros;
			for (i=0; i<NDOF; i++) M[i] = 0.0;
			for (i=0; i<NGRF; i++) G[i] = 0.0;
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
					dM_dqd[i] = 0.0;
				}	
				// internal moments due to passive limits,
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
					
				}
				dM_dx_jcs[NDOF] = NdM_dx;						// store final element number

				if (Nalloc_dM_dx < NdM_dx) {
					mexErrMsgTxt("gait3d_pelvis213_torquec: Insufficient memory allocation for dM_dx.");	
				}
				#if JACOBIANSIZING
					printf("dM_dx: allocated %d, actual %d\n", Nalloc_dM_dx, NdM_dx);
				#endif

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
				int jc = NCVAR*i;					// index within the contact variables in x for this contact element: : Fx,Fy,Fz,xc,yc,zc,xf,yf,zf
				int j = 2*NDOF + jc;	// index within x for variables of this contact element: : Fx,Fy,Fz,xc,yc,zc,xf,yf,zf
				
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
								dgrf_dx_irs[Ndgrf_dx] = 2*NDOF + k;	// store row number of this matrix element
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
						mexErrMsgTxt("gait3d_pelvis213_torque: Incorrect memory allocation for dgrf_dx.");	
					}
				}
			}
					
			return;	
		}	

		// For dynamics : call the C function that was generated by Autolev
		// TODO: we can create separate functions for results with and without Jacobians
		if (cmdDynamics) {
			if (derivatives) {
				gait3d_pelvis213_torque_al(&parameters, q, qd, qdd, G, zero, dz_dq, dz_dqd, dz_dqdd, dz_dG,
					fk, dfk_dq, fkdot, dfkdot_dq);
            } else {
				gait3d_pelvis213_torque_NoDer_al(&parameters, q, qd, qdd, G, zero,
					fk, dfk_dq, fkdot, dfkdot_dq); // derivatives of fk could also be removed in future
			}
		}

		// For forward kinematics: call the C function that was generated by Autolev
		// TODO: we can create separate functions for results with and without Jacobians
		if (cmdFkin) {
			gait3d_pelvis213_torque_FK_al(&parameters, q, qd, qdd, fk, dfk_dq, fkdot, dfkdot_dq);
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
                Nalloc_df_du = 0;
				plhs[3] = mxCreateSparse(0, nf, Nalloc_df_du, mxREAL);
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
					
					
					// derivatives of ZERO[i] with respect to ground contact forces should go into rows 2*NDOF+1 to end
					for (k=0; k<NCVAR*ncontacts; k++) {
						int kk;
						double tmp = 0.0;
						// add the contributions dz/dG * dG/dx
						for (kk=0; kk<NGRF; kk++) {
							tmp = tmp + dz_dG[i][kk]*dG_dxc[kk][k];
						}
						if (fabs(tmp) > 1e-10) {
							df_dx_irs[Ndf_dx] = 2*NDOF + k;		    // store row number of this matrix element
							df_dx[Ndf_dx++] = tmp;					// store the value of this matrix element
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

			// then there are NCF elements of f for each contact point
			for (i=0; i<ncontacts; i++) {
				double df_dfk[NCF][12];
				double df_dfkdot[NCF][12];
				double df_dxc[NCF][NCVAR];
				double df_dxcdot[NCF][NCVAR];
				
				// index to the NCVAR contact variables of this contact point
				int ixc = 2*NDOF + NCVAR*i;			

				// index to the NCF contact residuals of this contact point
				int ifc = 2*NDOF + NCF*i;			

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
					mexErrMsgTxt("gait3d_pelvis213_torque: Incorrect memory allocation for df_dx.");	
				}
				if (Nalloc_df_dxdot < Ndf_dxdot) {
					mexErrMsgTxt("gait3d_pelvis213_torque: Incorrect memory allocation for df_dxdot.");	
				}
				if (Nalloc_df_du < Ndf_du) {
					mexErrMsgTxt("gait3d_pelvis213_torque: Incorrect memory allocation for df_du.");		
				}
				if (nlhs > 4) if (Nalloc_df_dM < Ndf_dM) {
						mexErrMsgTxt("gait3d_pelvis213_torque: Incorrect memory allocation for df_dM.");	
				}
                if (nlhs > 5) if (Nalloc_df_dG < Ndf_dG) {
						mexErrMsgTxt("gait3d_pelvis213_torque: Incorrect memory allocation for df_dG.");	
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
		mexErrMsgTxt("gait3d_pelvis213_torque: Command not recognized.");
	}	
		

}
