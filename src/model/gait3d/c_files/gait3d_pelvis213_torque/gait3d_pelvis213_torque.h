// gait3d_pelvis213_torque.h
// This file defines the data structure that holds model parameters for the gait3d_pelvis213_torque model

#include "gait3d_pelvis213_torque_multibody.h"
#include "../contact/gait3d_contact.h"

// M_PI is known in gcc but not in Visual C++ 2008
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


// structure that holds parameters for one joint
typedef struct {
	double K1;						// linear stiffness (small), in N/rad
	double K2;						// quadratic stiffness, outside range of motion, in N/rad^2
	double B;						// damping, in Ns/rad
	double qneutral;		        // neutral position to compute passive joint moments
    double qmin_passiveMoment, qmax_passiveMoment; // range of dof for computation of passive moments
	
} jointprop;
