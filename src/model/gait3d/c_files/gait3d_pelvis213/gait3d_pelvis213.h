// gait3d_pelvis213.h
// This file defines the data structure that holds model parameters for the gait3d_pelvis213 model

#define MAXMUSDOF 4		// maximum number of DOFs spanned by one muscle
#define MAXPOLTERMS 40	// maximum number of polynomial terms in a muscle path model
#define NAMELENGTH 40	// number of characters in muscle name

#include "gait3d_pelvis213_multibody.h"
#include "../contact/gait3d_contact.h"

// M_PI is known in gcc but not in Visual C++ 2008
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

// structure that holds muscle properties for one muscle
typedef struct {
    char name[NAMELENGTH];    				// Name of muscle 
	double Lceopt; 							// Optimal length of CE (m)
    double pennopt;         				// Pennation angle at optimal length of CE, in radians
	double width;							// Width of CE force-length relationship relative to Lceopt
	double Fmax;							// Maximal isometric force of CE (N)
	double Vmax;							// Max. contraction velocity in Lceopt/s
	double Tact, Tdeact;					// Activation and deactivation time constants
	double gmax;							// Maximum eccentric force
	double SEEslack;						// Slack length of the SEE (m)
	double PEEslack;						// Slack length of the PEE, relative to Lceopt
	double umax;							// Strain of SEE at Fmax load
	double krel;                            // Stiffness of PEE, force/Fmax at elongation of Width*Lceopt
	double HillA;							// Normalized Hill parameter a/Fmax for F-v relationship (usually 0.25)
	int    nmusdof;              			// Number of DOFs between origin and insertion 
    int    musdof[MAXMUSDOF];    			// List of DOFs between origin and insertion 
    int    npolterms;            			// Number of terms in polynomial 
    double  polcoeff[MAXPOLTERMS];        	// Polynomial coefficients 
    int    expon[MAXPOLTERMS][MAXMUSDOF];   // Polynomial exponents 
	int spans[NDOF];						// true/false flags for all DOFs to indicate if the muscle spans it
	
	// the following parameters are derived from other parameters during initialization
	double c3;								// Continuity parameter for eccentric force-velocity relationship
	double kPEE;							// Stiffness parameter of PEE, relative to Fmax/Lceopt^2
	double kSEE;							// Stiffness parameter of SEE, in Fmax/m^2	
} muscleprop;

// structure that holds parameters for one joint
typedef struct {
	double K1;						// linear stiffness (small), in N/rad
	double K2;						// quadratic stiffness, outside range of motion, in N/rad^2
	double B;						// damping, in Ns/rad
	double qneutral;		        // neutral position to compute passive joint moments
    double qmin_muscleMoment, qmax_muscleMoment;   // range of dof for computation of muscle moments
    double qmin_passiveMoment, qmax_passiveMoment; // range of dof for computation of passive moments
	
} jointprop;
