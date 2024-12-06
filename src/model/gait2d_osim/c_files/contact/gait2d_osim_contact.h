// gait2d_osim_contact.h
// This file defines things needed in the contact model

#define NCVAR 4						// number of state variables for each contact element
#define NCEQ 4						// number of equality constraints (f=0) for each contact element
#define NCINEQ 0					// number of inequality constraints (f>0) for each contact element
#define NCF 4               		// total number of f's for each contact element

typedef struct {
	int segment;			// Which body segment it is attached to (referring to my forward kinematic segment list)
	int grfseg;				// Which GRF segment this contact point contributes to (referring to my GRF array)
	double xp,yp;		// Position of contact element on body segment
} contactprop;

// prototype for the Autolev C function for contact deformation
void contact_al(contactprop* contact,										// parameters for one contact point
	double fk[6], double fkdot[6], double x[NCVAR], double xdot[NCVAR],		// inputs for one contact point
	double f[NCF], double df_dfk[NCF][6], double df_dfkdot[NCF][6],			// outputs
	double df_dx[NCF][NCVAR], double df_dxdot[NCF][NCVAR]);							
