// gait3d_smoothsphere_contact.h
// This file defines things needed in the 3D smoothsphere contact model

#define NCVAR 6						// number of state variables for each contact element
#define NCF 6               		// total number of f's for each contact element

typedef struct {
	int segment;			// Which body segment it is attached to (referring to forward kinematic segment list)
	int grfseg;				// Which GRF segment this contact point contributes to (referring to GRF array)
	double xp,yp,zp;		// Position of contact element on body segment
} contactprop;

// prototype for the C function for contact deformation
void contact_al(contactprop* contact,										// parameters for one contact point
	double fk[12], double fkdot[12], double x[NCVAR], double xdot[NCVAR],		// inputs for one contact point
	double f[NCF], double df_dfk[NCF][12], double df_dfkdot[NCF][12],			// outputs
	double df_dx[NCF][NCVAR], double df_dxdot[NCF][NCVAR]);							
