% Autolev source file for the contact model

% The autolev code will generate C code for the contact model, in the
% form: fmin <= f(x,xdot) <= fmax, 
% where x are the 4 variables of the contact model:
%   Fx,Fy    force acting on the contact point, in global coordinate system
%   xc,yc    global coordinates of the contact point
% and f are the 4 contact equations coded in this file
% fmin = zeros(4,1);
% fmax = zeros(4,1);
% 
% After cleaning by clean_contact.c, the C function will look like this:
% void contact_al(contactprop* contact,									// parameters for one contact point
%	double fk[6], double fkdot[6], double x[4], double xdot[4],			// inputs for one contact point
%	double f[4], double df_dfk[4][6], double df_dfkdot[4][6],			// outputs
%	double df_dx[4][4], double df_dxdot[4][4]);							

PAUSE 0
AUTOZ ON
OVERWRITE ALL

% constants
constants contact__xp,contact__yp		% position of contact element on body segment

% input arrays
variables fk{6}, fkdot{6}		% forward kinematic variables for body segment: XO,YO,R11,R12,R21,R22
variables x{4},xdot{4}			% state variables of the contact element: Fx,Fy,xc,yc

% rename the input variables for more readable code
variables Fx,Fy,xc,yc
% force acting on the contact point, in global coordinate system
Fx = x1
Fy = x2
% global coordinates of the contact point, and velocity
xc = x3
yc = x4
xcdot = xdot3
ycdot = xdot4

% transform force from global to local coordinates
% Floc = R' * Fglb
variables Fxloc, Fyloc, Fzloc
Fxloc = fk3*Fx + fk5*Fy
Fyloc = fk4*Fx + fk6*Fy

%========================================
% equations 1-3: contact points
%========================================
% transform contact point coordinates xc,yc,zc from global to local coordinates
% Ploc = R' * (Pglb - p)
% p = fk1,fk2 and R = [fk3,fk4;fk5,fk6]
variables xcloc, ycloc, zcloc
xcloc = fk3*(xc-fk1) + fk5*(yc-fk2)
ycloc = fk4*(xc-fk1) + fk6*(yc-fk2)

% constrain the location of the contact points and contact states
f1 = xcloc - contact__xp				
f2 = ycloc - contact__yp

%=========================================
% equations 4-6: contact equations
%=========================================

% simple penetration-based contact model, stiffness 100 BW per meter.
f3 = (Fy + 100 * 0.5 * (yc - sqrt(yc^2 + 0.001^2)) * (1 - 0.75*ycdot))

% simple friction model
f4 = (Fx + Fy * xcdot / sqrt(xcdot^2 + 1e-4))


% combine all residuals into one vector f
f = [f1;f2;f3;f4]

% generate expressions for the derivatives of f:
fk = [fk1,fk2,fk3,fk4,fk5,fk6];
fkdot = [fkdot1,fkdot2,fkdot3,fkdot4,fkdot5,fkdot6];
x = [x1,x2,x3,x4];
xdot = [xdot1,xdot2,xdot3,xdot4];
df_dfk 	= zee(d(f,fk))
df_dfkdot	= zee(d(f,fkdot))
df_dx 		= zee(d(f,x))
df_dxdot 	= zee(d(f,xdot))

encode f, df_dfk, df_dfkdot, df_dx, df_dxdot
code Algebraic() contact_2d_raw.c

EXIT


