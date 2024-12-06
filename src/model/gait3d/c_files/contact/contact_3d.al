% Autolev source file for the contact model

% The autolev code will generate C code for the contact model, in the
% form: fmin <= f(x,xdot) <= fmax, 
% where x are the 6 variables of the contact model:
%   Fx,Fy,Fz	force acting on the contact point, in global coordinate system
%   xc,yc,zc    global coordinates of the contact point
% and f are the 6 contact equations coded in this file
% fmin = zeros(6,1);
% fmax = zeros(6,1);
% 
% After cleaning by clean_contact.c, the C function will look like this:
% void contact_al(contactprop* contact,										// parameters for one contact point
%	double fk[12], double fkdot[12], double x[6], double xdot[6],			// inputs for one contact point
%	double f[6], double df_dfk[6][12], double df_dfkdot[6][12],			// outputs
%	double df_dx[6][6], double df_dxdot[6][6]);							

PAUSE 0
AUTOZ ON
OVERWRITE ALL

% constants
constants contact__xp,contact__yp, contact__zp			% position of contact element on body segment

% input arrays
variables fk{12}, fkdot{12}		% forward kinematic variables for body segment: XO,YO,ZO,R11,R12,R13,R21,R22,R23,R31,R32,R33
variables x{6},xdot{6}			% state variables of the contact element: Fx,Fy,Fz,xc,yc,zc

% rename the input variables for more readable code
variables Fx,Fy,Fz,xc,yc,zc
% force acting on the contact point, in global coordinate system
Fx = x1
Fy = x2
Fz = x3
% global coordinates of the contact point, and velocity
xc = x4
yc = x5
zc = x6
xcdot = xdot4
ycdot = xdot5
zcdot = xdot6

% transform force from global to local coordinates
% Floc = R' * Fglb
variables Fxloc, Fyloc, Fzloc
Fxloc = fk4*Fx + fk7*Fy + fk10*Fz
Fyloc = fk5*Fx + fk8*Fy + fk11*Fz
Fzloc = fk6*Fx + fk9*Fy + fk12*Fz

%========================================
% equations 1-3: contact points
%========================================
% transform contact point coordinates xc,yc,zc from global to local coordinates
% Ploc = R' * (Pglb - p)
% p = fk1,fk2,fk3 and R = [fk4,fk5,fk6;fk7,fk8,fk9;fk10,fk11,fk12]
variables xcloc, ycloc, zcloc
xcloc = fk4*(xc-fk1) + fk7*(yc-fk2) + fk10*(zc-fk3)
ycloc = fk5*(xc-fk1) + fk8*(yc-fk2) + fk11*(zc-fk3)
zcloc = fk6*(xc-fk1) + fk9*(yc-fk2) + fk12*(zc-fk3)

% constrain the location of the contact points and contact states
f1 = xcloc - contact__xp				
f2 = ycloc - contact__yp
f3 = zcloc - contact__zp

%=========================================
% equations 4-6: contact equations
%=========================================

% simple penetration-based contact model, stiffness 100 BW per meter.
f4 = (Fy + 100 * 0.5 * (yc - sqrt(yc^2 + 0.001^2)) * (1 - 0.75*ycdot))

% simple friction model
f5 = (Fx + Fy * xcdot / sqrt(xcdot^2 + 1e-4))
f6 = (Fz + Fy * zcdot / sqrt(zcdot^2 + 1e-4))


% combine all residuals into one vector f
f = [f1;f2;f3;f4;f5;f6]

% generate expressions for the derivatives of f:
fk = [fk1,fk2,fk3,fk4,fk5,fk6,fk7,fk8,fk9,fk10,fk11,fk12];
fkdot = [fkdot1,fkdot2,fkdot3,fkdot4,fkdot5,fkdot6,fkdot7,fkdot8,fkdot9,fkdot10,fkdot11,fkdot12];
x = [x1,x2,x3,x4,x5,x6];
xdot = [xdot1,xdot2,xdot3,xdot4,xdot5,xdot6];
df_dfk 	= zee(d(f,fk))
df_dfkdot	= zee(d(f,fkdot))
df_dx 		= zee(d(f,x))
df_dxdot 	= zee(d(f,xdot))

encode f, df_dfk, df_dfkdot, df_dx, df_dxdot
code Algebraic() contact_3d_raw.c

EXIT


