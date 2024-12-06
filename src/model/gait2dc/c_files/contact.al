% deform.al: Autolev code to generate equations for the contact deformation with damping, for one contact point

PAUSE 0
AUTOZ ON
OVERWRITE ALL

% This autolev code will generate C code for the contact model.
% After cleaning by autolevclean.exe, the C function will look like this:

% void contact_al(contactprop* contact,										// parameters for one contact point
%	double fk[6], double fkdot[6], double xc[4], double xcdot[2],			// inputs for one contact point
%	double f[4], double df_dfk[4][6], double df_dfkdot[4][6],
%	double df_dxc[4][4], double df_dxcdot[4][2]);							// outputs

% constants
I = [1,0;0,1]
constants contact__kxx, contact__kxy, contact__kyy
constants contact__b, contact__c
constants contact__x,contact__y
constants contact__LambdaX, contact__LambdaY

% input arrays
variables fk{6}, fkdot{6}
variables xc{4},xcdot{2}

% some C code will be generated to convert the four input arrays into scalars with consecutive names

% friction equation for horizontal contact surface:			
% sqrt((cFy(vx-1) - Fx)^2 + Lambda) - sqrt((cFy(vx+1) - Fx)^2 + Lambda) - 2Fx = 0
variables d1,d2,sd1,sd2

scale1 = 10;
d1 = contact__c*xc4*(scale1*xcdot1-1) - xc3;
sd1 = sqrt(d1*d1 + contact__LambdaX);
d2 = contact__c*xc4*(scale1*xcdot1+1) - xc3;
sd2 = sqrt(d2*d2 + contact__LambdaX);
Ffriction = sd1 - sd2 - 2*xc3;

% smoothed Fischer-Burmeister equation for normal force: Fy + y - sqrt(Fy^2 + y^2 + LambdaY) = 0 
variables d

% Fy is in BW and y is in m, so we need to rescale to make them somewhat equal
scale2 = 200;
d = sqrt(xc4*xc4 + scale2*scale2*xc2*xc2 + contact__LambdaY);
Fnormal = xc4 + scale2*xc2 - d;

% extract input variables for deformation model
R = [fk3 , fk4 ; fk5 , fk6];
Rdot = [fkdot3 , fkdot4 ; fkdot5 , fkdot6];
p = [fk1; fk2];
pdot = [fkdot1; fkdot2];
pg = [xc1;xc2];					% global x,y of contact point
pgdot = [xcdot1;xcdot2];		% global x,y velocity of contact point
Fg = [xc3;xc4];					% global Fx,Fy at contact point, acting on segment

% extract model parameters from struct
K = [contact__kxx,contact__kxy; contact__kxy,contact__kyy];
B = contact__b
p0 = [contact__x;contact__y];

% transform global to local
Rinv = Transpose(R)
Fl = Rinv*Fg
ploc = Rinv*(pg-p)
v = Rinv*(pgdot-pdot) + transpose(Rdot)*(pg-p)
x = ploc - p0					% the amount of deformation along local x and y

% linear elasticity with additive damping
F = K*x + B*v

% the four model residuals
f = zee([Ffriction; Fnormal; Fl-F])

% generate expressions for the derivatives of f:
fk = [fk1,fk2,fk3,fk4,fk5,fk6];
fkdot = [fkdot1,fkdot2,fkdot3,fkdot4,fkdot5,fkdot6];
xc = [xc1,xc2,xc3,xc4];
xcdot = [xcdot1,xcdot2];
df_dfk 		= zee(d(f,fk))
df_dfkdot	= zee(d(f,fkdot))
df_dxc 		= zee(d(f,xc))
df_dxcdot 	= zee(d(f,xcdot))

encode f, df_dfk, df_dfkdot, df_dxc, df_dxcdot
code Algebraic() contact_al_raw.c

EXIT
