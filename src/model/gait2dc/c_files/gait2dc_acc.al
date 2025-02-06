% gait2d_acc.al
%
% This is an input file for Autolev to generate linear accelerometer models
% for the gait2d model.
% Sensor signals Sx and Sy are in segment reference frame, and 
% modeled as a function of q,qdot,qdotdot,segment,xpos,ypos
%
% Sx = a1 + a2*xpos + a3*ypos
% Sy = a4 + a5*xpos + a6*ypos
%
% Autolev will generate code for the coefficients a1..a6, for all seven segments
% in a vector a.
% And we also will need the Jacobians da/dq, da/dqdot, da/dqdotdot
%
% Author: Ton van den Bogert
% Last revised: 03/24/2011

RUN gait2dc_kin.al

% gravity is needed in the accelerometer model
Constants par__gravity					% Local gravitational acceleration value

%--------------------------------------------------------------------------
% Create the accelerometer models for each segment, using the
% "accelerometer.a" macro which creates a six-element array a1..a6
% for a body segment.
%--------------------------------------------------------------------------
accelerometer(Trunk,Hip)		% generates accelerometer coefficients for Trunk, placement relative to Hip
accelerometer(RThigh,Hip)
accelerometer(RShank,RKnee)
accelerometer(RFoot,RAnkle)
accelerometer(LThigh,Hip)
accelerometer(LShank,LKnee)
accelerometer(LFoot,LAnkle)

% Put all of them into a single array with 42 elements
acc = [acc_Trunk; acc_RThigh; acc_RShank; acc_RFoot; acc_LThigh; acc_LShank; acc_LFoot]

% We also need Jacobians w.r.t. q, qd, qdd	
dacc_dq = D(acc,q); 
dacc_dqd = D(acc,qd); 
dacc_dqdd = D(acc,qdd); 
encode acc, dacc_dq, dacc_dqd, dacc_dqdd

% Generate C code
Code Algebraic() gait2dc_acc_al_raw.c

EXIT
