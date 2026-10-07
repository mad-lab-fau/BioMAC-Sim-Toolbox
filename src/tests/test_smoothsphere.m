% test_smoothsphere.m
% Verification script for Gait2d_osim_smoothsphere dynamics derivatives

addpath(genpath('/home/rzlin/ri94mihu/phd/BiomechPriorVAE/BioMAC-Sim-Toolbox/'));

% Instantiate model
model = Gait2d_osim_smoothsphere('gait2d.osim');

% Generate random states and controls
rng(1);
nstates = model.nStates;
ncontrols = model.nControls;

% Neutral state with some random perturbation
x = model.states.xneutral + 0.02 * randn(nstates, 1);
xdot = 0.05 * randn(nstates, 1);
u = 0.5 * rand(ncontrols, 1);

% Evaluate analytical derivatives
[f, dfdx_anal, dfdxdot_anal, dfdu_anal] = model.getDynamics(x, xdot, u);

% Convert analytical transposes to standard Jacobians
dfdx_anal = dfdx_anal';
dfdxdot_anal = dfdxdot_anal';
dfdu_anal = dfdu_anal';

% Finite differences
dh = 1e-7;
nconstraints = length(f);

% dfdx numerical
dfdx_num = zeros(nconstraints, nstates);
for i = 1:nstates
    x_plus = x;
    x_plus(i) = x_plus(i) + dh;
    f_plus = model.getDynamics(x_plus, xdot, u);
    dfdx_num(:, i) = (f_plus - f) / dh;
end

% dfdxdot numerical
dfdxdot_num = zeros(nconstraints, nstates);
for i = 1:nstates
    xdot_plus = xdot;
    xdot_plus(i) = xdot_plus(i) + dh;
    f_plus = model.getDynamics(x, xdot_plus, u);
    dfdxdot_num(:, i) = (f_plus - f) / dh;
end

% dfdu numerical
dfdu_num = zeros(nconstraints, ncontrols);
for i = 1:ncontrols
    u_plus = u;
    u_plus(i) = u_plus(i) + dh;
    f_plus = model.getDynamics(x, xdot, u_plus);
    dfdu_num(:, i) = (f_plus - f) / dh;
end

% Compare
err_dfdx = max(max(abs(dfdx_anal - dfdx_num)));
err_dfdxdot = max(max(abs(dfdxdot_anal - dfdxdot_num)));
err_dfdu = max(max(abs(dfdu_anal - dfdu_num)));

fprintf('Max dfdx derivative error: %g\n', err_dfdx);
fprintf('Max dfdxdot derivative error: %g\n', err_dfdxdot);
fprintf('Max dfdu derivative error: %g\n', err_dfdu);

if err_dfdx < 2e-3 && err_dfdxdot < 2e-3 && err_dfdu < 2e-3
    disp('DERIVATIVE TEST PASSED!');
else
    error('DERIVATIVE TEST FAILED!');
end
