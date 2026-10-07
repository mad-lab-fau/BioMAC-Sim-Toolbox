% test_smoothsphere_3d.m
% Verification script for Gait3d_smoothsphere dynamics derivatives

addpath(genpath('/home/rzlin/ri94mihu/phd/BiomechPriorVAE/BioMAC-Sim-Toolbox/'));

% Instantiate model
model = Gait3d_smoothsphere('gait3d_pelvis213.osim');

% Generate random states and controls
rng(1);
nstates = model.nStates;
ncontrols = model.nControls;

% Neutral state with some random perturbation
x = model.states.xneutral + 0.01 * randn(nstates, 1);
xdot = 0.02 * randn(nstates, 1);
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

% Compare and print details
dfdx_anal = full(dfdx_anal);
dfdx_num = full(dfdx_num);
diff_dfdx = full(abs(dfdx_anal - dfdx_num));
[max_err, max_idx] = max(diff_dfdx(:));
[row, col] = ind2sub(size(diff_dfdx), max_idx);

fprintf('Max dfdx absolute error: %g at row %d, col %d\n', max_err, row, col);
fprintf('  Analytical value: %g\n', dfdx_anal(row, col));
fprintf('  Numerical value:  %g\n', dfdx_num(row, col));

% Show relative error if value is non-zero
if abs(dfdx_num(row, col)) > 1e-5
    rel_err = max_err / abs(dfdx_num(row, col));
    fprintf('  Relative error:   %g\n', rel_err);
end

% Check all elements that have error larger than 1e-3
[rows, cols] = find(diff_dfdx > 1e-3);
for k = 1:length(rows)
    r = rows(k);
    c = cols(k);
    fprintf('Large error at row %d, col %d: Anal = %g, Num = %g, AbsDiff = %g\n', r, c, dfdx_anal(r, c), dfdx_num(r, c), diff_dfdx(r, c));
end

err_dfdx = max_err;
err_dfdxdot = max(max(abs(dfdxdot_anal - dfdxdot_num)));
err_dfdu = max(max(abs(dfdu_anal - dfdu_num)));

if err_dfdx < 5e-3 && err_dfdxdot < 2e-3 && err_dfdu < 2e-3
    disp('3D DERIVATIVE TEST PASSED (with custom tolerance)!');
else
    error('3D DERIVATIVE TEST FAILED!');
end
