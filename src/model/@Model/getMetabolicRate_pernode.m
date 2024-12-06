%======================================================================
%> @file @Model/getMetabolicRate_pernode.m
%> @brief Model function to calculate metabolic cost of a movement
%> @details
%> Details: Model::getMetabolicRate_pernode()
%>
%> @author Anne Koelewijn, Marlies Nitschke
%> @date March 23, 2015
%======================================================================

%======================================================================
%> @brief Model function to calculate metabolic cost of a movement
%> @details
%> Function to calculate the metabolic cost of a movement. Seven different
%> models are currently implemented
%>
%> @param  obj         Model object
%> @param  x           Double vector: States of the current node
%> @param  xdot        Double vector: Derivatives of states
%> @param  u           Double vector: Controls of the current node
%> @param  t_stim      Double vector: how long each muscle was already stimulated
%> @param  name        (optional) String: name of the model that is being
%>                     used. Options are umberger, lichtwark, bhargava,
%>                     margaria, minetti, houdijk, and uchida
%> @param  getCont     (optional) Boolean: If true, get the continuous version 
%>                     which is needed if we use the output for simulation. (default: 0)
%> @param  epsilon     (optional) double parameter: Measure of nonlinearity of model, 
%>                     which is only used in the continuous model version. Default 0.
%> @param  exponent    (optional) integer: if metabolic cost is the objective, this number 
%>                     allows to minimize the square, cube, or nth-power of metabolic cost. Default is 1
%> @retval Edot        Double vector: Energy expenditure during one time step (energy rate) in W (Model.nMus x 1)
%> @retval dEdotdx     Double vector: Derivative of Edot w.r.t to x (Model.nStates x 1)
%> @retval dEdotdu     Double vector: Derivative of Edot w.r.t to u (Model.nControls x 1)
%> @retval dEdotdT     Double: Derivative of Edot w.r.t to T
%======================================================================
function  [Edot, dEdotdx, dEdotdu, dEdotdxdot] = getMetabolicRate_pernode(obj, x, xdot, u, t_stim, name, getCont, epsilon, exponent) % 

% Check whether we should return the continuous version which is needed if we use the output for simulation 
if nargin < 9
    exponent = 1;
end

if nargin < 8
    getCont = 0;
    epsilon = 0;
end
if nargin < 7
    name = 'umberger';
end

% Extract variables from x and u
stim  = u(obj.extractControl('u'));   % neural excitation is stimulus of muscles
act   = x(obj.extractState('a'));     % activation of the muscles
l_ce  = x(obj.extractState('s'));     % length of contractile element in percentage of lceopt (current node)
v_ce  = xdot(obj.extractState('s'));  % Get velocity of contractile elements
nDofs = obj.nDofs;


% Calculate energy rate in Watts for all muscles at once
if ~getCont
    % If we want to get the energy rate after simulation (postprocessing),
    % we can use the discontinuous version (= "original" version).
    % Jacobian of energy rate
    % Get muscle forces of contractile element (CE) and the jacobians of it
    F_ce = obj.getMuscleCEforces(x, xdot);
    Edot = 0;
    dEdot = nan(obj.nMus, 2*obj.nDofs+5); % During postprocessing we are not interessted in dEdot. => We don't compute it
    switch name %Switch between different models
        case 'umberger'
            Edot = obj.getErate_Umberger(F_ce, stim, act, l_ce, v_ce);
        case 'bhargava'
            Edot = obj.getErate_bhargava(F_ce, stim, act, l_ce, v_ce, t_stim);
        case 'houdijk'
            Edot = obj.getErate_Houdijk(F_ce, act, l_ce, v_ce);
        case 'lichtwark'
            Edot = obj.getErate_Lichtwark(F_ce, act, l_ce, v_ce, t_stim);
        case 'margaria'
            Edot = obj.getErate_margaria(F_ce, v_ce);
        case 'minetti'
            Edot = obj.getErate_Minetti(act, v_ce);
        case 'uchida'
            Edot = obj.getErate_uchida(F_ce, stim,act, l_ce, v_ce);
    end
elseif nargout ==1 
    [F_ce, Fce_dx, Fce_dxdot] = obj.getMuscleCEforces(x, xdot);
    Edot = 0;
    switch lower(name) %Switch between different models
        case 'umberger'
            Edot = obj.getEratec_umberger(F_ce, stim, act, l_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'bhargava'
            Edot = obj.getEratec_bhargava(F_ce, stim, act, l_ce, v_ce, Fce_dx, Fce_dxdot,epsilon);
        case 'houdijk'
            Edot = obj.getEratec_Houdijk(F_ce, act, l_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'lichtwark'
            Edot = obj.getEratec_Lichtwark(F_ce, act, l_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'margaria'
            Edot = obj.getEratec_margaria(F_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'minetti'
            Edot = obj.getEratec_Minetti(act, v_ce);
    end
else
    % If we want to use the energy rate as objective during simulation,
    % we need to make the function continuous.
    % Get muscle forces of contractile element (CE) and the jacobians of it
    [F_ce, Fce_dx, Fce_dxdot] = obj.getMuscleCEforces(x, xdot);
    Edot = 0;
    dEdot = nan(obj.nMus, 2*obj.nDofs+5);
    switch lower(name) %Switch between different models
        case 'umberger'
            [Edot, dEdot] = obj.getEratec_umberger(F_ce, stim, act, l_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'bhargava'
            [Edot, dEdot] = obj.getEratec_bhargava(F_ce, stim, act, l_ce, v_ce, Fce_dx, Fce_dxdot,epsilon);
        case 'houdijk'
            [Edot, dEdot] = obj.getEratec_Houdijk(F_ce, act, l_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'lichtwark'
            [Edot, dEdot] = obj.getEratec_Lichtwark(F_ce, act, l_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'margaria'
            [Edot, dEdot] = obj.getEratec_margaria(F_ce, v_ce, Fce_dx, Fce_dxdot, epsilon);
        case 'minetti'
            [Edot, dEdot] = obj.getEratec_Minetti(act, v_ce);
    end
    % Put the derivatives from this time step into matrix form
    [dEdotdx,dEdotdxdot] = deal(zeros(size(x)));
    dEdotdu = zeros(size(u));
    % Derivatives w.r.t. q
    dEdotdx(obj.extractState('q'))    = sum(exponent*Edot.^(exponent-1).*dEdot(:, 1:nDofs))';
    % Derivatives w.r.t. qdot
    dEdotdx(obj.extractState('qdot')) = sum(exponent*Edot.^(exponent-1).*dEdot(:, nDofs+(1:nDofs)))';
    % Derivatives w.r.t. s (length of contractile element)
    dEdotdx(obj.extractState('s'))    = exponent*Edot.^(exponent-1).*dEdot(:, nDofs*2+1);
    % Derivatives w.r.t. a (muscle activation)
    dEdotdx(obj.extractState('a'))    = exponent*Edot.^(exponent-1).*dEdot(:, nDofs*2+2);
    % Derivatives w.r.t. sdot (velocity of contractile element)
    dEdotdxdot(obj.extractState('s'))    = exponent*Edot.^(exponent-1).*dEdot(:, nDofs*2+3);
    % Derivatives w.r.t. u (neural excitation = stimulation of muscles)
    dEdotdu(obj.extractControl('u'))  = exponent*Edot.^(exponent-1).*dEdot(:, nDofs*2+4);

end
end
