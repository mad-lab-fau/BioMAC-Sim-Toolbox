%======================================================================
%> @file @Model/getEratec_Minetti.m
%> @brief Model function to calculate the energy rate of a single time
%> step with the continuous version of Minetti and Alexander's model
%> @details
%> Details: Model::getEratec_Minetti()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using Minetti's model. No discsontinuities here, so it is the same
%> as the original model.
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> Minetti's model. Based on: Minetti & Alexander, J Theor Biol, 1997. Uses
%> muscle-level work instead of torque-level
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
%> @param  obj         Model object
%> @param act          activation level of the muscle (between 0 and 1)
%> @param v_ce         current velocity of the contractile element
%>
%> The output of getEratec_Minetti.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval             Edot is the energy expenditure during the current time step
%> @retval             dEdot is the gradient of Edot with respect to all states and
%>                     inputs
%======================================================================

function [Edot,dEdot] = getEratec_Minetti(obj, act, v_ce)

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
fmax     = obj.muscles.fmax;    % Max muscle force
if ismember('v_max', obj.muscles.Properties.VariableNames)
    v_max = obj.muscles.v_max;  % Maximum shortening velocity
else
    v_max = obj.muscles.vmax;  % Maximum shortening velocity
end
v_max = v_max.*l_ceopt;

v_ce = -v_ce.*l_ceopt./v_max; %shortening positive

phi_num = (0.054+0.506*v_ce+2.46*v_ce.^2);
phi_den = (1-1.13*v_ce+12.8*v_ce.^2-1.64*v_ce.^3);
phi = phi_num./phi_den;

Edot = act.*fmax.*v_max.*phi;

if nargout > 1

    dphi_dpnum = 1./phi_den;
    dphi_dpden = -phi_num./phi_den.^2;

    dpnum_dvce = 0.506+2*2.46*v_ce;
    dpden_dvce = -1.13+12.8*2*v_ce-1.64*3*v_ce.^2;

    dphi_dvce = dphi_dpnum.*dpnum_dvce + dphi_dpden.*dpden_dvce;

    dEdot_dxi = zeros(obj.nMus,obj.nDofs*2);
    dEdot_dact = fmax.*v_max.*phi;
    dEdot_dvce = -act.*dphi_dvce.*fmax.*v_max.*l_ceopt./v_max;

    dEdot = [dEdot_dxi zeros(obj.nMus,1) dEdot_dact dEdot_dvce zeros(obj.nMus,1)];
end