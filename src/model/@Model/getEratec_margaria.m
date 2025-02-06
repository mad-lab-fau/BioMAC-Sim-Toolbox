%======================================================================
%> @file @Model/getEratec_margaria.m
%> @brief Model function to calculate the energy rate of a single time
%> step with the continuous version of Margaria's model
%> @details
%> Details: Model::getEratec_margaria()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using a continuous version of Margaria's model
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> a continuous version of Margaria's model: 25 % efficient during shortening, 120% during lengthening.
%> Margaria, Int Z Angew Physiol Einschl Arbeitsphysiol, 1968
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
%> @param  obj         Model object
%> @param F_ce         Force in the contractile element
%> @param v_ce         Normalized velocity of the contractile element
%> @param dFdx         Derivative of CE forces with respect to the state
%> @param dFdxdot      Derivative of CE forces with respect to the state
%>                     derivative
%> @param epsilon      Measure of the nonlinearity of the continuous functions
%>
%> The output of getEratec_margaria.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot is the energy expenditure during the current time step
%> @retval dEdot is the gradient of Edot with respect to all states and
%> inputs
%======================================================================

function [Edot, dEdot] = getEratec_margaria(obj, F_ce, v_ce, dFdx,dFdxdot,epsilon)

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length

v_ce_l = 1/2*(v_ce+sqrt(v_ce.^2+epsilon^2)); %lengthening velocity
v_ce_s = 1/2*(v_ce-sqrt((-v_ce).^2+epsilon^2)); %shortening velocity

Edot = (-F_ce/0.25.*v_ce_s+F_ce/1.2.*v_ce_l).*l_ceopt;

if nargout > 1
    dFce_dxi = dFdx([obj.extractState('q'); obj.extractState('qdot')] ,:)';
    dFce_dlce = diag(dFdx(obj.extractState('s'),:));
    dFce_dact = diag(dFdx(obj.extractState('a'),:));
    dFce_dvce = diag(dFdxdot(obj.extractState('s'),:));

    dvcel_dvce = 1/2.*(1+v_ce./sqrt(v_ce.^2+epsilon^2));
    dvces_dvce = 1/2.*(1-v_ce./sqrt((-v_ce).^2+epsilon^2));

    dEdot_dxi = bsxfun(@times,dFce_dxi,-v_ce_s/0.25)+bsxfun(@times,dFce_dxi, v_ce_l/1.2);
    dEdot_dact = -dFce_dact/0.25.*v_ce_s+dFce_dact/1.2.*v_ce_l;
    dEdot_dlce = -dFce_dlce/0.25.*v_ce_s+dFce_dlce/1.2.*v_ce_l;
    dEdot_dvce = -F_ce/0.25.*dvces_dvce+F_ce/1.2.*dvcel_dvce-dFce_dvce/0.25.*v_ce_s+dFce_dvce/1.2.*v_ce_l;

    dEdot = [dEdot_dxi dEdot_dlce dEdot_dact dEdot_dvce zeros(obj.nMus,1)].*l_ceopt;
end