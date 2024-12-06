%======================================================================
%> @file @Model/getErate_Minetti.m
%> @brief Model function to calculate the energy rate of a single time
%> step with Minetti and Alexander's model
%> @details
%> Details: Model::getErate_Minetti()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%step using Minetti's model
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> Minetti's model. Based on: Minetti & Alexander, J Theor Biol, 1997. Uses
%> muscle-level work instead of torque-level
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
%> @param obj          Model object
%> @param act          activation level of the muscle (between 0 and 1)
%> @param v_ce         current velocity of the contractile element
%>
%> The output of getErate_Minetti is a double which is equal to the energy rate
%of the time step in W
%> @retval Edot is the energy expenditure during the current time step
%======================================================================

function Edot = getErate_Minetti(obj, act, v_ce)

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
fmax     = obj.muscles.fmax;    % Max muscle force

if ismember('v_max', obj.muscles.Properties.VariableNames)
    v_max = obj.muscles.v_max;  % Maximum shortening velocity
else
    v_max = obj.muscles.vmax;  % Maximum shortening velocity
end
v_max = v_max.*l_ceopt;

%Normalize v_ce to v_max:
v_ce = -v_ce.*l_ceopt./v_max; %shortening positive
phi = (0.054+0.506*v_ce+2.46*v_ce.^2)./(1-1.13*v_ce+12.8*v_ce.^2-1.64*v_ce.^3);

Edot = act.*fmax.*v_max.*phi;