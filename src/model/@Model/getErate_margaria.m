%======================================================================
%> @file @Model/getErate_margaria.m
%> @brief Model function to calculate the energy rate of a single time
%> step with Margaria's model
%> @details
%> Details: Model::getErate_margaria()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using Margaria's model
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> Margaria's model: 25 % efficient during shortening, 120% during lengthening.
%> Margaria, Int Z Angew Physiol Einschl Arbeitsphysiol, 1968
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
%> @param  obj         Model object
%> @param F_ce         Force in the contractile element
%> @param v_ce         Normalized velocity of the contractile element
%>
%> The output of getErate_margaria.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot is the energy expenditure during the current time step
%======================================================================

function Edot = getErate_margaria(obj, F_ce, v_ce)

l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
W = F_ce.*v_ce.*l_ceopt;
Edot = -W/0.25.*(v_ce < 0)+W/1.2.*(v_ce>=0);

