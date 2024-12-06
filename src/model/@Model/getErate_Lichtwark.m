%======================================================================
%> @file @Model/getErate_Lichtwark.m
%> @brief Model function to calculate the energy rate of a single time
%> step with Lichtwark and Wilson's model
%> @details
%> Details: Model::getErate_Lichtwark()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using a continuous version of Lichtwark & Wilson's model
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> Lichtwark & Wilson's model, described in the supplement of Lichtwark & Wilson,
%> J Biomech, 2007. The exception in gamma (line 68), which is according to
%> Lichtwark & Wilson, J Exp Biol, 2005, because results were better in the
%> 2019 comparison paper
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
%> @param  obj         Model object
%> @param F_ce         Force in the contractile element
%> @param act          Activation level of the muscle (between 0 and 1)
%> @param l_ce         Normalized length of the contractile element
%> @param v_ce         Normalized velocity of the contractile element
%> @param t_stim       Duration of stimulation of muscle
%>
%> The output of getErate_Lichtwark.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot        The energy expenditure during the current time step
%======================================================================

function Edot = getErate_Lichtwark(obj, F_ce, act, l_ce, v_ce,t_stim)

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
if ismember('width', obj.muscles.Properties.VariableNames)
    width    = obj.muscles.width;   % Width of hill curve
else
    width    = sqrt(obj.muscles.kactive);  % Width of hill curve
end
fmax     = obj.muscles.fmax;    % Max muscle force

if ismember('arel', obj.muscles.Properties.VariableNames)
    Ahill = obj.muscles.arel;
else
    Ahill = 0.25;  % Maximum shortening velocity
%     warning('Using default number for Hill parameter A')
end
if ismember('gmax', obj.muscles.Properties.VariableNames)
    gmax = obj.muscles.gmax;  % Maximum shortening velocity
else
    gmax = obj.muscles.flen;  % Maximum shortening velocity
end
if ismember('v_max', obj.muscles.Properties.VariableNames)
    v_max = obj.muscles.v_max;  % Maximum shortening velocity
else
    v_max = obj.muscles.vmax;  % Maximum shortening velocity
end

F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

% Force velocity relationship
c = v_max.*Ahill.*(gmax-1)./(Ahill+1);
g = ((v_max+v_ce)./(v_max-v_ce./Ahill)).*(v_ce < 0)+((gmax.*v_ce+c)./(v_ce+c)).*(v_ce >=0);

v_ce = -v_ce; %shortening positive

G = 4; % number described in supplement
gamma = 0.8*exp(-0.72*t_stim)+0.175*exp(-0.022*t_stim);
Hm = (gamma.*v_max/G^2).*(v_ce>0)+(0.3*gamma.*v_max/G^2+0.7*gamma.*v_max/G^2.*exp(-7*v_max.*(g-1))).*(v_ce<=0);
Hs = (v_ce/G).*(v_ce>0)+(-0.5*g.*v_ce).*(v_ce<=0);

%polo units -> isometric force and optimal fiber length
Hdot = 0.3*act.*Hm+0.7*act.*F_iso.*Hm+act.*F_iso.*Hs;

Hdot = Hdot.*l_ceopt.*fmax;

w = F_ce.*v_ce.*l_ceopt; %max(0,F_ce.*v_ce.*l_ceopt);%

Edot = Hdot+w;