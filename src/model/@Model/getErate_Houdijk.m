%======================================================================
%> @file @Model/getErate_Houdijk.m
%> @brief Model function to calculate the energy rate of a single time
%> step with Houdijk et al.'s model
%> @details
%> Details: Model::getErate_Houdijk()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%step using Houdijk et al.'s model
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> Houdijk et al.'s model. Houdijk et al., J Biomech, 2006.
%> The paper desription is not optimal, see comments
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
%>Inputs are muscle states of the specific muscle that is considered.
%> @param  obj         Model object
%> @param F_ce         Force in the contractile element
%> @param act          Activation level of the muscle (between 0 and 1)
%> @param l_ce         Normalized length of the contractile element
%> @param v_ce         Normalized velocity of the contractile element
%>
%> The output of getErate_Houdijk.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot        The energy expenditure during the current time step
%> @retval w           Muscle work during the current time step
%======================================================================


function [Edot,w] = getErate_Houdijk(obj, F_ce, act, l_ce, v_ce)

% Define variables
rho   = 1059.7;   % Muscle density
sigma = 250e3;    % Max muscle stress

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
FT       = obj.muscles.FT;      % Percentage of fast twitch fibers
if ismember('width', obj.muscles.Properties.VariableNames)
    width    = obj.muscles.width;   % Width of hill curve
else
    width    = sqrt(obj.muscles.kactive);  % Width of hill curve
end
fmax     = obj.muscles.fmax;    % Max muscle force

mmass = (fmax/sigma)*rho.*l_ceopt;

ST = 1-FT;
bargs = 0.45*24.4;
bargf = 0.35*150;
barhs = 24.4-bargs;
barhf = 150-bargf;

barg = bargs*ST+bargf*FT;
barh = barhs*ST+barhf*FT;
k3 = 6*ST+12*FT;
k4 = 8*ST+14*FT;
bara = 0.16*fmax.*ST+0.28*fmax.*FT;

nu = act.^2; %I (Anne) am 95% sure that I got this and the next equation from correspondence with Ross Miller. Can't find it however
numax = k3+k4.*act; %This model was also accurately checked against Ross Miller's 2015 paper so at least the same as his implementation

F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

g = mmass.*barg .* nu .*(1-exp(-0.25-18.2./(nu.*numax)))./(1-exp(-0.25-18.2./numax)); %This equation is a bit ambiguous in the paper
h = mmass.*(barg+barh).*act.*(F_iso-barg./(barg+barh));
s = (bara.*act.*F_iso.*abs(v_ce).*l_ceopt).*(v_ce<=0);

w =  -F_ce.*v_ce.*l_ceopt;

Edot = g+h+s+w; %[W]
% Edot = max(0, Edot);
