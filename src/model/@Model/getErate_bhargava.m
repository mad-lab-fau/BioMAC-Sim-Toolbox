%======================================================================
%> @file @Model/getErate_bhargava.m
%> @brief Model function to calculate the energy rate of a single time
%> step with Bhargava et al.'s model
%> @details
%> Details: Model::getErate_bhargava()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using Bhargava et al.'s model
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> a continuous version of Bhargava et al.'s model. Model paper: 
%> Bhargava et al., J Biomech, 2004
%> 
%>Inputs are muscle states of the specific muscle that is considered.
%> @param  obj         Model object
%> @param F_ce         Force in the contractile element
%> @param stim         Stimulation of the muscle (between 0 and 1)
%> @param act          Activation level of the muscle (between 0 and 1)
%> @param l_ce         Normalized length of the contractile element
%> @param v_ce         Normalized velocity of the contractile element
%> @param t_stim       Duration of stimulation of muscle
%>
%> The output of getEratec_bhargava.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot        The energy expenditure during the current time step
%======================================================================

function Edot = getErate_bhargava(obj, F_ce, stim, act, l_ce, v_ce,t_stim)

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

mmass = (fmax/sigma).*l_ceopt*rho;

ST = 1-FT;

F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

A_f = 133; %W/kg
A_s = 40;  %W/kg
M_f = 111; %W/kg
M_s = 74;  %W/kg

tau_phi = 45/1000; %s

u_f = 1-cos(pi/2*stim);
u_s = sin(pi/2*stim);
phi = 0.06+exp(-t_stim.*stim/tau_phi);

A = phi.*mmass.*(FT*A_f.*u_f+ST*A_s.*u_s);

lM = 0.5.*(l_ce<0.5)+l_ce.*(l_ce<1).*(l_ce>0.5)+(-2*l_ce+3).*(l_ce<1.5).*(l_ce>1);

M = lM.*mmass.*(FT*M_f.*u_f+ST*M_s.*u_s);

alpha = (0.16*act.*F_iso.*fmax+0.18*F_ce).*(v_ce <=0)+(0.157*F_ce).*(v_ce>0);
S = -alpha.*v_ce.*l_ceopt;

w =  F_ce.*-v_ce.*l_ceopt;%max(0,F_ce.*-v_ce.*l_ceopt);%

Edot = A+M+S+w; %[W]
% Edot = max(0, Edot);