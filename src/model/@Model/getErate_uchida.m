%======================================================================
%> @file @Model/getErate_uchida.m
%> @brief Model function to calculate the energy rate of a single time
%> step with Uchida et al.'s model
%> @details
%> Details: Model::getErate_uchida()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%step for post processing. This is the discontinuous version of the code
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> Uchida's model: Uchida et al., PLOS ONE, 2016
%> 
%>Inputs are muscle states of the specific muscle that is considered
%> @param  obj    Model class object
%> @param  F_ce   Double vector: Muscle force of CE returned by obj.getMuscleCEforces(x) (Model.nMus x 1)
%> @param  stim   Double vector: Stimulation of muscles = neural excitation u (Model.nMus x 1)
%> @param  act    Double vector: Activation of muscles a (Model.nMus x 1)
%> @param  l_ce   Double vector: Length of contractile element of muscles s (Model.nMus x 1)
%> @param  v_ce   Double vector: Velocity of contractile element of muscles (Model.nMus x 1)
%>
%> The output of getErate_uchida.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot   Double vector: Energy expenditure during one time step (energy rate) in W (Model.nMus x 1)
%> @retval w_ce   Double vector: Mechanical work rate of contractile element 
%>                in W/kg (normalized to muscle mass) (Model.nMus x 1)
%> @retval h_sl   Double vector: Heat rate due to shortening and lengthening of muscles 
%>                in W/kg (normalized to muscle mass) (Model.nMus x 1)
%> @retval h_am   Double vector: Heat rate from the activation of muscles and its maintenance 
%>                in W/kg (normalized to muscle mass) (Model.nMus x 1)
%======================================================================

function [Edot,w_ce, h_sl, h_am] = getErate_uchida(obj, F_ce, stim, act, l_ce, v_ce)

% Define variables
rho   = 1059.7;   % Muscle density
sigma = 250e3;    % Max muscle stress
S = 1.5; % for aerobic activity

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
FT       = obj.muscles.FT;      % Percentage of fast twitch fibers
width    = obj.muscles.width;   % Width of hill curve
fmax     = obj.muscles.fmax;    % Max muscle force

mmass = (fmax/sigma)*rho.*l_ceopt;

%Calculate % of FT and ST fibers used
stim_f = 1-cos(pi/2*stim);
stim_s = sin(pi/2*stim);

P_FT = FT.*stim_f./((1-FT).*stim_s+FT.*stim_f);

% Calculate other variables
A = stim.*(stim>act)+1/2.*(stim+act).*(stim <= act);
F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

% Calculate velocity of contractile element
vbar_cemaxft = obj.muscle.v_max;  % Maximum shortening velocity
vbar_cemaxst = vbar_cemaxft/2.5; 

% Activation and maintenance
A_am = A.^0.6;

% Nominal value
Nh_am = 25.*(act<=(1-P_FT))+(128*P_FT+25).*(act>(1-P_FT));
h_am = (Nh_am.*A_am*S).*(l_ce <= 1) + ((0.4*Nh_am+0.6*Nh_am.*F_iso).*A_am*S).*(l_ce > 1);

% Shortening-Lengthening energy
alpha_ST = 100/vbar_cemaxst;
alpha_FT = 153/vbar_cemaxft*(act > 1-P_FT); %Zero when activation is less than %ST fibers
alpha_L = 4*alpha_ST; 

% Nominal value
Nh_sl = alpha_L.*v_ce.*(v_ce > 0) + (100*(alpha_ST*vbar_cemaxst<-alpha_ST*v_ce.*(1-P_FT))-alpha_ST*v_ce.*(1-P_FT).*(alpha_ST*vbar_cemaxst > -alpha_ST*v_ce.*(1-P_FT))-alpha_FT.*v_ce.*P_FT).*(v_ce <= 0);

A_sl = A.*(v_ce > 0)+ A.^(2.0).*(v_ce <=0);
h_sl = Nh_sl.*A_sl*S;
h_sl = h_sl.*F_iso.*(l_ce > 1) + h_sl.*(l_ce <= 1);

% Contractile element work
w_ce = -F_ce.*v_ce./mmass;

Edot = (h_am+h_sl+w_ce).*mmass; %Multiply by muscle mass to get W
Edot = max(0, Edot);