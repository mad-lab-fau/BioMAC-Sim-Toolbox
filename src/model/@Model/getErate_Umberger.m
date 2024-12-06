%======================================================================
%> @file @Model/getErate_Umberger.m
%> @brief Model function to calculate the energy rate of a single time
%> step with Umberger et al.'s model
%> @details
%> Details: Model::getErate_Umberger()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step with Umberger et al.'s model. Ross Miller's code was also used as a reference.
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> Umberger's model. It uses
%> the 2010 version (no negative work): Umberger, J R Soc Interface, 2010
%> 2003 (original) Version: Umberger et al., "A model of human muscle 
%> energy expenditure," Computer methods in biomechanics
%> and biomedical engineering, vol. 6, no. 2, pp. 99–111, 2003.
%>
%>
%>Inputs are muscle states of the specific muscle that is considered
%> @param  obj    Gait2dc class object
%> @param  F_ce   Double vector: Muscle force of CE returned by Model.getMuscleCEforces(x) (Model.nMus x 1)
%> @param  stim   Double vector: Stimulation of muscles = neural excitation u (Model.nMus x 1)
%> @param  act    Double vector: Activation of muscles a (Model.nMus x 1)
%> @param  l_ce   Double vector: Length of contractile element of muscles s (Model.nMus x 1)
%> @param  v_ce   Double vector: Velocity of contractile element of muscles (Model.nMus x 1)
%>
%> The output of getErate_Umberger.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot   Double vector: Energy expenditure during one time step (energy rate) in W (Model.nMus x 1)
%> @retval w_ce   Double vector: Mechanical work rate of contractile element 
%>                in W/kg (normalized to muscle mass) (Model.nMus x 1)
%> @retval h_sl   Double vector: Heat rate due to shortening and lengthening of muscles 
%>                in W/kg (normalized to muscle mass) (Model.nMus x 1)
%> @retval h_am   Double vector: Heat rate from the activation of muscles and its maintenance 
%>                in W/kg (normalized to muscle mass) (Model.nMus x 1)
%======================================================================
function [Edot, w_ce, h_sl, h_am] = getErate_Umberger(obj, F_ce, stim, act, l_ce, v_ce)

% Define variables
rho   = 1059.7;   % Muscle density
sigma = 250e3;    % Max muscle stress
S     = 1.5;      % For aerobic activity, S = 1.0 for anaerobic activity

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
FT       = obj.muscles.FT;      % Percentage of fast twitch fibers
if ismember('width', obj.muscles.Properties.VariableNames)
    width    = obj.muscles.width;   % Width of hill curve
else
    width    = sqrt(obj.muscles.kactive);  % Width of hill curve
end
fmax     = obj.muscles.fmax;    % Max muscle force

% Calculate other variables
A = stim .* (stim > act) + 1/2 .* (stim + act).* (stim <= act);
F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship
musMass = (fmax/sigma)*rho.*l_ceopt;

% Calculate velocity of contractile element
if ismember('v_max', obj.muscles.Properties.VariableNames)
    v_cemaxft = obj.muscles.v_max;  % Maximum shortening velocity
else
    v_cemaxft = obj.muscles.vmax;  % Maximum shortening velocity
end
v_cemaxst = v_cemaxft/2.5; 

% Nominal value
Nh_am = 25.*(act<=(1-FT))+(128*FT+25).*(act>(1-FT));
A_am = A.^0.6; % Activation and maintenance
h_am = (Nh_am.*A_am*S).*(l_ce <= 1) + ((0.4*Nh_am + 0.6*Nh_am.*F_iso).*A_am*S).*(l_ce > 1);

% Shortening-Lengthening energy
alpha_ST = 100./v_cemaxst;
alpha_FT = 153./v_cemaxft.*(act > 1-FT); % Zero when activation is less than %ST fibers
alpha_L = 0.3*alpha_ST; %4*alpha_ST; %Anne advised to use 0.3 like it was done in a version of 2010

% Nominal value
Nh_sl_ifSmaller = (alpha_ST.*v_cemaxst < -alpha_ST.*v_ce.*(1-FT));
Nh_sl_ifGreater = (alpha_ST.*v_cemaxst > -alpha_ST.*v_ce.*(1-FT));
Nh_sl = alpha_L .* v_ce .* (v_ce > 0) ... % if v_ce > 0
      + (100 * Nh_sl_ifSmaller - alpha_ST .* v_ce .* (1-FT) .* Nh_sl_ifGreater - alpha_FT .* v_ce .* FT) .*(v_ce <= 0); % if v_ce >= 0  
A_sl = A.*(v_ce > 0)+ A.^(2.0).*(v_ce <=0);
h_sl = Nh_sl.*A_sl*S;
h_sl = h_sl.*F_iso.*(l_ce > 1) + h_sl.*(l_ce <= 1);

% Contractile element work
% => Anne: "The original model subtracts negative work, but this is physically impossible, 
%    so there is an update in 2010 where one coefficient has changed and the negative work 
%    is not subtracted. I have been using both versions, but I would advise changing to the 
%    2010 version."
w_ce = max(0, -F_ce.*v_ce.*l_ceopt./musMass); 

% Energy rate
Edot = (h_am+h_sl+w_ce).*musMass; % Multiply by muscle mass to get W


end
