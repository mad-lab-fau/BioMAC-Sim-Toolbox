%======================================================================
%> @file @Model/getEratec_Houdijk.m
%> @brief Model function to calculate the energy rate of a single time
%> step with the continuous version of Houdijk et al.'s model
%> @details
%> Details: Model::getEratec_Houdijk()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using a continuous version of Houdijk et al.'s model
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> a continuous version of Houdijk et al.'s model. Houdijk et al., J
%> Biomech, 2006.
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
%>Inputs are muscle states of the specific muscle that is considered.
%> @param  obj         Model object
%> @param F_ce         Force in the contractile element
%> @param act          Activation level of the muscle (between 0 and 1)
%> @param l_ce         Normalized length of the contractile element
%> @param v_ce         Normalized velocity of the contractile element
%> @param dFdx         Derivative of CE forces with respect to the state
%> @param dFdxdot      Derivative of CE forces with respect to the state
%>                     derivative
%> @param epsilon      Measure of the nonlinearity of the continuous functions
%>
%> The output of getEratec_Houdijk.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot        The energy expenditure during the current time step
%> @retval dEdot       The gradient of Edot with respect to all states and
%>                     inputs%
%======================================================================


function [Edot,dEdot] = getEratec_Houdijk(obj, F_ce, act, l_ce, v_ce, dFdx, dFdxdot, epsilon)

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
bargs = 0.45*24.4;
bargf = 0.35*150;
barhs = 24.4-bargs;
barhf = 150-bargf;

barg = bargs*ST+bargf*FT;
barh = barhs*ST+barhf*FT;
k3 = 6*ST+12*FT;
k4 = 8*ST+14*FT;
bara = 0.16*fmax.*ST+0.28*fmax.*FT;

nu = act.^2;
numax = k3+k4.*act;

F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

g = mmass.*barg .* nu .*(1-exp(-0.25-18.2./(nu.*numax)))./(1-exp(-0.25-18.2./numax));
h = mmass.*(barg+barh).*act.*(F_iso-barg./(barg+barh));
v_ce_s = 1/2*(v_ce-sqrt((-v_ce).^2+epsilon^2)).*l_ceopt; %shortening velocity
s = bara.*act.*F_iso.*-v_ce_s;

% Contractile element work
w_cebar = -F_ce.*v_ce.*l_ceopt; %max(0,-F_ce.*v_ce./mmass);
w = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));

Edot = g+h+s+w; %[W]

if nargout > 1

    %Derivatives
    dFiso_dxi = zeros(obj.nMus,obj.nDofs*2);
    dFiso_dlce = F_iso.*(-2*(l_ce-1)./width.^2); %-F_iso.*(2/W.^2*(l_ce./l_ceopt-1)); %-2.0*x*F1 / m->Width;
    dFiso_dact = zeros(obj.nMus,1);

    dFce_dxi = dFdx([obj.extractState('q'); obj.extractState('qdot')] ,:)';
    dFce_dlce = diag(dFdx(obj.extractState('s'),:));
    dFce_dact = diag(dFdx(obj.extractState('a'),:));
    dFce_dvce = diag(dFdxdot(obj.extractState('s'),:));

    dnu_dact = 2*act;
    dnumax_dact = k4;

    % g = mmass.*barg .* nu .*(1-exp(-0.25-18.2./(nu.*numax)))./(1-exp(-0.25-18.2./numax));
    dg_dxi = zeros(obj.nMus,obj.nDofs*2);
    dg_dlce = zeros(obj.nMus,1);
    dg_dact = mmass.*barg.*dnu_dact.*(1-exp(-0.25-18.2./(nu.*numax)))./(1-exp(-0.25-18.2./numax)) ...
        -mmass.*barg.*(18.2*exp(-0.25-18.2./(nu.*numax))./(nu.*numax)./(1-exp(-0.25-18.2./numax))).*dnu_dact ...
        -mmass.*barg.*nu.*(18.2*exp(-0.25-18.2./(nu.*numax))./(nu.*numax.^2)./(1-exp(-0.25-18.2./numax))).*dnumax_dact ...
        +mmass.*barg.*nu*18.2.*(1-exp(-0.25-18.2./(nu.*numax))).*exp(-0.25-18.2./numax)./(numax.^2.*(1-exp(-0.25-18.2./numax)).^2).*dnumax_dact;

    % h = mmass.*(barg+barh).*a.*(F_iso-barg./(barg+barh));
    dh_dxi = bsxfun(@times,mmass.*(barg+barh).*act,dFiso_dxi);
    dh_dlce = mmass.*(barg+barh).*act.*dFiso_dlce;
    dh_dact = mmass.*(barg+barh).*(F_iso-barg./(barg+barh))+mmass.*(barg+barh).*act.*dFiso_dact;

    dvces_dvce = 1/2.*(1-v_ce./sqrt(v_ce.^2+epsilon^2)).*l_ceopt;

    %s = bara.*a.*F_iso.*-v_ce_s;
    ds_dxi = zeros(obj.nMus,obj.nDofs*2);
    ds_dact = bara.*F_iso.*-v_ce_s+bara.*act.*dFiso_dact.*-v_ce_s;
    ds_dlce = bara.*act.*dFiso_dlce.*-v_ce_s;
    ds_dvce = bara.*act.*F_iso.*-dvces_dvce;

    dwbar_dxi = -bsxfun(@times,dFce_dxi,v_ce.*l_ceopt);
    dwbar_dact = -v_ce.*l_ceopt.*dFce_dact;
    dwbar_dlce = -v_ce.*l_ceopt.*dFce_dlce;
    dwbar_dvce = -F_ce.*l_ceopt - v_ce.*l_ceopt.*dFce_dvce;

    %w_ce = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));
    dwce_dwcebar = 1/2*(1+w_cebar./sqrt(w_cebar.^2+epsilon^2));
    dw_dxi = bsxfun(@times,dwce_dwcebar,dwbar_dxi);
    dw_dact = dwce_dwcebar.*dwbar_dact;
    dw_dlce = dwce_dwcebar.*dwbar_dlce;
    dw_dvce = dwce_dwcebar.*dwbar_dvce;

    dEdot_dxi = ds_dxi+dw_dxi+dg_dxi+dh_dxi;
    dEdot_dact = ds_dact+dw_dact+dg_dact+dh_dact;
    dEdot_dlce = ds_dlce+dw_dlce+dg_dlce+dh_dlce;
    dEdot_dvce = ds_dvce+dw_dvce;
    dEdot_dstim = zeros(obj.nMus,1);
    dEdot = [dEdot_dxi dEdot_dlce dEdot_dact dEdot_dvce dEdot_dstim];
end