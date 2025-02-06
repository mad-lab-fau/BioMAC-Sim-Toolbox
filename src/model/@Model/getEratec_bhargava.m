%======================================================================
%> @file @Model/getEratec_bhargava.m
%> @brief Model function to calculate the energy rate of a single time
%> step with the continuous version of Bhargava et al.'s model
%> @details
%> Details: Model::getEratec_bhargava()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using a continuous version of Bhargava et al.'s model
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> a continuous version of Bhargava et al.'s model. Line 60, phi, was
%> chosen based on graphs of phi as a function of t_stim, since t_stim 
%> cannot be differentiated. Model paper: bhargava et al., J Biomech, 2004
%> 
%>Inputs are muscle states of the specific muscle that is considered.
%> @param  obj         Model object
%> @param F_ce         Force in the contractile element
%> @param stim         Stimulation of the muscle (between 0 and 1)
%> @param act          Activation level of the muscle (between 0 and 1)
%> @param l_ce         Normalized length of the contractile element
%> @param v_ce         Normalized velocity of the contractile element
%> @param dFdx         Derivative of CE forces with respect to the state
%> @param dFdxdot      Derivative of CE forces with respect to the state
%>                     derivative
%> @param epsilon      Measure of the nonlinearity of the continuous functions
%>
%> The output of getEratec_bhargava.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot        The energy expenditure during the current time step
%> @retval dEdot       The gradient of Edot with respect to all states and
%>                     inputs

%======================================================================

function [Edot,dEdot] = getEratec_bhargava(obj, F_ce, stim, act, l_ce, v_ce, dFdx, dFdxdot,epsilon)

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

A_f = 133; %W/kg
A_s = 40;  %W/kg
M_f = 111; %W/kg
M_s = 74;  %W/kg

u_f = 1-cos(pi/2*stim);
u_s = sin(pi/2*stim);
phi = 0.2;%0.06+exp(-t_stim.*u/tau_phi);
F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

A = phi.*mmass.*(FT*A_f.*u_f+ST*A_s.*u_s);

%lM = 0.5.*(l_ce<0.5)+l_ce.*(l_ce<1).*(l_ce>=0.5)+(-2*l_ce+3).*(l_ce<1.5).*(l_ce>=1);
lm1 = 0.5+1/2*(l_ce-0.5+sqrt((l_ce-0.5).^2+epsilon^2));
lm2 = 1/2*((-2*l_ce+3)+sqrt((-2*l_ce+3).^2+epsilon^2));
lM = lm2-1/2*(lm2-lm1+sqrt((lm2-lm1).^2+epsilon^2));

M = lM.*mmass.*(FT*M_f.*u_f+ST*M_s.*u_s);

%alpha = (0.16*a.*F_iso.*F_max+0.18*F_ce).*(v_ce <=0)+(0.157*F_ce).*(v_ce>0);
v_ce_s = 1/2*(v_ce-sqrt((v_ce).^2+epsilon^2)).*l_ceopt; %shortening velocity
v_ce_l = 1/2*(v_ce+sqrt((v_ce).^2+epsilon^2)).*l_ceopt; %lengthening velocity

alpha_l = 0.16*act.*F_iso.*fmax+0.18*F_ce; %for shortening
alpha_s = 0.157*F_ce; %for lengthening
S = -alpha_l.*v_ce_l-alpha_s.*v_ce_s;

% Contractile element work
w_cebar = -F_ce.*v_ce.*l_ceopt;
w = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));

Edot = A+M+S+w; %[W]

if nargout > 1
    % Derivatives
    %F_iso = exp(-(l_ce-l_ceopt).^2./(W.*l_ceopt).^2); % Force length relationship
    dFiso_dlce = F_iso.*(-2*(l_ce-1)./width.^2);

    dFce_dxi = dFdx([obj.extractState('q'); obj.extractState('qdot')] ,:)';
    dFce_dlce = diag(dFdx(obj.extractState('s'),:));
    dFce_dact = diag(dFdx(obj.extractState('a'),:));
    dFce_dvce = diag(dFdxdot(obj.extractState('s'),:));
    dFce_dstim = zeros(obj.nMus,1);

    %u_f = 1-cos(pi/2*u);
    duf_dstim = pi/2*sin(pi/2*stim);
    duf_dlce = zeros(obj.nMus,1);

    %u_s = sin(pi/2*u);
    dus_dstim = pi/2*cos(pi/2*stim);
    dus_dlce = zeros(obj.nMus,1);

    % phi = 0.2;%0.06+exp(-t_stim.*a/tau_phi);
    dphi_dact = zeros(obj.nMus,1);%-t_stim/tau_phi.*exp(-t_stim.*a/tau_phi);

    %A = phi.*mmass.*(FT*A_f.*u_f+ST*A_s.*u_s);
    dA_dxi = zeros(obj.nMus,obj.nDofs*2);%phi.*mmass.*(FT*A_f.*duf_dxi+ST*A_s.*dus_dxi);
    dA_dact = dphi_dact.*mmass.*(FT*A_f.*u_f+ST*A_s.*u_s);
    dA_dlce = phi.*mmass.*(FT*A_f.*duf_dlce+ST*A_s.*dus_dlce);
    dA_dstim = phi.*mmass.*(FT*A_f.*duf_dstim+ST*A_s.*dus_dstim);

    dlm1_dlce = 1/2*(1+(l_ce-0.5)./sqrt((l_ce-0.5).^2+epsilon^2));
    dlm2_dlce = -(1+(-2*l_ce+3)./sqrt((-2*l_ce+3).^2+epsilon^2));
    dlM_dlce = dlm2_dlce-1/2*(1+(lm2-lm1)./+sqrt((lm2-lm1).^2+epsilon^2)).*(dlm2_dlce-dlm1_dlce);

    %M = lM.*mmass.*(FT*M_f.*u_f+ST*M_s.*u_s);
    dM_dxi = zeros(obj.nMus,obj.nDofs*2);%lM.*mmass.*(FT*M_f.*duf_dxi+ST*M_s.*dus_dxi);
    dM_dstim = lM.*mmass.*(FT*M_f.*duf_dstim+ST*M_s.*dus_dstim);
    dM_dlce = lM.*mmass.*(FT*M_f.*duf_dlce+ST*M_s.*dus_dlce)+dlM_dlce.*mmass.*(FT*M_f.*u_f+ST*M_s.*u_s);

    %alpha_l = 0.16*a.*F_iso.*F_max+0.18*F_ce; %for lengthening
    dalphal_dxi = 0.18*dFce_dxi;
    dalphal_dact = 0.16*F_iso.*fmax+0.18*dFce_dact;
    dalphal_dlce = 0.16*act.*dFiso_dlce.*fmax+0.18*dFce_dlce;
    dalphal_dvce = 0.18*dFce_dvce;
    dalphal_dstim = 0.18*dFce_dstim;

    %alpha_s = 0.157*F_ce; %for shortening
    dalphas_dxi = 0.157*dFce_dxi;
    dalphas_dact = 0.157*dFce_dact;
    dalphas_dlce = 0.157*dFce_dlce;
    dalphas_dvce = 0.157*dFce_dvce;
    dalphas_dstim = 0.157*dFce_dstim;

    %v_ce_s = 1/2*(v_ce-sqrt((v_ce).^2+epsilon^2)); %shortening velocity
    % v_ce_l = 1/2*(v_ce+sqrt((v_ce).^2+epsilon^2)); %lengthening velocity
    dvce_ldvce = 1/2.*(1+v_ce./sqrt(v_ce.^2+epsilon^2)).*l_ceopt;
    dvce_sdvce = 1/2.*(1-v_ce./sqrt(v_ce.^2+epsilon^2)).*l_ceopt;

    % S = -alpha_l.*v_ce_l-alpha_s*v_ce_s;
    dS_dxi = -bsxfun(@times,dalphal_dxi,v_ce_l)-bsxfun(@times,dalphas_dxi,v_ce_s);
    dS_dact = -dalphal_dact.*v_ce_l -dalphas_dact.*v_ce_s;
    dS_dlce = -dalphal_dlce.*v_ce_l-dalphas_dlce.*v_ce_s;
    dS_dvce = -alpha_l.*dvce_ldvce - alpha_s.*dvce_sdvce-dalphal_dvce.*v_ce_l - dalphas_dvce.*v_ce_s;
    dS_dstim = -dalphal_dstim.*v_ce_l-dalphas_dstim.*v_ce_s;

    % w_cebar = -F_ce.*v_ce;
    dwbar_dxi = -bsxfun(@times,dFce_dxi,v_ce.*l_ceopt);
    dwbar_dact = -v_ce.*l_ceopt.*dFce_dact;
    dwbar_dlce = -v_ce.*l_ceopt.*dFce_dlce;
    dwbar_dvce = -F_ce.*l_ceopt - v_ce.*l_ceopt.*dFce_dvce;
    dwbar_dstim = -v_ce.*l_ceopt.*dFce_dstim;

    %w_ce = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));
    dwce_dwcebar = 1/2*(1+w_cebar./sqrt(w_cebar.^2+epsilon^2));
    dW_dxi = bsxfun(@times,dwce_dwcebar,dwbar_dxi);
    dW_dact = dwce_dwcebar.*dwbar_dact;
    dW_dlce = dwce_dwcebar.*dwbar_dlce;
    dW_dvce = dwce_dwcebar.*dwbar_dvce;
    dW_dstim = dwce_dwcebar.*dwbar_dstim;

    dEdot_dxi = dA_dxi+dM_dxi+dS_dxi+dW_dxi;
    dEdot_dact = dA_dact+dS_dact+dW_dact;
    dEdot_dlce = dA_dlce+dM_dlce+dS_dlce+dW_dlce;
    dEdot_dstim = dA_dstim+ dM_dstim+dS_dstim+dW_dstim;
    dEdot_dvce = dS_dvce+dW_dvce;
    dEdot = [dEdot_dxi dEdot_dlce dEdot_dact dEdot_dvce dEdot_dstim];
end