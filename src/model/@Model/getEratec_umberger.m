%======================================================================
%> @file @Model/getEratec_umberger.m
%> @brief Model function to calculate the energy rate of a single time
%> step with the continuous version of Umberger et al.'s model
%> @details
%> Details: Model::getEratec_umberger()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step with the continuous version of Umberger et al.'s model
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> a continuous version of Umberger's model. This one is published in
%> Koelewijn, Dorscky et al., Comp Meth Biomed Biomech Eng, 2016, and uses
%> the 2010 version (no negative work): Umberger, J R Soc Interface, 2010
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
%> The output of getEratec_umberger.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot        The energy expenditure during the current time step
%> @retval dEdot       The gradient of Edot with respect to all states and
%>                     inputs
%======================================================================

function [Edot, dEdot] = getEratec_umberger(obj, F_ce, stim, act, l_ce, v_ce, dFdx, dFdxdot, epsilon)

% Define variables
rho   = 1059.7;   % Muscle density
sigma = 250e3;    % Max muscle stress
S = 1.5;          % For aerobic movement. Use 1.0 for anaerobic

% Get muscle variables
l_ceopt  = obj.muscles.lceopt;  % Optimal fiber length
FT       = obj.muscles.FT;      % Percentage of fast twitch fibers
if ismember('width', obj.muscles.Properties.VariableNames)
    width    = obj.muscles.width;   % Width of hill curve
else
    width    = sqrt(obj.muscles.kactive);  % Width of hill curve
end
fmax     = obj.muscles.fmax;    % Max muscle force

% Speed improvements by not loading these
nDofs = obj.nDofs;
nMus = obj.nMus;

mmass = (fmax/sigma).*l_ceopt*rho;

% Calculate other variables
ACST = 1/2*(act+stim);
A = (stim+1/2*((ACST-stim)+sqrt((ACST-stim).^2+epsilon^2))); %STIM.*(STIM>ACT)+1/2.*(STIM+ACT).*(STIM <= ACT);

F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

% Calculate velocity of contractile element
v_cemaxft = obj.muscles.vmax;  % Maximum shortening velocity
v_cemaxst = v_cemaxft/2.5; 

v_ce_l = 1/2*(v_ce+sqrt(v_ce.^2+epsilon^2)); %lengthening velocity, small when shortening
v_ce_s = 1/2*(v_ce-sqrt((-v_ce).^2+epsilon^2)); %shortening velocity, small when lengthening

% Activation and maintenance
A_am = real(A.^0.6);

% Ross Fiso-thing
F_iso1 = (l_ce <= 1)+F_iso.*(l_ce > 1);

% Nominal value
Nh_am = 25.*(act<=(1-FT))+(128*FT+25).*(act>(1-FT));
h_am = (0.4*Nh_am+0.6*Nh_am.*F_iso1).*A_am*S;

% Shortening-Lengthening energy
alpha_ST = 100./v_cemaxst;
alpha_FT = 153./v_cemaxft.*(act > 1-FT); %Zero when activation is less than %ST fibers
alpha_L = 0.3*alpha_ST;

% Nominal value
Nh_sl = alpha_L.*v_ce_l + (100*(alpha_ST.*v_cemaxst<-alpha_ST.*v_ce_s.*(1-FT))-alpha_ST.*v_ce_s.*(1-FT).*(alpha_ST.*v_cemaxst > -alpha_ST.*v_ce_s.*(1-FT))-alpha_FT.*v_ce_s.*FT);

A2 = A.^(2.0);%A2 = A.*(v_ce > 0)+ A.^(2.0).*(v_ce <=0);%
% vbarfunc = tanh(v_ce);%v_ce./sqrt(1+v_ce.^2);
A_sl = A2;%1/2*(A-A2).*vbarfunc+1/2*(A2+A);%A.;%A.*(v_ce > 0)+ A.^(2.0).*(v_ce <=0);

h_sl = Nh_sl.*A_sl*S;
h_sl = h_sl.*F_iso1;

% Contractile element work
w_cebar = -F_ce.*v_ce.*l_ceopt./mmass;
w_ce = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));

Edot = (h_am+h_sl+w_ce).*mmass;
% Edot = (1+1/2*((Edot1-1)+sqrt((Edot1-1).^2+epsilon^2))).*mmass;% (1*(1> h_am+h_sl+w_ce)+ (h_am+h_sl+w_ce).*(1 <= h_am+h_sl+w_ce)).*mmass; %Multiply by muscle mass to get W

if nargout > 1
    % Derivatives
    % exp(-(l_ce-l_ceopt).^2./(W.*l_ceopt).^2);
    dFiso_dxi = zeros(nMus,nDofs*2);
    dFiso_dlce = F_iso.*(-2*(l_ce-1)./width.^2); %-F_iso.*(2/W.^2*(l_ce./l_ceopt-1)); %-2.0*x*F1 / m->Width;
    dFiso_dact = zeros(nMus,1);
    dFiso_dstim = zeros(nMus,1);

    % ACST = 1/2*(ACT+STIM);
    %A = (STIM+1/2*((ACST-STIM)+sqrt((ACST-STIM).^2+epsilon^2)));
    dA_dxi = zeros(nMus,nDofs*2);
    dA_dlce = zeros(nMus,1);
    dA_dact = 1/4*(1+(ACST-stim)./sqrt((ACST-stim).^2+epsilon^2));%ones(nMus,1);%
    dA_dstim = ones(nMus,1)-1/4*(1+(ACST-stim)./sqrt((ACST-stim).^2+epsilon^2));zeros(nMus,1);%

    dAam_dxi = zeros(nMus,nDofs*2);
    dAam_dlce = zeros(nMus,1);
    dAam_dact = real(0.6*A.^(-0.4).*dA_dact);
    dAam_dstim = real(0.6*A.^(-0.4).*dA_dstim);

    % Ross Fiso-thing
    % F_iso1 = (l_ce <= 1)+F_iso.*(l_ce > 1);
    dFiso1_dlce = dFiso_dlce.*(l_ce > 1);
    dFiso1_dact = dFiso_dact.*(l_ce > 1);
    dFiso1_dstim = dFiso_dstim.*(l_ce > 1);
    dFiso1_dxi = zeros(nMus,nDofs*2);

    dham_dxi = zeros(nMus,nDofs*2);
    for i = 1:nDofs*2
        dFiso1_dxi(:,i) = dFiso1_dxi(:,i)+(dFiso_dxi(:,i)).*(l_ce > l_ceopt);
        dham_dxi(:,i) = 0.6*S*dFiso1_dxi(:,i).*(A_am.*Nh_am)+S*dAam_dxi(:,i).*(0.4*Nh_am+0.6*F_iso1.*Nh_am);
    end
    dham_dlce = A_am*S*0.6.*Nh_am.*dFiso1_dlce+(0.4*Nh_am+0.6*F_iso1.*Nh_am)*S.*dAam_dlce;
    dham_dact = A_am*S*0.6.*Nh_am.*dFiso1_dact+(0.4*Nh_am+0.6*F_iso1.*Nh_am)*S.*dAam_dact;
    dham_dstim = A_am*S*0.6.*Nh_am.*dFiso1_dstim+(0.4*Nh_am+0.6*F_iso1.*Nh_am)*S.*dAam_dstim;

    dFce_dxi = dFdx([obj.extractState('q'); obj.extractState('qdot')] ,:)';
    dFce_dlce = diag(dFdx(obj.extractState('s'),:));
    dFce_dact = diag(dFdx(obj.extractState('a'),:));
    dFce_dvce = diag(dFdxdot(obj.extractState('s'),:));
    dFce_dstim = zeros(nMus,1);

    % Wrt lengthening and shortening velocity
    % v_ce_l = 1/2*(v_ce+sqrt(v_ce.^2+epsilon^2)); %lengthening velocity, small when shortening
    % v_ce_s = 1/2*(v_ce-sqrt((-v_ce).^2+epsilon^2)); %shortening velocity, small when lengthening
    dv_ce_ldvce = 1/2.*(1+v_ce./sqrt(v_ce.^2+epsilon^2));
    dv_ce_sdvce = 1/2.*(1-v_ce./sqrt((-v_ce).^2+epsilon^2));

    % Nominal value
    % Nh_sl = alpha_L.*v_ce_l + (100*(alpha_ST.*v_cemaxst<-alpha_ST.*v_ce_s.*(1-FT))-alpha_ST.*v_ce_s.*(1-FT).*(alpha_ST.*v_cemaxst > -alpha_ST.*v_ce_s.*(1-FT))-alpha_FT.*v_ce_s.*FT);
    Nh_sl = alpha_L.*v_ce_l + (100*(alpha_ST.*v_cemaxst<-alpha_ST.*v_ce_s.*(1-FT))-alpha_ST.*v_ce_s.*(1-FT).*(alpha_ST.*v_cemaxst > -alpha_ST.*v_ce_s.*(1-FT))-alpha_FT.*v_ce_s.*FT);
    dNhsl_dvce = alpha_L.*dv_ce_ldvce - alpha_ST.*dv_ce_sdvce.*(1-FT).*(alpha_ST.*v_cemaxst > -alpha_ST.*v_ce_s.*(1-FT)) - alpha_FT.*FT.*dv_ce_sdvce;

    %A2 = A.^(2.0);
    dA2_dxi = 2*bsxfun(@times,dA_dxi,A);
    dA2_dact = 2*dA_dact.*A;
    dA2_dlce = 2*dA_dlce.*A;
    dA2_dstim = 2*dA_dstim.*A;

    %A_sl = A2;
    dAsl_dxi = dA2_dxi;
    dAsl_dact = dA2_dact;
    dAsl_dlce = dA2_dlce;
    dAsl_dstim = dA2_dstim;

    % h_sl = Nh_sl.*A_sl*S;
    % h_sl = h_sl.*F_iso1;
    dhsl_dxi = zeros(nMus,nDofs*2);
    for i=1:(nDofs*2)
        dhsl_dxi(:,i) = dAsl_dxi(:,i).*Nh_sl.*F_iso1+dFiso1_dxi(:,i).*Nh_sl.*A_sl*S;
    end
    dhsl_dact = Nh_sl.*dAsl_dact*S.*F_iso1+Nh_sl.*A_sl*S.*dFiso1_dact;
    dhsl_dlce = Nh_sl.*dAsl_dlce*S.*F_iso1+Nh_sl.*A_sl*S.*dFiso1_dlce;
    dhsl_dvce = A_sl*S.*F_iso1.*dNhsl_dvce;
    dhsl_dstim = Nh_sl.*dAsl_dstim*S.*F_iso1+Nh_sl.*A_sl*S.*dFiso1_dstim;

    % w_cebar = -F_ce.*v_ce.*l_ceopt./mmass
    dwcebar_dxi = -bsxfun(@times,dFce_dxi,v_ce.*l_ceopt./mmass);
    dwcebar_dact = -v_ce.*l_ceopt./mmass.*dFce_dact;
    dwcebar_dlce = -v_ce.*l_ceopt./mmass.*dFce_dlce;
    dwcebar_dvce = -F_ce.*l_ceopt./mmass-v_ce.*l_ceopt./mmass.*dFce_dvce;
    dwcebar_dstim = -v_ce.*l_ceopt./mmass.*dFce_dstim;

    %w_ce = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));
    dwce_dwcebar = 1/2*(1+w_cebar./sqrt(w_cebar.^2+epsilon^2));
    dwce_dxi = bsxfun(@times,dwce_dwcebar,dwcebar_dxi);
    dwce_dact = dwce_dwcebar.*dwcebar_dact;
    dwce_dlce = dwce_dwcebar.*dwcebar_dlce;
    dwce_dvce = dwce_dwcebar.*dwcebar_dvce;
    dwce_dstim = dwce_dwcebar.*dwcebar_dstim;

    dEdot_dxi = zeros(nMus,nDofs*2);
    for i = 1:nDofs*2
        dEdot_dxi(:,i) = ((dwce_dxi(:,i)+dhsl_dxi(:,i)+dham_dxi(:,i)).*mmass);
    end
    dEdot_dact = (mmass.*(dwce_dact+dhsl_dact+dham_dact));
    dEdot_dlce = (mmass.*(dwce_dlce+dhsl_dlce+dham_dlce));
    dEdot_dvce = (mmass.*(dwce_dvce+dhsl_dvce));
    dEdot_dstim = (mmass.*(dwce_dstim+dhsl_dstim+dham_dstim));
    dEdot = [dEdot_dxi dEdot_dlce dEdot_dact dEdot_dvce dEdot_dstim];

end
end