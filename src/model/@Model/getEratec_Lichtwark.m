%======================================================================
%> @file @Model/getEratec_Lichtwark.m
%> @brief Model function to calculate the energy rate of a single time
%> step with the continuous version of Lichtwark and Wilson's model
%> @details
%> Details: Model::getEratec_Lichtwark()
%>
%> @author Anne Koelewijn
%> @date July, 2021
%======================================================================
%======================================================================
%> @brief Model function to calculate the energy rate of a single time
%> step using a continuous version of Lichtwark & Wilson's model
%>
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%> a continuous version of Lichtwark & Wilson's model, described in the 
%> supplement of Lichtwark & Wilson, J Biomech, 2007. 
%> 
%>Inputs are muscle parameters of the specific muscle that is considered.
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
%> The output of getEratec_Lichtwark.m is a double which is equal to the energy rate
%> of the time step in W/kg
%> @retval Edot        The energy expenditure during the current time step
%> @retval dEdot       The gradien
%======================================================================

function [Edot,dEdot] = getEratec_Lichtwark(obj, F_ce, act, l_ce, v_ce, dFdx,dFdxdot, epsilon)

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
    Ahill = 0.25;%hard coded in 3d model
end
if ismember('gmax', obj.muscles.Properties.VariableNames)
    gmax = obj.muscles.gmax;
else
    gmax = obj.muscles.flen;
end

v_max = obj.muscles.vmax;

v_ce = -v_ce; %shortening positive

F_iso = exp(-(l_ce-1).^2./width.^2); % Force length relationship

% Force velocity relationship
c = v_max.*Ahill.*(gmax-1)./(Ahill+1);
%differentiable
v_ce_g = -v_ce; % should use shortening negative
g = ((v_max+v_ce_g)./(v_max-v_ce_g./Ahill)).*(v_ce_g < 0)+((gmax.*v_ce_g+c)./(v_ce_g+c)).*(v_ce_g >=0);

v_ce_s = 1/2*(v_ce+sqrt(v_ce.^2+epsilon^2)); %shortening velocity, small when lengthening
v_ce_l = 1/2*(v_ce-sqrt(v_ce.^2+epsilon^2)); %lengthening velocity, small when shortening

gamma = 1.5; % number described in supplement
G = 4; % number described in supplement

Hms = gamma*v_max/G^2;
Hml = 0.3*gamma*v_max/G^2+0.7*gamma*v_max/G^2.*exp(-7*v_max.*(g-1));
%sinusoidal approach
x = 1/epsilon;
Hm = (Hms).*(v_ce>0)+Hml.*(v_ce < -pi/x)+(1/2*(Hms+Hml)-(Hml-Hms)/2.*cos(x*v_ce)).*(v_ce >= -pi/x).*(v_ce<=0);
%Hm = 1/2*(Hms+Hml)-(Hml-Hms)/pi.*atan(100*v_ce);

Hs = (v_ce_s/G)+(-0.5*g.*v_ce_l);

%polo units -> isometric force and optimal fiber length
Hdot = 0.3*act.*Hm+0.7*act.*F_iso.*Hm+act.*F_iso.*Hs;
Hdot = Hdot.*l_ceopt.*fmax;

% Contractile element work
w_cebar = F_ce.*v_ce.*l_ceopt;
w = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));

Edot = Hdot+w; %Not sure about the units

if nargout > 1

    %derivatives
    dg_dvce = (1./(v_max-v_ce_g./Ahill) + (v_max+v_ce_g)./(v_max-v_ce_g./Ahill).^2./Ahill).*(v_ce_g < 0) +...
        (gmax./(v_ce_g+c)-(gmax.*v_ce_g+c)./(v_ce_g+c).^2).*(v_ce_g >= 0);

    dvcel_dvce = -1/2.*(1-v_ce./sqrt(v_ce.^2+epsilon^2));
    dvces_dvce = -1/2.*(1+v_ce./sqrt(v_ce.^2+epsilon^2));

    dHml_dvce = -0.7*gamma*v_max/G^2.*(-7*v_max).*exp(-7*v_max.*(g-1)).*dg_dvce;

    % Hm = (Hms).*(v_ce>0)+Hml.*(v_ce < -pi/x)+(1/2*(Hms+Hml)-(Hml-Hms)/2.*cos(x*v_ce)).*(v_ce >= -pi/x).*(v_ce<=0);
    dHm_dvce = -dHml_dvce.*(v_ce < -pi/x)-(1/2*dHml_dvce-dHml_dvce/2.*cos(x*v_ce)+x*(Hml-Hms)/2.*sin(x*v_ce)).*(v_ce>= -pi/x).*(v_ce<=0);

    %Hs = (vbar_ce_s/G)+(-0.5*g.*v_ce);
    dHs_dvce = 1/G*dvces_dvce-0.5*g.*dvcel_dvce-0.5*v_ce_l.*dg_dvce;

    % Hdot = 0.3*a.*Hm+0.7*a.*F_iso.*Hm+a.*F_iso.*Hs;
    dFiso_dxi = zeros(obj.nMus,obj.nDofs*2);
    dFiso_dlce = F_iso.*(-2*(l_ce-1)./width.^2);
    dFiso_dact = zeros(obj.nMus,1);

    dFce_dxi = dFdx([obj.extractState('q'); obj.extractState('qdot')] ,:)';
    dFce_dlce = diag(dFdx(obj.extractState('s'),:));
    dFce_dact = diag(dFdx(obj.extractState('a'),:));
    dFce_dvce = diag(dFdxdot(obj.extractState('s'),:));

    % Hdot = 0.3*a.*Hm+0.7*a.*F_iso.*Hm+a.*F_iso.*Hs;
    dHdot_dxi = bsxfun(@times,0.7*act.*Hm,dFiso_dxi)+bsxfun(@times,act.*Hs,dFiso_dxi);
    dHdot_dact = 0.7*act.*dFiso_dact.*Hm+0.3*Hm+0.7*F_iso.*Hm+act.*dFiso_dact.*Hs+F_iso.*Hs;
    dHdot_dlce = 0.7*act.*dFiso_dlce.*Hm+act.*dFiso_dlce.*Hs;
    dHdot_dvce = 0.3*act.*dHm_dvce+0.7*act.*F_iso.*dHm_dvce+act.*F_iso.*dHs_dvce;

    dHdot_dxi = bsxfun(@times,l_ceopt.*fmax,dHdot_dxi);
    dHdot_dact = dHdot_dact.*l_ceopt.*fmax;
    dHdot_dlce = dHdot_dlce.*l_ceopt.*fmax;
    dHdot_dvce = dHdot_dvce.*l_ceopt.*fmax;

    % w_cebar = F_ce.*v_ce;
    dwbar_dxi = -bsxfun(@times,dFce_dxi,v_ce.*l_ceopt);
    dwbar_dact = v_ce.*l_ceopt.*dFce_dact;
    dwbar_dlce = v_ce.*l_ceopt.*dFce_dlce;
    dwbar_dvce = -F_ce.*l_ceopt + v_ce.*l_ceopt.*dFce_dvce;

    %w_ce = 1/2*(w_cebar+sqrt(w_cebar.^2+epsilon^2));
    dw_dwbar = 1/2*(1+w_cebar./sqrt(w_cebar.^2+epsilon^2));
    dw_dxi = bsxfun(@times,dw_dwbar,dwbar_dxi);
    dw_dact = dw_dwbar.*dwbar_dact;
    dw_dlce = dw_dwbar.*dwbar_dlce;
    dw_dvce = dw_dwbar.*dwbar_dvce;

    dEdot_dxi = dw_dxi+dHdot_dxi;
    dEdot_dact = dw_dact+dHdot_dact;
    dEdot_dlce = dw_dlce+dHdot_dlce;
    dEdot_dvce = dw_dvce+dHdot_dvce;
    dEdot_dstim = zeros(obj.nMus,1);
    dEdot = [dEdot_dxi dEdot_dlce dEdot_dact dEdot_dvce dEdot_dstim];
end