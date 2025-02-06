%======================================================================
%> @file comFixConstraint.m
%> @brief Collocation function to fix center of mass motion
%> @details
%> Details: Collocation::comFixConstraint()
%>
%> @author Anne Koelewijn
%> @date April, 2020
%======================================================================

%======================================================================
%> @brief Function to compute constraint violation of center of mass motion
%>
%> @details
%> This function was never actually used. If you would like to, make sure
%> that it is tested before.
%> 
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least speed of the movement
%> @param targetMotion  Double: target motion of the center of mass
%======================================================================
function output = comFixConstraint(obj,option,X, targetMotion)

% compute demanded output
if strcmp(option,'confun')
    CoM = obj.getComMotion(X);
    %calculate vertical motion
    output = CoM - targetMotion;
elseif strcmp(option,'jacobian')
    output = getComMotionDeriv(obj,X);
else
    error('Unknown option.');
end
end

