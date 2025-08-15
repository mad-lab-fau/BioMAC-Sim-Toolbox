%======================================================================
%> @file speedConstraint.m
%> @brief Collocation function to compute speed constraint
%> @details
%> Details: Collocation::speedConstraint()
%>
%> @author Marlies Nitschke
%> @date January, 2018
%======================================================================

%======================================================================
%> @brief Computes constraint violation demanding periodic movement
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least speed of the movement
%> @param targetSpeed   Double: target speed to reach
%======================================================================
function output = speedConstraint_IMU(obj,option,X,targetSpeed)
%% check input parameter
if  ~isfield(obj.idx,'speed') % check whether controls are stored in X
    error('State vector X does not contain speed.')
end
    
%% compute demanded output
% forward translation in x direction
speed = X(obj.idx.speed);

if strcmp(option,'confun') %constraints of periodicity constraint
    % output =    % TODO1. Add code here. We want to constrain the difference between speed and targetSpeed
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(1,length(X),obj.Jnnz);
    %output(1, ) = ; % TODO 2. Add the constraint's derivative at the correct index. Make sure to use the derivative of the function in line 33    
else
    error('Unknown option');
end
end

