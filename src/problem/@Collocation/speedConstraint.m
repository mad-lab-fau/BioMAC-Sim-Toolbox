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
function output = speedConstraint(obj,option,X,targetSpeed)
%% check input parameter
if  ~isfield(obj.idx,'speed') % check whether controls are stored in X
    error('State vector X does not contain speed.')
end
    
%% compute demanded output
% forward translation in x direction
speed = X(obj.idx.speed);
% if  isa(targetSpeed, 'TrackingData')
%     idxVars = find(ismember(targetSpeed.variables.type, 'speed')); %indices of tracking variables
%     targetSpeed = targetSpeed.variables.mean{idxVars};
% end
if strcmp(option,'confun') %constraints of periodicity constraint
    output =  speed - targetSpeed;   
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(1,length(X),obj.Jnnz);
    output(1, obj.idx.speed) = 1;    
else
    error('Unknown option');
end
end
