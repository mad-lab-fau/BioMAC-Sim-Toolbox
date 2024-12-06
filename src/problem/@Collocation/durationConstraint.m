%======================================================================
%> @file durationConstraint.m
%> @brief Collocation function to compute duration constraint
%> @details
%> Details: Collocation::speedConstraint()
%>
%> @author Anne Koelewijn
%> @date October, 2019
%======================================================================

%======================================================================
%> @brief Computes constraint violation demanding periodic movement
%>
%> @param obj               Collocation class object
%> @param option            String parsing the demanded output
%> @param X                 Double array: State vector containing at least speed of the movement
%> @param targetDuration    Double: target duration to reach
%======================================================================
function output = durationConstraint(obj,option,X,targetDuration)
%% check input parameter
if  ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('State vector X does not contain duration.')
end
    
%% compute demanded output
% forward translation in x direction
duration = X(obj.idx.dur);
% if  isa(targetDuration, 'TrackingData')
%     idxVars = find(ismember(targetDuration.variables.type, 'duration')); %indices of tracking variables
%     targetDuration = targetDuration.variables.mean{idxVars};
% end
if strcmp(option,'confun') %constraints of periodicity constraint
    output =  duration - targetDuration;   
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(1,length(X),obj.Jnnz);
    output(1, obj.idx.dur) = 1;    
else
    error('Unknown option');
end
end

