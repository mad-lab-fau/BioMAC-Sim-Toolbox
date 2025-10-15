%======================================================================
%> @file periodicityConstraint.m
%> @brief Collocation function to compute periodicity constraint for muscle
%> states only
%> @details
%> Details: Collocation::musclePeriodicityConstraint()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding periodic movement when
%> only evaluating muscle dynamics
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = musclePeriodicityConstraint(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states_mus')
    error('State vector X does not contain required states.')
end
    
%% compute demanded output
nStates = size(obj.idx.states_mus,1);
nControls = size(obj.idx.controls,1);


if strcmp(option,'init') %constraints of periodicity constraint
    output = nan;
elseif strcmp(option,'confun') %constraints of periodicity constraint
    output = zeros(nStates+nControls,1);
    output(1:nStates) = X(obj.idx.states_mus(:,end)) - X(obj.idx.states_mus(:,1));
    output(nStates+1:nStates+nControls) = X(obj.idx.controls(:,end)) - X(obj.idx.controls(:,1));
    
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(nStates+nControls,length(X),obj.Jnnz); %> @todo where to get Jnnz
    
    output(1:nStates,obj.idx.states_mus(:,end))   = speye(nStates);
    output(1:nStates,obj.idx.states_mus(:,1))     = -speye(nStates);

    output(nStates+1:nStates+nControls,obj.idx.controls(:,end))   = speye(nControls);
    output(nStates+1:nStates+nControls,obj.idx.controls(:,1))     = -speye(nControls);       
    
else
    error('Unknown option');
end
end

