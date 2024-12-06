%======================================================================
%> @file fromStandingConstraint.m
%> @brief Collocation function to find walking/running from standing
%> @details
%> Details: Collocation::fromStandingConstraint()
%>
%> @author Anne Koelewijn
%> @date July, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to find walking/running from standing
%>
%> @todo Generalize this function and remove hard coded indices. E.g. we do
%> not assume fixed indices for the states. The function should also work
%> for 2D.
%>
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = fromStandingConstraint(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'speed') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('State vector X does not contain required states.')
end

if size(X(obj.idx.states),2) ~= (obj.nNodes+1) || size(X(obj.idx.controls),2) ~= (obj.nNodes+1) %> @todo should we add another field to the state vector instead of expanding states to obj.nNodes+1
    error('Model states and controls need to be optimized at collocation nodes + 1')
end
    
%% compute demanded output
nStates = size(obj.idx.states,1);
nControls =  size(obj.idx.controls,1);
icx = 1:nStates; % first indices are for the states
% icu = (nStates+1):(nStates+nControls); % next indices are for the controls

% forward translation in x direction
speed = X(obj.idx.speed);
dur = X(obj.idx.dur);
unitdisplacement = zeros(nStates,1);
unitdisplacement(obj.model.idxForward) = 1;

if strcmp(option,'confun') %constraints of periodicity constraint
    output = zeros(nStates+1,1);
    
    % state must start at zeros, upright position
    output(icx) = X(obj.idx.states(:,1));
    output(icx(5)) = output(icx(5))-0.9549; %y should be at right height from standing solution
    output(end) = X(obj.idx.states(3,end)) - 1; %1 meter displacement
    
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(nStates+1,length(X),obj.Jnnz); %> @todo where to get Jnnz
    
    output(icx,obj.idx.states(:,1)) = speye(nStates);
    output(end,obj.idx.states(3,end)) = 1;    
else
    error('Unknown option');
end
end

