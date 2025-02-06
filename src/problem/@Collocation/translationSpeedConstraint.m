%======================================================================%======================================================================
%> @file translationSpeedConstraint.m
%> @brief Collocation function to compute a constraint based on speed and translation
%> @details
%> Details: Collocation::translationSpeedConstraint()
%>
%> @author Marlies Nitschke
%> @date March, 2019
%======================================================================

%======================================================================
%> @brief Computes constraint based on speed and translation
%>
%> @details
%> Constraint: deltaTranslation/duration - speed = 0
%> 
%> If two entries of type speed are contained in X, it assumes the first
%> one to be forward (x-direction) and the second one to be sideward (z-direction).
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least states, speed and duration of the movement
%======================================================================
function output = translationSpeedConstraint(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'speed') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('State vector X does not contain required states.')
end

if ~obj.isPeriodic
    error('Only implemented for periodic movements.')
end

%% compute demanded output
% get speed and duration
speed = X(obj.idx.speed);
nSpeeds = numel(speed); 
duration = X(obj.idx.dur);

% forward translation in x direction
idxForward = intersect(obj.model.idxForward, obj.model.extractState('q'));
% sideward translation in z direction
idxSideward = intersect(obj.model.idxSideward, obj.model.extractState('q'));

if strcmp(option,'confun') %constraints of periodicity constraint
    output = zeros(nSpeeds,1);
    
    output(1) = (X(obj.idx.states(idxForward, end)) - X(obj.idx.states(idxForward, 1))) / duration - speed(1);
    
    if nSpeeds > 1
        output(2) = (X(obj.idx.states(idxSideward, end)) - X(obj.idx.states(idxSideward, 1))) / duration - speed(2);
    end
    
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(nSpeeds,length(X),obj.Jnnz); 
    
    output(1, obj.idx.states(idxForward, end)) = 1/duration;
    output(1, obj.idx.states(idxForward, 1)) = -1/duration;
    output(1, obj.idx.dur) = -(X(obj.idx.states(idxForward, end)) - X(obj.idx.states(idxForward, 1))) /duration^2;
    output(1, obj.idx.speed(1)) = -1;
    
    if nSpeeds > 1
        output(2, obj.idx.states(idxSideward, end)) = 1/duration;
        output(2, obj.idx.states(idxSideward, 1)) = -1/duration;
        output(2, obj.idx.dur) = -(X(obj.idx.states(idxSideward, end)) - X(obj.idx.states(idxSideward, 1))) /duration^2;
        output(2, obj.idx.speed(2)) = -1;
    end
    
else
    error('Unknown option');
end
end
