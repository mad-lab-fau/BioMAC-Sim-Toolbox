%======================================================================
%> @file headStabTerm.m
%> @brief Head Stabilistation Term from Veerkamp, 2021
%> @details
%> Details: Collocation::headStabTerm()
%> Modified Veerkamps Objective to make it continuously diff'able twice
%>
%> @author Markus Gambietz
%> @date March, 2021
%======================================================================


%======================================================================
%> @brief Matlab function to compute Veerkamp's HeadStab
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'states' and 'dur' of the model
%======================================================================

function output = headStabTerm(obj,option,X)

fctname = 'headStabTerm';

%% initialization: Return dummy value
if strcmp(option,'init')
    output = NaN;
    % check inputs
    if ~isfield(obj.idx, 'states')
        error('Model states are not stored in the state vector')
    end
    %Throw an error if the used model isn't the 2dc model, since hardcoded
    if isa(obj.model,'Gait2dc')
        obj.objectiveInit.(fctname).length = obj.model.segments{'pelvis','length'};           
    elseif isa(obj.model,"Gait2d_osim") && obj.model.nDofs == 9
        % length is implemented as twice the CoM y distance from the pelvis
        obj.objectiveInit.(fctname).length = 2*obj.model.segments{'torso','mass_center'}(2)+obj.model.joints{'back','location'}(2);           
    else 
        error('headStabTerm is not yet implemented for non-lumbarlock Gait2dc/osim models')
    end   
    return;
end

nNodesDur = obj.nNodesDur;
nDims = 2;
dxIdx = 10;
dyIdx = 11;
tiltIdx = 3;
dtiltIdx = 12;
len = obj.objectiveInit.(fctname).length;

% Head accelerations in x & y
dq_x = diff(X(obj.idx.states(dxIdx,:)));
dq_y = diff(X(obj.idx.states(dyIdx,:)));
dtilt = diff(X(obj.idx.states(dtiltIdx,:)));
tilt = X(obj.idx.states(tiltIdx,1:end-1));

a_x = dq_x - cos(tilt) .* dtilt * len;
a_y = dq_y - sin(tilt) .* dtilt * len;


duration = X(obj.idx.dur);
factor = (nNodesDur-1)/nDims/duration^2;

%% Objectives & Gradients; Copycat from RegTerm
if strcmp(option,'objval')
    output = 0;
    % Veerkamp is originally without the squared, uses abs() instead
    F = a_x.^2 + a_y.^2;
    output = output + sum(F)*factor;
    
elseif strcmp(option,'gradient')
    output = zeros(size(X));
    % x - coordinates
    output(obj.idx.states(dxIdx,1:end-1)) =  - 2*factor*(dq_x- cos(tilt) .* dtilt * len);
    output(obj.idx.states(dxIdx,2:end)) = output(obj.idx.states(dxIdx,2:end))   + 2*factor*(dq_x- cos(tilt) .* dtilt * len);
    % y - coordinates
    output(obj.idx.states(dyIdx,1:end-1)) =  -  2*factor*(dq_y- sin(tilt) .* dtilt * len);
    output(obj.idx.states(dyIdx,2:end)) = output(obj.idx.states(dyIdx,2:end))   + 2*factor*(dq_y- sin(tilt) .* dtilt * len);
    
    % tilt: coordinates
    output(obj.idx.states(tiltIdx,1:end-1)) =  2*factor*((dq_x - cos(tilt) .* dtilt * len)*len.*dtilt.*sin(tilt) - ...
        (dq_y - sin(tilt) .* dtilt * len)*len.*dtilt.*cos(tilt));
    
    % d_tilt: coordinates
    output(obj.idx.states(dtiltIdx,1:end-1)) = 2*factor*((dq_x - cos(tilt) .* dtilt * len)*len.*cos(tilt) + ...
        (dq_y - sin(tilt) .* dtilt * len)*len.*sin(tilt));
    output(obj.idx.states(dtiltIdx,2:end)) = output(obj.idx.states(dtiltIdx,2:end)) - ...
        2*factor*((dq_x - cos(tilt) .* dtilt * len )*len.*cos(tilt) + ...
        (dq_y - sin(tilt) .* dtilt * len)*len.*sin(tilt));
    
    % duration:
    output(obj.idx.dur) = output(obj.idx.dur)-2*(nNodesDur-1)/duration^3*(sum(a_x.^2+a_y.^2))/nDims;
else
    error('Unknown option')
end

end