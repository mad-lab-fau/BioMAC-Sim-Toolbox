%======================================================================
%> @file impactTerm.m
%> @brief Foot-Ground Impact Term from Veerkamp, 2021
%> @details
%> Details: Collocation::impactTerm()
%> Modified Veerkamps Objective to make it continuously diff'able twice
%>
%> @author Markus Gambietz
%> @date March, 2021
%======================================================================


%======================================================================
%> @brief Matlab function to compute Veerkamp's FGImpact
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'states' and 'dur' of the model
%======================================================================

function output = impactTerm(obj,option,X)

fctname = 'impactTerm';

%% initialization: Return dummy value
if strcmp(option,'init')
    output = NaN;
    % check inputs
    if ~isfield(obj.idx, 'states')
        error('Model states are not stored in the state vector')
    end
    if ~isfield(obj.idx, 'dur')
        error('Model duration is not stored in the state vector')
    end
    %Throw an error if the used model isn't the 2dc model, since hardcoded
    if isa(obj.model,'Gait2dc')
        obj.objectiveInit.(fctname).grf_idx = [53 54 57 58 61 62 65 66];
    elseif isa(obj.model,'Gait2d_osim') 
        if  obj.model.nDofs == 9
            obj.objectiveInit.(fctname).grf_idx = [55 56 61 62 67 68 73 74];
        elseif obj.model.nDofs == 10
            obj.objectiveInit.(fctname).grf_idx = [55 56 61 62 67 68 73 74]+2;
        end           
    else
        error('impactTerm is hardcoded for the Gait2dc/osim models')
    end
    return;
end

% Rows of the grf values in standard Gait2dc model
grf_idx = obj.objectiveInit.(fctname).grf_idx;
%h = X(obj.idx.dur)/obj.nNodes;
duration = X(obj.idx.dur);
nNodesDur = obj.nNodesDur;
factor = (nNodesDur-1)/length(grf_idx)/duration^2;
%% Objectives & Gradients
if strcmp(option,'objval')
    output = 0;
    % Veerkamp is originally without the squared, uses abs() instead
    for i = grf_idx
        F = diff(X(obj.idx.states(i,:)));
        output = output + sum(F.^2)*factor;%*X(obj.idx.dur)*X(obj.idx.dur));
    end
    
elseif strcmp(option,'gradient')
    output = zeros(size(X));
    for i = grf_idx
        dGRFdt = diff(X(obj.idx.states(i,:)));
        output(obj.idx.states(i,1:end-1)) = output(obj.idx.states(i,1:end-1)) - 2*dGRFdt*factor;
        output(obj.idx.states(i,2:end)) = output(obj.idx.states(i,2:end)) + 2*dGRFdt*factor;
        output(obj.idx.dur) = output(obj.idx.dur)-2*(nNodesDur-1)/X(obj.idx.dur)^3*(sum(dGRFdt.^2))/length(grf_idx);
    end
    %output = output/obj.nNodes;
else
    error('Unknown option')
end

end