%% 
%======================================================================
%> @file kneeExtTerm.m
%> @brief Knee extension / kneeinj Term from Veerkamp, 2021
%> @details
%> Details: Collocation::kneeExtTerm()
%> Modified Veerkamps Objective to make it continuously diff'able twice
%>
%> @author Markus Gambietz
%> @date March, 2021
%======================================================================


%======================================================================
%> @brief Matlab function to compute Veerkamp's KneeExt
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'states' of the model
%======================================================================

function output = kneeExtTerm(obj,option,X)

fctname = 'kneeExtTerm';

%% initialization: Return dummy value
if strcmp(option,'init')
    output = NaN;
    % check inputs
    if ~isfield(obj.idx, 'states')
        error('Model states are not stored in the state vector')
    end
    %Throw an error if the used model isn't the 2dc model, since hardcoded
    if isa(obj.model,'Gait2dc')
        % For a 3D implementation, just add the 
    
        obj.objectiveInit.(fctname).idx = [5 8]; %Having right and then left
        obj.objectiveInit.(fctname).rangeR = obj.model.dofs{'knee_angle_r','range'};
        obj.objectiveInit.(fctname).rangeL = obj.model.dofs{'knee_angle_l','range'};
        % phi_0 is 2 degrees, which makes it 2*pi/180, it doesnt matter, approx:
        % simplified: This term uses the same constraints for all joints
        obj.objectiveInit.(fctname).phi_2 = mean([(obj.model.joints.jointD(2)*pi/180)^2; (obj.model.joints.jointD(5)*pi/180)^2]);
    elseif isa(obj.model,'Gait2d_osim')
        obj.objectiveInit.(fctname).idx = [5 8]; %Having right and then left
        obj.objectiveInit.(fctname).rangeR = obj.model.dofs{'knee_angle_r','range_passiveMoment'};
        obj.objectiveInit.(fctname).rangeL = obj.model.dofs{'knee_angle_l','range_passiveMoment'};
        % Transition region is hardcoded as two degrees?
        obj.objectiveInit.(fctname).phi_2 = 2*pi/180;
    else
        error('kneeExtTerm indices not valid non-Gait2dc/osim models')

    end   
    return;
end

kidx = obj.objectiveInit.(fctname).idx;
r_r  = obj.objectiveInit.(fctname).rangeR;
l_r  = obj.objectiveInit.(fctname).rangeL;
phi_2 = obj.objectiveInit.(fctname).phi_2;
nNodes = obj.nNodes;

% Continuous formulation: Punish values greater / smaller than 0
% X_n contains violation values for [r_lower r_upper l_lower l_upper]
X_n = [r_r(1)-X(obj.idx.states(kidx(1),1:nNodes))...
    X(obj.idx.states(kidx(1),1:nNodes))-r_r(2)...
    l_r(1)-X(obj.idx.states(kidx(2),1:nNodes))...
    X(obj.idx.states(kidx(2),1:nNodes))-l_r(2)]';

%% Objectives & Gradients; Copycat from RegTerm
if strcmp(option,'objval')
    % Veerkamp is originally uses abs()
    % Ton's model uses f(x) = 0.5 * (x + sqrt(x^2+x_0^2)); x_0 = 2 degrees
    % x^2 is not applicable, as it isnt twice diff'able in piecewise
    % formulations
    output = sum(sum(0.5 * (X_n + sqrt(X_n.^2+phi_2))))/nNodes;
        
elseif strcmp(option,'gradient')
    output = zeros(size(X));

    grad = 0.5 * (X_n./sqrt(X_n.^2+phi_2) + 1);
    % Lower boundaries are flipped
    grad([1 3],:) = - grad([1 3],:);
    % Assign the gradient to the correct joints
    output(obj.idx.states(kidx(1),1:nNodes)) = sum(grad([1 2],:))/nNodes;
    output(obj.idx.states(kidx(2),1:nNodes)) = sum(grad([3 4],:))/nNodes;
else
    error('Unknown option')
end

end