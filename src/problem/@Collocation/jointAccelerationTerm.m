%======================================================================
%> @file jointAccelerationTerm.m
%> @brief Collocation function to compute the joint acceleration from states vector
%> @details
%> Details: Collocation::jointAccelerationTerm()
%>
%> @author Anne Koelewijn
%> @date May, 2020
%======================================================================

%======================================================================
%> @brief Matlab function to compute effort of the muscles from states vector
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'states' of the model
%======================================================================
function output = jointAccelerationTerm(obj,option,X)
fctname = 'jointAccelerationTerm';

%% check input parameter
if strcmp(option,'init')
    if ~isfield(obj.idx,'states')% check whether states are stored in X
        error('Model states need to be stored in state vector X.')
    end
    
    % Get indices of qdot of all DoFs of joints
    globalDoFNames = {'pelvis_tx', 'pelvis_ty', 'pelvis_tz', 'pelvis_tilt', 'pelvis_list','pelvis_rotation','pelvis_obliquity'};
    DoFNames = obj.model.dofs.Properties.RowNames;
    jointDoFsNames = DoFNames(~ismember(DoFNames, globalDoFNames));
    obj.objectiveInit.(fctname).idxqdotJoint = obj.model.extractState('qdot', jointDoFsNames);
   
    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
idxqdotJoint = obj.objectiveInit.(fctname).idxqdotJoint;
xd = diff(X(obj.idx.states(idxqdotJoint,:)),1,2);

% Normalisation factor
duration = X(obj.idx.dur);
nDims = length(idxqdotJoint);
factor = obj.nNodes/nDims/duration^2;


if strcmp(option,'objval') 
    output = sum(sum(xd.^2)) * factor; % make it the average
elseif strcmp(option,'gradient')
    output = zeros(size(X));
    % Coordinates
    output(obj.idx.states(idxqdotJoint,1:end-1))   = output(obj.idx.states(idxqdotJoint,1:end-1))   - 2*factor*xd;
    output(obj.idx.states(idxqdotJoint,2:end))     = output(obj.idx.states(idxqdotJoint,2:end))     + 2*factor*xd;

    % Duration
    output(obj.idx.dur) = output(obj.idx.dur)-2*obj.nNodes/duration^3*(sum(sum(xd.^2)))/nDims;
else
    error('Unknown option.');
end