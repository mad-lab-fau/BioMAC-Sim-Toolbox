%======================================================================
%> @file @Collocation/muscleDynamicConstraints.m
%> @brief Collocation function to compute dynamic constraint
%> @details
%> Details: Collocation::muscleDynamicConstraints()
%>
%> @author Anne Koelewijn
%> @date September, 2025
%======================================================================

%======================================================================
%> @brief Computes constaint violation of dynamics violation for muscle
%> dynamics only. This function is intended to be used for a dynamic
%> optimization of the muscle states, using joint angles and joint moments
%> from an experiment (e.g., created with inverse kinematics and dynamics)
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%> @param angleData     TrackingData: the full dynamics state will be
%>                      created using the 'angles' from the data
%> @param dur           Double: motion duration used to calculate time step
%> @param isPeriodic    Boolean: 1 is motion is periodic, 0 otherwise
%> @param momData       (optional) TrackingData: joint moments if they are to be directly applied to the model (e.g., gait3d_pelvis213)
%======================================================================
function output = muscleDynamicConstraints(obj,option,X,angleData, dur, isPeriodic, momData)
%% check input parameter
if ~isfield(obj.idx,'states_mus') || ~isfield(obj.idx,'controls') % check whether controls are stored in X
    error('Model states and controls and duration need to be stored in state vector X.')
end

fctname = 'muscleDynamicConstraints';

if nargin < 6
    isPeriodic = 0;
end

if strcmp(obj.Euler,'SIE')
    %% state variable indices for semi-implicit Euler method
    % https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
    % Discretization is defined by ixSIE1 (state variables from the node before the time step)
    % and ixSIE2 (state variables from the node after the time step).
    % Two versions: SIE A and SIE B. Both conserve energy.
    % SIE A appeared easier to solve in tests.  This makes sense because it uses backward
    % Euler for muscle equations. SIE B code is commented out and provided for comparison.
    % SIE A: 
    ixSIE2 = sort([obj.model.extractState('qdot'); obj.model.extractState('s'); obj.model.extractState('a')]);
    ixSIE1 = setdiff(1:obj.model.nStates, ixSIE2);
    % SIE B
    % ixSIE1 = sort([obj.model.extractState('q')]);
    % ixSIE2 = setdiff(1:obj.model.nStates, ixSIE1);
end

%% initalization
if strcmp(option,'init')
    % initialize some variables (faster to get it once)
    angles = angleData.variables(strcmp(angleData.variables.type, 'angle'), :);
    nNodesDur = obj.nNodesDur;
    obj.constraintInit.(fctname).h = dur/(nNodesDur-1);
    for iAng = 1:height(angles)
        obj.constraintInit.(fctname).angles(iAng,:) = angles.mean{iAng};
    end
    if isPeriodic
        obj.constraintInit.(fctname).angles(:,end+1) = obj.constraintInit.(fctname).angles(:,1);
    end
    obj.constraintInit.(fctname).angularvelocities = diff(obj.constraintInit.(fctname).angles,[],2)/obj.constraintInit.(fctname).h;
    if isPeriodic
        obj.constraintInit.(fctname).angularvelocities(:,end+1) = obj.constraintInit.(fctname).angularvelocities(:,1);
    else
        obj.constraintInit.(fctname).angularvelocities = [zeros(height(angles),1) obj.constraintInit.(fctname).angularvelocities];
    end
    obj.constraintInit.(fctname).idxInModel_ang = obj.model.extractState('q', angles.name); % model indices of joint angles
    obj.constraintInit.(fctname).idxInModel_avl = obj.model.extractState('qdot', angles.name); % model indices of angular velocities
    obj.constraintInit.(fctname).idxInModel_mus = sort([obj.model.extractState('a');obj.model.extractState('s')]); % model indices of muscle variables

    % get factor to convert simulated data to unit of tracking data
    % (assuming that all angles and angular velocities have the same unit)
    if strcmp(angles.unit{1}, 'rad')
        obj.constraintInit.(fctname).factorToMeas = 1;
    elseif strcmp(angles.unit{1}, 'deg')
        obj.constraintInit.(fctname).factorToMeas = pi/180;
    end

    if nargin == 7
        obj.constraintInit.(fctname).useMom = 1;
        moments = momData.variables(strcmp(momData.variables.type, 'moment'), :);
        for iMom = 1:height(moments)
            obj.constraintInit.(fctname).moments(iMom,:) = moments.mean{iMom}/obj.model.mExtraScaleFactor;
        end
        if isPeriodic
            obj.constraintInit.(fctname).moments(:,end+1) = obj.constraintInit.(fctname).moments(:,1);
        end
        obj.constraintInit.(fctname).idxInControl_mom = obj.model.extractControl('torque', moments.name); % model indices of joint moments
    else
        obj.constraintInit.(fctname).useMom = 0;
    end

    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
% get variables from initalization (faster)
angles = obj.constraintInit.(fctname).angles;
angularvelocities = obj.constraintInit.(fctname).angularvelocities;
factorToMeas = obj.constraintInit.(fctname).factorToMeas;
idxInModel_ang = obj.constraintInit.(fctname).idxInModel_ang;
idxInModel_avl = obj.constraintInit.(fctname).idxInModel_avl;
idxInModel_mus = obj.constraintInit.(fctname).idxInModel_mus;

nNodesDur = obj.nNodesDur;
h = obj.constraintInit.(fctname).h;
nconstraintspernode = length(idxInModel_mus);

useMom = obj.constraintInit.(fctname).useMom;
if useMom
    moments = obj.constraintInit.(fctname).moments;
    idxInControl_mom = obj.constraintInit.(fctname).idxInControl_mom;
end

if isPeriodic
    nNodesEnd = nNodesDur;
else
    nNodesEnd = nNodesDur - 1;
end
if strcmp(option,'confun')
    output = zeros(nconstraintspernode*nNodesEnd,1);
        
    % dynamic equations must be zero
   
    for iNode=1:nNodesEnd
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        x1 = zeros(obj.model.nStates,1);
        x1(idxInModel_ang) = angles(:,iNode)*factorToMeas;
        x1(idxInModel_avl) = angularvelocities(:,iNode)*factorToMeas;
        x1(idxInModel_mus) = X(obj.idx.states_mus(:,iNode));
        
        x2 = zeros(obj.model.nStates,1);
        x2(idxInModel_ang) = angles(:,iNode+1)*factorToMeas;
        x2(idxInModel_avl) = angularvelocities(:,iNode+1)*factorToMeas;
        x2(idxInModel_mus) = X(obj.idx.states_mus(:,iNode+1));

        xd =(x2-x1)/h;
        u = X(obj.idx.controls(:,iNode+1));
        if useMom
            u(idxInControl_mom) = moments(:,iNode+1);
        end
        
        if strcmp(obj.Euler,'BE')
            f = obj.model.getDynamics(x2,xd,u);	% backward Euler discretization
        elseif strcmp(obj.Euler,'ME')
            % we're using u2 instead of (u1+u2)/2 because u is
            % open loop and it makes no difference except u2
            % will converge more easily
            f = obj.model.getDynamics((x1+x2)/2,xd,u);
        elseif strcmp(obj.Euler,'SIE')
            % semi-implicit Euler is energy neutral
			% take specified states ixSIE2 from x2, all other states from x1
			x1(ixSIE2) = x2(ixSIE2);
			f = obj.model.getDynamics(x1,xd,u);
        end
        output(ic) = f(idxInModel_mus);
    end
elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode*nNodesEnd,length(X),obj.Jnnz);
    
    for iNode = 1:nNodesEnd
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        
        ix1 = obj.idx.states_mus(:,iNode);
        ix2 = obj.idx.states_mus(:,iNode+1);
        iu = obj.idx.controls(:,iNode+1);
        
        x1 = zeros(obj.model.nStates,1);
        x1(idxInModel_ang) = angles(:,iNode)*factorToMeas;
        x1(idxInModel_avl) = angularvelocities(:,iNode)*factorToMeas;
        x1(idxInModel_mus) = X(obj.idx.states_mus(:,iNode));
        
        x2 = zeros(obj.model.nStates,1);
        x2(idxInModel_ang) = angles(:,iNode+1)*factorToMeas;
        x2(idxInModel_avl) = angularvelocities(:,iNode+1)*factorToMeas;
        x2(idxInModel_mus) = X(obj.idx.states_mus(:,iNode+1));

        xd =(x2-x1)/h;
        u = X(obj.idx.controls(:,iNode+1));
        if useMom
            u(idxInControl_mom) = moments(:,iNode+1);
        end

        if strcmp(obj.Euler,'BE')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x2,xd,u);
            output(ic,ix1) = -dfdxdot(idxInModel_mus,idxInModel_mus)'/h;
            output(ic,ix2) = dfdx(idxInModel_mus,idxInModel_mus)' + dfdxdot(idxInModel_mus,idxInModel_mus)'/h;
        elseif strcmp(obj.Euler,'ME')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics((x1+x2)/2,xd,u);
            output(ic,ix1) = dfdx(idxInModel_mus,idxInModel_mus)'/2 - dfdxdot(idxInModel_mus,idxInModel_mus)'/h;
            output(ic,ix2) = dfdx(idxInModel_mus,idxInModel_mus)'/2 + dfdxdot(idxInModel_mus,idxInModel_mus)'/h;
        elseif strcmp(obj.Euler,'SIE')
            error('Derivates have not been tested yet')
			x1(ixSIE2) = x2(ixSIE2);
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x1,xd,u);
            output(ic,ix1) = -dfdxdot(idxInModel_mus,idxInModel_mus)'/h;
            output(ic,ix2) = dfdxdot(idxInModel_mus,idxInModel_mus)'/h;
            output(ic,ix1(ixSIE1)) = output(ic,ix1(ixSIE1)) + dfdx(ixSIE1,:)';
            output(ic,ix2(ixSIE2)) = output(ic,ix2(ixSIE2)) + dfdx(ixSIE2,:)';
        end
        if useMom
            dfdu(idxInControl_mom,:) = []; %only use muscle states for derivative.
        end
        output(ic,iu) = dfdu(:,idxInModel_mus)';
        
    end
else
    error('Unknown option.');
end
end
