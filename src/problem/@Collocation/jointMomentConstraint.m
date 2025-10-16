%======================================================================
%> @file jointMomentConstraint.m
%> @brief Collocation function to constraint joint moments 
%> @details
%> Details: Collocation::jointMomentConstraint()
%>
%> @author Anne Koelewijn
%> @date September, 2025
%======================================================================

%======================================================================
%> @brief Constraint to ensure joint moments in dynamics are the same as
%> measured/estimated joint moments in trackingData
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'confun' or 'Jacobian'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' of the model
%> @param   data    TrackingData: All rows with type 'moment' will be
%>                  constrained
%>
%> @retval  output  vector with differences between joint moments from 
%>                  simulation and measurement, or corresponding Jacobian
%======================================================================
function output = jointMomentConstraint(obj,option,X,momData, angData, dur, isPeriodic,torque_names)

fctname = 'jointMomentConstraint';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states_mus') % check whether model states are stored in X
        error('Only muscle states should be stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    obj.constraintInit.(fctname).variables = momData.variables;
    obj.constraintInit.(fctname).nVars = height(momData.variables); % total number of variables
    obj.constraintInit.(fctname).measMean = cell2mat(momData.variables.mean'); % mean measured data
    nNodesDur = obj.nNodesDur;
    h = dur/(nNodesDur-1);
    for iAng = 1:height(angData.variables)
        obj.constraintInit.(fctname).angles(iAng,:) =  angData.variables.mean{iAng};
    end
    if isPeriodic
        obj.constraintInit.(fctname).angles(:,end+1) = obj.constraintInit.(fctname).angles(:,1);
    end
    obj.constraintInit.(fctname).angularvelocities = diff(obj.constraintInit.(fctname).angles,[],2)/h;
    if isPeriodic
        obj.constraintInit.(fctname).angularvelocities(:,end+1) = obj.constraintInit.(fctname).angularvelocities(:,1);
    else
        obj.constraintInit.(fctname).angularvelocities = [zeros(height(angData.variables),1) obj.constraintInit.(fctname).angularvelocities];
    end
    obj.constraintInit.(fctname).idxInModel_ang = obj.model.extractState('q', angData.variables.name); % model indices of joint angles
    obj.constraintInit.(fctname).idxInModel_avl = obj.model.extractState('qdot', angData.variables.name); % model indices of angular velocities
    obj.constraintInit.(fctname).idxInModel_mus = sort([obj.model.extractState('a');obj.model.extractState('s')]); % model indices of muscle variables
    iDx = false(obj.model.nDofs,1);
    for iVar = 1:height(momData.variables)
        iDx = or(iDx,strcmp(obj.model.dofs.Properties.RowNames,momData.variables.name{iVar}));
    end
    obj.constraintInit.(fctname).idxInM = find(iDx == 1); % model indices of muscle variables

    % get factor to convert simulated data to unit of tracking data
    % (assuming that all angles and angular velocities have the same unit)
    if strcmp(angData.variables.unit{1}, 'rad')
        obj.constraintInit.(fctname).factorToMeas = 1;
    elseif strcmp(angData.variables.unit{1}, 'deg')
        obj.constraintInit.(fctname).factorToMeas = pi/180;
    end

    if nargin == 8
        obj.constraintInit.(fctname).useMom = 1;
        torqueData = momData.extractData('moment',torque_names);
        input_moments = torqueData.variables;
        %torque_names);
        for iMom = 1:height(input_moments)
            obj.constraintInit.(fctname).input_moments(iMom,:) = input_moments.mean{iMom}/obj.model.mExtraScaleFactor;
        end
        if isPeriodic
            obj.constraintInit.(fctname).input_moments(:,end+1) = obj.constraintInit.(fctname).input_moments(:,1);
        end
        obj.constraintInit.(fctname).idxInControl_mom = obj.model.extractControl('torque', torque_names); % model indices of joint moments
    else
        obj.constraintInit.(fctname).useMom = 0;
    end
    
    % Return a dummy value
    output = NaN;
    return;
end

%% compute demanded output
nNodesDur = obj.nNodesDur;

% get variables from initialization (faster)
nVars = obj.constraintInit.(fctname).nVars;
measMean = obj.constraintInit.(fctname).measMean;

angles = obj.constraintInit.(fctname).angles;
angularvelocities = obj.constraintInit.(fctname).angularvelocities;
factorToMeas = obj.constraintInit.(fctname).factorToMeas;
idxInModel_ang = obj.constraintInit.(fctname).idxInModel_ang;
idxInModel_avl = obj.constraintInit.(fctname).idxInModel_avl;
idxInModel_mus = obj.constraintInit.(fctname).idxInModel_mus;
idxInM = obj.constraintInit.(fctname).idxInM;

if isPeriodic
    nNodesEnd = nNodesDur;
else
    nNodesEnd = nNodesDur - 1;
end

useMom = obj.constraintInit.(fctname).useMom;
if useMom
    input_moments = obj.constraintInit.(fctname).input_moments;
    idxInControl_mom = obj.constraintInit.(fctname).idxInControl_mom;
end


norm_factor = obj.model.bodymass*norm(obj.model.gravity);

if strcmp(option,'confun')
    output = zeros(nVars*nNodesEnd,1);
        
    % joint moments must equal measured/estimated ones
    for iNode=1:nNodesEnd
        ic = (1:nVars) +  (iNode-1)*nVars; %indices of constraints of iNode in c
        x = zeros(obj.model.nStates,1);
        x(idxInModel_ang) = angles(:,iNode)*factorToMeas;
        x(idxInModel_avl) = angularvelocities(:,iNode)*factorToMeas;
        x(idxInModel_mus) = X(obj.idx.states_mus(:,iNode));
        u = X(obj.idx.controls(:,iNode+1));
        if useMom
            u(idxInControl_mom) = input_moments(:,iNode); %% iNode or iNode+1
        end
        M = obj.model.getJointmoments(x, u);
        
        output(ic) = (M(idxInM)-measMean(iNode,:)')/norm_factor;
    end
elseif strcmp(option,'jacobian')
    output = spalloc(nVars*nNodesEnd,length(X),obj.Jnnz);
    
    for iNode = 1:nNodesEnd
        ic = (1:nVars) +  (iNode-1)*nVars; %indices of constraints of iNode in c
        x = zeros(obj.model.nStates,1);
        x(idxInModel_ang) = angles(:,iNode)*factorToMeas;
        x(idxInModel_avl) = angularvelocities(:,iNode)*factorToMeas;
        x(idxInModel_mus) = X(obj.idx.states_mus(:,iNode));
        u = X(obj.idx.controls(:,iNode+1));
        if useMom
            u(idxInControl_mom) = input_moments(:,iNode); %% iNode or iNode+1
        end

        [~, dMdx, dMdu] = obj.model.getJointmoments(x, u);
        
        output(ic, obj.idx.states_mus(:,iNode)) = dMdx(idxInModel_mus,idxInM)'/norm_factor;
        if useMom
            dMdu(idxInControl_mom,:) = []; %only use muscle states for derivative.
        end
        output(ic, obj.idx.controls(:,iNode)) = dMdu(:,idxInM)'/norm_factor;
    end
else
    error('Unknown option.');
end
end
