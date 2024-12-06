%======================================================================
%> @file trackGRF.m
%> @brief Collocation function to track ground reaction forces
%> @details
%> Details: Collocation::trackGRF()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Matlab function to track ground reaction forces
%>
%> @details 
%> Computes difference between simulated and measured ground
%> reaction forces.
%> Supports data in bodyweight ('BW'), bodyweight percentage ('BW%'), and newton ('N').
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' of the model
%> @param   data    TrackingData: All rows with type 'GRF' will be tracked
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient for input option 'gradient'
%======================================================================
function output = trackGRF(obj,option,X,data)

fctname = 'trackGRF';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states') % check whether controls are stored in X
        error('Model states are not stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    variables = data.variables(strcmp(data.variables.type, 'GRF'), :); % data with the correct type
    obj.objectiveInit.(fctname).nVars = height(variables); % total number of variables
    obj.objectiveInit.(fctname).measMean = cell2mat(variables.mean'); % mean measured data
    if data.nTrials == 1
        obj.objectiveInit.(fctname).measVar = ones(size(obj.objectiveInit.(fctname).measMean)); % use 1 to operate as if we are not dividing by the variance
    else
        obj.objectiveInit.(fctname).measVar = cell2mat(variables.var'); % variance of measured data
    end
    obj.objectiveInit.(fctname).names = variables.name; % names of variables

    % get factor to convert simulated data to unit of tracking data
    % (assuming that all angles have the same unit)
    if strcmp(variables.unit{1}, 'BW')
        obj.objectiveInit.(fctname).factorToMeas = 1;
    elseif strcmp(variables.unit{1}, 'BW%')
        obj.objectiveInit.(fctname).factorToMeas = 100;
    elseif strcmp(variables.unit{1}, 'N')
        obj.objectiveInit.(fctname).factorToMeas = obj.model.bodymass*norm(obj.model.gravity);
    end

    % Return a dummy value
    output = NaN;
    return;
end


%% compute demanded output
% get variables from initalization (faster)
nVars = obj.objectiveInit.(fctname).nVars;
measMean = obj.objectiveInit.(fctname).measMean;
measVar = obj.objectiveInit.(fctname).measVar;
factorToMeas = obj.objectiveInit.(fctname).factorToMeas;
names = obj.objectiveInit.(fctname).names;

% initialize some more variables
x = X(obj.idx.states); % extract states

if strcmp(option,'objval') %objective value for tracking GRFs  
    output = 0;
    GRF = zeros(obj.nNodes, 12);
    for iNode = 1:obj.nNodes %accumulate tracking variables for nNodes
        GRF(iNode,:) = obj.model.getGRF(x(:,iNode));
    end
    
    for iVar = 1:nVars
        switch names{iVar}
            case 'GRF_x_r'
                simVar = GRF(:, 1);
            case 'GRF_y_r'
                simVar = GRF(:, 2);
            case 'GRF_z_r'
                simVar = GRF(:, 3);
            case 'GRF_x_l'
                simVar = GRF(:, 7);
            case 'GRF_y_l'
                simVar = GRF(:, 8);
            case 'GRF_z_l'
                simVar = GRF(:, 9);
        end
        output = output + ...
            sum((simVar*factorToMeas - measMean(:, iVar)).^2./measVar(:, iVar)) / nVars /obj.nNodes;
      
    end

elseif strcmp(option,'gradient') %gradient for tracking GRFs 
    % track GRF
    output = zeros(size(X));
    if nVars % calculate it only ones
        GRF = zeros(obj.nNodes, 12);
        dGRFdx = zeros(obj.model.nStates,obj.nNodes,12);
        for iNode = 1:obj.nNodes
            [GRF(iNode,:),dGRFdxtmp] = obj.model.getGRF(x(:,iNode));
            dGRFdx(:,iNode,:) = full(dGRFdxtmp); % it is only sparse for 3D; not sure if we need this
        end
    end
    for iVar = 1:nVars
        switch names{iVar}
            case 'GRF_x_r'
                simVar = GRF(:, 1);
                dsimVar_dx = dGRFdx(:,:,1);
            case 'GRF_y_r'
                simVar = GRF(:, 2);
                dsimVar_dx = dGRFdx(:,:,2);
            case 'GRF_z_r'
                simVar = GRF(:, 3);
                dsimVar_dx = dGRFdx(:,:,3);
            case 'GRF_x_l'
                simVar = GRF(:, 7);
                dsimVar_dx = dGRFdx(:,:,7);
            case 'GRF_y_l'
                simVar = GRF(:, 8);
                dsimVar_dx = dGRFdx(:,:,8);
            case 'GRF_z_l'
                simVar = GRF(:, 9);
                dsimVar_dx = dGRFdx(:,:,9);
        end
        df4_dsimVar = 2*(simVar*factorToMeas - measMean(:, iVar))./measVar(:, iVar) / nVars / obj.nNodes;
        output(obj.idx.states(:,1:obj.nNodes)) = output(obj.idx.states(:,1:obj.nNodes)) + repmat(df4_dsimVar',obj.model.nStates,1).*dsimVar_dx*factorToMeas; %chain rule
        
    end
   
else
    error('Unknown option');
end

end