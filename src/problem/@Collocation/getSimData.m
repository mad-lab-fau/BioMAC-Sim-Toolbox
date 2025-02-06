%======================================================================
%> @file getSimData.m
%> @brief Collocation function to extract simulated variables
%> @details
%> Details: Collocation::getSimData()
%>
%> @author Eva Dorschky, Marlies Nitschke
%> @date June, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to extract simulated variables
%>
%> @details
%> Warning: This function was not tested completely!
%>
%> Currently, this function can return simulated data for the following
%> data types:
%>   - translation
%>   - angle
%>   - qdot
%>   - moment
%>   - GRF
%>   - CoP
%>   - marker
%>   - acc
%>   - gyro
%>   - a
%>   - s
%>   - LMTU
%>   - LCE
%>   - LSEE
%>   - LdotMTU
%>   - LdotCE
%>   - LdotSEE
%>   - muscleForce
%>   - CEForce
%>   - u
%>   - speed
%>   - duration
%>   - footAngle
%>   - musclePower
%>   - CEPower
%>   - SEEPower
%>   - torque
%>   - sdot
%>   - muscleMetRate
%>   - CoM
%>   - jointPower
%>
%> @todo Check which unit is requested and convert the data if needed.
%> However this is currently only done for angles which are automatically
%> converted into degrees if the unit is 'deg'.
%>
%> @todo Use model.GRFNAMES instead of hard coded indices. Eventally change
%> the names in GRFNAMES
%>
%> @todo Check again if we use always the correct discretization method. If
%> we use BE, shouldn't we use states(idxSimVar, 2:nNodes+1) and controls(idxSimVar, 2:nNodes+1)?
%>
%> @param  obj        Collocation class object
%> @param  X          Double matrix: State vector (i.e. result) of the problem
%> @param  simVar     Table: Information on variables from which we want to
%>                    get the simulated data. This table has to contain
%>                    the columns: type, name, and for IMU also position and direction. 
%> @param  modelName  For metabolic rate: name of the model that should be
%>                    used
%> @retval simVar     Table: Input table with an additional column "sim" 
%>                    containing the simulated data for this X.
%======================================================================
function simVar = getSimData(obj,X, simVar, modelName)

%% Check whether the simVar is not emtpy
if isempty(simVar)
    msg = ['simVar is empty. Please, provide a table with informtion on simulated data. ', ...
        'This has to contain at least the columns type, name, unit. ', ...
        'If IMU data is requested, also position and direction have to be given.'];
    error('Collocation:getSimData', msg);
end

%% Extract parameters which are needed
nNodes = obj.nNodes;
nNodesDur = obj.nNodesDur;
model  = obj.model;

%% Check if this function supportes all the requested types in simVar and what data we need to extract
% Define which types are supported by this function
supportedTypes = {'translation', 'angle', 'moment', 'GRF', 'acc', 'gyro', ...
                  'a', 'muscleForce', 'u', 'speed', 'duration', 's', 'qdot', ...
                  'footAngle', 'CEPower', 'torque', 'sdot', 'muscleMetRate', ...
                  'CoM', 'jointPower', 'CEForce', 'SEEPower', 'musclePower', ...
                  'LdotMTU', 'LdotCE', 'LdotSEE', 'LCE', 'LMTU', 'LSEE', 'marker', 'CoP'};

% Define which data has to be computed to get the simulated data for each
% of the supported types
% -> Indices of types in supportedTyped which need the derivative of the entire state vector; 
idxTypesNeedStated = [5, 15, 17, 18, 21, 22, 25, 26];
% -> Indices of types in supportedTypes which need the state vector;
% This is also needed to get q, qd, a, s
idxTypesNeedState = unique([1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 31, idxTypesNeedStated]);
% -> Indices of types in supportedTyped which need the control vector
idxTypesNeedControl = [3, 9, 16, 18, 20];
% -> Indices of types in supportedTypes which need duration; Duration is
% needed to compute xd
idxTypesNeedDur   = unique([11, idxTypesNeedStated]);

% Check which types are contained in the simulated data and how often
simTypes = unique(simVar.type);
nEachType = zeros(size(simTypes));
for iType = 1 : length(simTypes)
   nEachType(iType) = numel(find(strcmp(simVar.type, simTypes{iType})));
end
idxSimTypes = find(ismember(supportedTypes, simTypes)); % index of simTypes in supportedTypes

% Check whether there are unsupported types in simVar
unsupportedTypes = simTypes(~ismember(simTypes, supportedTypes));
if ~isempty(unsupportedTypes)
    msgUnSupTypes = sprintf('''%s'', ', unsupportedTypes{:});
    msgSupTypes   = sprintf('''%s'', ', supportedTypes{:});
    msg = ['The types ' msgUnSupTypes(1:end-2) ' are unsupported by this function. ', ...
           'This function currently supports the types ', msgSupTypes(1:end-2), '.'];
    error('Collocation:getSimData', msg);
end
        
%% Extract data from X which is needed for multiple types in simVar
% => Avoid computing this multiple times
% -> Extract states
if any(ismember(idxSimTypes, idxTypesNeedState))
    if ~isfield(obj.idx,'states') % check whether states are stored in X
        error('Collocation:getSimData', 'States are not stored in state vector X.')
    end
    % We are not interested in the states at nNodes+1 for
    % getJointMoments(), etc., but we need it to compute q, qd, qdd
    states = X(obj.idx.states); 
end
% -> Extract duration
if any(ismember(idxSimTypes, idxTypesNeedDur))
    if ~isfield(obj.idx,'dur') % check whether duration is stored in X
        warning('Collocation:getSimData', 'Duration is not stored in state vector X.')
        dur = 0; % probably static standing
        h = nan;
    else
        dur = X(obj.idx.dur);
        h = dur/(nNodesDur-1); % time step to compute qdd
    end
end
% -> Extract states derivatives
if any(ismember(idxSimTypes, idxTypesNeedStated))
    % Get xd 
    if ~isnan(h)
        statesd = ( states(:, 2:nNodesDur) - states(:, 1:(nNodesDur-1)) ) / h;
    else
        statesd = zeros(obj.model.nStates, nNodesDur-1);
    end
end
% -> Extract controls
if any(ismember(idxSimTypes, idxTypesNeedControl))
    if ~isfield(obj.idx,'controls') % check whether controls are stored in X
        error('Collocation:getSimData', 'Controls are not stored in state vector X.')
    end
    controls = X(obj.idx.controls(:, 1:nNodes)); % We are not interested in the controls at nNodes+1
end


%% Get simVar
% Add column sim to table if it does not exist already
if ~any(ismember(simVar.Properties.VariableNames, 'sim'))
  simVar.sim = cell(height(simVar), 1);  
end
% Go trough all types which are requested in simTypes
for iType = 1 : length(simTypes) 
    % Get row indices were this type is used in simVar
    idxVar = find(ismember(simVar.type, simTypes{iType}));

    % Get data according to the current type
    switch simTypes{iType}
        case {'translation', 'angle'}
            % Indices in state vector for all simVar names of the current type
            idxSimVar = model.extractState('q',simVar.name(idxVar)); 
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(states(idxSimVar, 1:nNodes)', 1)';
            % Convert from radian to degree if requested unit is degree
            simVar(idxVar, :) = convertUnit(simVar(idxVar, :), simTypes{iType}, 'deg', 'deg', 180/pi);
            
        case 'qdot'
            % Indices in state vector for all simVar names of the current type
            idxSimVar = model.extractState('qdot',simVar.name(idxVar));
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(states(idxSimVar, 1:nNodes)', 1)';
            % Convert from rad/s to deg/s if requested unit is deg/s
            simVar(idxVar, :) = convertUnit(simVar(idxVar, :), simTypes{iType}, 'deg/s', 'deg/s', 180/pi);
            
        case 'moment'
            % Get passive moments and muscle moments
            moments = zeros(model.nDofs,nNodes);
            for iNode = 1:nNodes 
                % Get passive joint moments and muscle moments
                moments(:,iNode) = model.getJointmoments(states(:, iNode), controls(:, iNode));
            end
            % Indices in dofs for all simVar names of the current type
            idxSimVar = model.extractState('q',simVar.name(idxVar));
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(moments(idxSimVar, 1:nNodes)', 1)';
            
        case 'jointPower'
            % Get passive moments and muscle moments
            moments = zeros(model.nDofs,nNodes);
            for iNode = 1:nNodes 
                % Get passive joint moments and muscle moments
                moments(:,iNode) = model.getJointmoments(states(:, iNode), controls(:, iNode));
            end
            % Get the joint velocities
            idxSimVarVel = model.extractState('qdot',simVar.name(idxVar));
            qdot = states(idxSimVarVel, 1:nNodes);
            
            %Get the power
            idxSimVarMom = model.extractState('q',simVar.name(idxVar));
            power = qdot.*moments(idxSimVarMom,:);
            
            simVar.sim(idxVar) = num2cell(power(:, 1:nNodes)', 1)';
            
            %Normalize to bodyweight if requested
            bodyweight = model.bodymass*norm(model.gravity);
            simVar(idxVar, :) = convertUnit(simVar(idxVar, :), simTypes{iType}, 'W/BW', 'W/BW', 1/bodyweight);
            
        case 'GRF'
            % Get GRFs data
            GRF = zeros(12,nNodes);
            for iNode = 1:nNodes 
                GRF(:,iNode) = model.getGRF(states(:, iNode));
            end
            % Put it into simVar for each row
            for iVar = idxVar'
                if isa(model, 'quad_11DOF') || isa(model, 'Quadruped')
                    switch simVar.name{iVar}
                        case 'GRF_x_r_f'
                            simVar.sim{iVar} = GRF(5,:)';
                        case 'GRF_y_r_f'
                            simVar.sim{iVar} = GRF(6,:)';
                        case 'GRF_x_l_f'
                            simVar.sim{iVar} = GRF(7,:)';
                        case 'GRF_y_l_f'
                            simVar.sim{iVar} = GRF(8,:)';
                        case 'GRF_x_r_b'
                            simVar.sim{iVar} = GRF(1,:)';
                        case 'GRF_y_r_b'
                            simVar.sim{iVar} = GRF(2,:)';
                        case 'GRF_x_l_b'
                            simVar.sim{iVar} = GRF(3,:)';
                        case 'GRF_y_l_b'
                            simVar.sim{iVar} = GRF(4,:)';
                    end
                else
                    switch simVar.name{iVar}
                        case 'GRF_x_r'
                            simVar.sim{iVar} = GRF(1,:)';
                        case 'GRF_y_r'
                            simVar.sim{iVar} = GRF(2,:)';
                        case 'GRF_z_r'
                            simVar.sim{iVar} = GRF(3,:)';
                        case 'GRF_x_l'
                            simVar.sim{iVar} = GRF(7,:)';
                        case 'GRF_y_l'
                            simVar.sim{iVar} = GRF(8,:)';
                        case 'GRF_z_l'
                            simVar.sim{iVar} = GRF(9,:)';
                    end
                end
            end
            
        case 'CoP'
            % It is not super efficient to do it separatly from the GRF, but in postprocessing we do not care much.

            % Get GRFs and CoP data
            GRF = zeros(12,nNodes);
            CoP_r = zeros(3,nNodes);
            CoP_l = zeros(3,nNodes);
            for iNode = 1:nNodes
                GRF(:,iNode) = model.getGRF(states(:, iNode));
                [CoP_r(:,iNode), CoP_l(:,iNode)] = model.getCoP(GRF(:, iNode));
            end
            % Put it into simVar for each row
            for iVar = idxVar'
                switch simVar.name{iVar}
                    case 'CoP_x_r'
                        simVar.sim{iVar} = CoP_r(1,:)';
                    case 'CoP_y_r'
                        simVar.sim{iVar} = CoP_r(2,:)';
                    case 'CoP_z_r'
                        simVar.sim{iVar} = CoP_r(3,:)';
                    case 'CoP_x_l'
                        simVar.sim{iVar} = CoP_l(1,:)';
                    case 'CoP_y_l'
                        simVar.sim{iVar} = CoP_l(2,:)';
                    case 'CoP_z_l'
                        simVar.sim{iVar} = CoP_l(3,:)';
                end
            end

        case 'marker'
            % Get simulated marker positions
            marker = zeros(length(idxVar), nNodes);
            qMarker = states(model.extractState('q'), 1:nNodes);
            for iNode = 1 : nNodes
                marker(:, iNode) = model.simuMarker(simVar,qMarker(:,iNode));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(marker(:, 1:nNodes)', 1)';
            % Convert from meter to millimeter if requested unit is mm
            simVar(idxVar, :) = convertUnit(simVar(idxVar, :), simTypes{iType}, 'mm', 'mm', 1000);
            
        case 'acc'
            % Get q, qd, and qdd
            iq = obj.model.extractState('q');
            iqd = obj.model.extractState('qdot');
            q = states(iq, :);
            qd = states(iqd, :);
            qdd = statesd(iqd, :);
            % Get simulated acc signal
            acc = zeros(length(idxVar), nNodesDur-1);
            for iNode = 1 : nNodesDur-1
                acc(:, iNode) = model.simuAccGyro(simVar,q(:,iNode),qd(:,iNode),qdd(:,iNode))';
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(acc(:, 1:nNodesDur-1)', 1)';
            % Convert from m/s^2 to mm/s^2 if requested unit is mm/s^2
            simVar(idxVar, :) = convertUnit(simVar(idxVar, :), simTypes{iType}, 'mm/s^2', 'mm/s^2', 1000);
            
        case 'gyro'
            % Get q, and qd
            iq = obj.model.extractState('q');
            iqd = obj.model.extractState('qdot');
            q = states(iq, :);
            qd = states(iqd, :);
            % Get simulated gyro signal
            gyro = zeros(length(idxVar), nNodes);
            for iNode = 1 : nNodes
                gyro(:, iNode) = model.simuAccGyro(simVar,q(:,iNode),qd(:,iNode))'; % don't need qdd for gyro
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(gyro(:, 1:nNodes)', 1)';
            % Convert from rad/s to deg/s if requested unit is deg/s
            simVar(idxVar, :) = convertUnit(simVar(idxVar, :), simTypes{iType}, 'deg/s', 'deg/s', 180/pi);
            
        case 'a'
            % Get indices of a in state vector
            idxSimVar = model.extractState('a',simVar.name(idxVar));
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(states(idxSimVar, 1:nNodes)', 1)';
            
        case 's'
            % Get indices of a in state vector
            idxSimVar = model.extractState('s',simVar.name(idxVar));
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(states(idxSimVar, 1:nNodes)', 1)';
            
        case 'LMTU'
            % Get length of MTU
            LMTU = zeros(model.nMus,nNodes);
            for iNode = 1:nNodes
                LMTU(:,iNode) = model.getLMTU(states(:,iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(LMTU(idxSimVar, 1:nNodes)', 1)';
            
        case 'LCE'
            % Get length of CE
            LCE = zeros(model.nMus,nNodes);
            for iNode = 1:nNodes
                LCE(:,iNode) = model.getLCE(states(:,iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(LCE(idxSimVar, 1:nNodes)', 1)';
            
        case 'LSEE'
            % Get length of SEE
            LSEE = zeros(model.nMus,nNodes);
            for iNode = 1:nNodes
                LSEE(:,iNode) = model.getLSEE(states(:,iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(LSEE(idxSimVar, 1:nNodes)', 1)';
            
        case 'LdotMTU'
            % Get length change of MTU
            LdotMTU = zeros(model.nMus,nNodesDur-1);
            for iNode = 1:nNodesDur-1
                LdotMTU(:,iNode) = model.getLdotMTU(states(:,iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(LdotMTU(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'LdotCE'
            % Get length change of CE
            LdotCE = zeros(model.nMus,nNodesDur-1);
            for iNode = 1:nNodesDur-1
                LdotCE(:,iNode) = model.getLdotCE(statesd(:, iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(LdotCE(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'LdotSEE'
            % Get length change of SEE
            LdotSEE = zeros(model.nMus,nNodesDur-1);
            for iNode = 1:nNodesDur-1
                LdotSEE(:,iNode) = model.getLdotSEE(states(:,iNode), statesd(:, iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(LdotSEE(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'muscleForce'
            % Get muscle force data
            muscleForce = zeros(model.nMus,nNodes);
            for iNode = 1:nNodes 
                muscleForce(:,iNode) = model.getMuscleforces(states(:,iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(muscleForce(idxSimVar, 1:nNodes)', 1)';
            
        case 'CEForce'
            % Get CE force data
            CEForce = zeros(model.nMus,nNodesDur-1);
            for iNode = 1:nNodesDur-1
                CEForce(:,iNode) = model.getMuscleCEforces(states(:,iNode), statesd(:, iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(CEForce(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'u'
            % Get indices of u in control vector
            idxSimVar = model.extractControl('u',simVar.name(idxVar));
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(controls(idxSimVar, 1:nNodes)', 1)';
            
        case 'torque'
            % Get indices of torque in control vector
            idxSimVar = model.extractControl('torque',simVar.name(idxVar));
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(controls(idxSimVar, 1:nNodes)' * model.mExtraScaleFactor, 1)';
            
        case 'duration'
            simVar.sim(idxVar) = num2cell(dur);
            
        case 'speed'
            if ~isfield(obj.idx,'speed') % check whether speed is stored in X
                error('Collocation:getSimData', 'Speed is not stored in state vector X.')
            end
            simVar.sim(idxVar) = num2cell(X(obj.idx.speed));
            
        case 'footAngle'
            % Get foot angle
            [angle_r, angle_l] = model.getFootAngle(states(:, 1:nNodes));
            % Indices in mus for all simVar names of the current type
            idxR = find(strcmp(simVar.name, 'angle_r'));
            if ~isempty(idxR); simVar.sim{idxR} = angle_r; end
            idxL = find(strcmp(simVar.name, 'angle_l'));
            if ~isempty(idxL); simVar.sim{idxL} = angle_l; end
            
        case 'musclePower'
            % Get power of entire MTU
            musclePower = zeros(model.nMus,nNodes);
            for iNode = 1:nNodes
                musclePower(:,iNode) = model.getMusclePower(states(:,iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(musclePower(idxSimVar, 1:nNodes)', 1)';
            
        case 'CEPower'
            % Get power of contractile element
            CEPower = zeros(model.nMus,nNodesDur-1);
            for iNode = 1:nNodesDur-1 
                CEPower(:,iNode) = model.getMuscleCEpower(states(:,iNode), statesd(:, iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(CEPower(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'SEEPower'
            % Get power of serial elastic element
            SEEPower = zeros(model.nMus,nNodesDur-1);
            for iNode = 1:nNodesDur-1
                SEEPower(:,iNode) = model.getMuscleSEEpower(states(:,iNode), statesd(:, iNode));
            end
            % Indices in mus for all simVar names of the current type
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(model.muscles.Properties.RowNames, simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(SEEPower(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'sdot'
            % Get indices of a in state vector
            idxSimVar = model.extractState('s',simVar.name(idxVar));
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(statesd(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'muscleMetRate'
            % Get indices of u in control vector
            idxSimVar = model.extractControl('u',simVar.name(idxVar));

            %Find t_stim, stimulation at beginning of gait cycle
            t_stim = obj.getStimTime(X);

            bodymass = model.bodymass;          % Bodymass in kg
            
            Edot = zeros(model.nMus,nNodesDur-1);
            for iNode = 1:nNodesDur-1
                curStates = states(:, iNode);
                Edot(:,iNode) = obj.model.getMetabolicRate_pernode(curStates, statesd(:, iNode), controls(:,iNode), t_stim(:,iNode), modelName, 0, 1);
            end
            Edot = Edot/bodymass;
            simVar.sim(idxVar) = num2cell(Edot(idxSimVar, 1:nNodesDur-1)', 1)';
            
        case 'CoM' 
            % Get center of mass
            dim = length(obj.model.getCoM(states(:, 1)));
            CoM = zeros(dim, nNodes);
            for iNode = 1:nNodes
                CoM(:, iNode) = obj.model.getCoM(states(:,iNode));
            end
            % Indices in supportedNames for all simVar names of the current type
            supportedNames = {'CoM_x', 'CoM_y', 'CoM_z'};
            idxSimVar = zeros(size(simVar.name(idxVar)));
            for iName = 1 : size(simVar.name(idxVar), 1)
                idxSimVar(iName) = find(strcmp(supportedNames(1:dim), simVar.name(idxVar(iName))));
            end
            % Put it into simVar for each row
            simVar.sim(idxVar) = num2cell(CoM(idxSimVar, 1:nNodes)', 1)';

    end

end


end

