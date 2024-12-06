%======================================================================
%> @file @Collocation/writeMotionToOsim.m
%> @brief Collocation function to write motion from a solution X to OpenSim files
%> @details
%> Details: Collocation::writeMotionToOsim()
%> 
%> @author Marlies Nitschke
%> @date November, 2018
%======================================================================

%======================================================================
%> @brief Function to write motion from a solution X to OpenSim files
%>
%> @details
%> Saves the motion to the following files:
%> 1. <filename>_kinematics.mot containing the joint angles
%> 2. <filename>_kinematics_subTrans.mot containing the joint angles while
%>    subtracting x and z translation
%> 3. <filename>_kinetics_GRFs.mot containing GRFs
%> 4. <filename>_kinetics_activation.sto containing muscle activation 
%> 5. <filename>_kinetics_moments.sto containing moments
%> 6. <filename>_kinetics_force.sto containing muscle forces
%> 7. <filename>_kinetics_controls.xml containing controls
%> 8. Optional: <filename>_markers.trc containing marker positions of markerTable specified
%> 9. Optional: <filename>_IMUs.csv containing IMU data of varTable (format definition of APDM)
%>
%> No resampling is applied.
%>
%> @param  obj           Collocation class object
%> @param  X             Double matrix: State vector (i.e. result) of the problem
%> @param  filename      String: Filename which is used to save the OpenSim files
%> @param  variableTable (optional) Table: Variables table specifying markers or IMUs with at least the columns:
%>                       type, name, segment, position, direction to call Model.simuMarker or Model.simuAccGyro.
%>                       Can be skipped when empty.
%> @param  rangeNodes    (optional) Array: Specifying first and last node
%>                       which should be exported. This might not work well for symmetric motions
%>                       which are automatically mirrored.
%> @param  t0            (optional) Double: Time at first time point which should be used to start the time vector.
%>                       If this value is not given, the time vector will start at 0.
%======================================================================
function writeMotionToOsim(obj, X, filename, variableTable, rangeNodes, t0)

%% Extract information
% Number of nodes
if nargin > 4
    nNodes = rangeNodes(2) - rangeNodes(1)+1;
    rangeNodes = rangeNodes(1):rangeNodes(2);
else
    nNodes = obj.nNodes;
    rangeNodes = 1 : nNodes;
end
% Start time
if nargin < 6
    t0 = 0;
end

% Model
model  = obj.model;
% Time stamps for each node
if nNodes > 1 % time depending movement
    if ~isfield(obj.idx,'dur')  % check whether speed is stored in X
        error('Model duration is not stored in state vector X and number of nodes is greater than 1.')
    end
    h = X(obj.idx.dur)/(obj.nNodesDur-1);
    duration = h * nNodes;
else % static pose
    h = 1;
end
% States
x = X(obj.idx.states(:, rangeNodes))';
% Controls 
u = X(obj.idx.controls(:, rangeNodes))';

%% Mirror movement if it is symmetric
if obj.isSymmetric && nNodes > 1
    % Get speed of movement
    if ~isfield(obj.idx,'speed')  % check whether speed is stored in X
        error('Model speed is not stored in state vector X.')
    end
    speed = X(obj.idx.speed);
    % Get unit displacement
    unitdisplacement = zeros(model.nStates,1);
    unitdisplacement(model.idxForward) = 1;
    
    % Mirror states
    flip =  x;
    flip = repmat(model.idxSymmetry.xsign,1,nNodes)'.*flip(:,model.idxSymmetry.xindex) + repmat(unitdisplacement,1,nNodes)'*duration*speed;
    x = [x;flip];
    
    % Mirror controls
    flip =  u;
    flip = repmat(model.idxSymmetry.usign,1,nNodes)'.*flip(:,model.idxSymmetry.uindex); % sign is always 1 and there is no displacement
    u = [u;flip];
    
    % Adapt number of nodes and duration
    nNodes = 2* nNodes;
end

%% Make two time points if nNodes == 1
% => OpenSim will not show the position otherwise
if nNodes == 1
    x = [x;x];
    u = [u;u];
    nNodes = 2;
end

%% Get times
times = (0:h:(nNodes-1)*h) +t0;

%% Write joint angles to <filename>_kinematics.mot
% Extract data
names = model.dofs.Properties.RowNames; 
q = x(:, model.extractState('q'));
pelvis_t_xyz = {'pelvis_tx', 'pelvis_ty', 'pelvis_tz'};
idxConvert = ~ismember(names, pelvis_t_xyz); % all which are not translation
q(:, idxConvert) = q(:, idxConvert) / pi * 180;

% Write to file
filenameKinem = [filename '_kinematics.mot'];
inDegrees = 1;
writeMotSto(times, q, names, filenameKinem, inDegrees);


%% Write joint angles to <filename>_kinematics_subTrans.mot
% Subtract translation
idxPelvisX = model.extractState('q', 'pelvis_tx'); % forwards
idxPelvisZ = model.extractState('q', 'pelvis_tz'); % sidewards
xSubTrans = x;
xSubTrans(:, model.idxForward)  = xSubTrans(:, model.idxForward)  - xSubTrans(:,idxPelvisX);
xSubTrans(:, model.idxSideward) = xSubTrans(:, model.idxSideward) - xSubTrans(:,idxPelvisZ);

% Extract data
qSubTrans = xSubTrans(:, model.extractState('q'));
qSubTrans(:, idxConvert) = qSubTrans(:, idxConvert) / pi * 180;
names = model.dofs.Properties.RowNames; 

% Write to file
filenameKinem = [filename '_kinematics_subTrans.mot'];
inDegrees = 1;
writeMotSto(times, qSubTrans, names, filenameKinem, inDegrees);


%% Write GRFs to <filename>_kinetics_GRFs.mot
% Extract data
names = {'ground_force_vx','ground_force_vy','ground_force_vz', ...       % right foot's GRF vector 
         'ground_force_px','ground_force_py','ground_force_pz', ...       % right foot's center of pressure
         'l_ground_force_vx','l_ground_force_vy','l_ground_force_vz', ... % left  foot's GRF vector 
         'l_ground_force_px','l_ground_force_py','l_ground_force_pz', ... % left  foot's center of pressure
         'ground_torque_x','ground_torque_y','ground_torque_z', ...       % right foot's moment vector 
         'l_ground_torque_x','l_ground_torque_y','l_ground_torque_z'};    % left  foot's moment vector 
dataGRF = nan(nNodes, length(names));
convert_BW_to_N = norm(model.gravity) * model.bodymass;  % convert BW to N
for iNode = 1 : nNodes
    [grf, ~] = model.getGRF(x(iNode,:)');
    [CoP_r, CoP_l] = model.getCoP(grf);   
    grf_SI_unit = convert_BW_to_N * grf;
    
    dataGRF(iNode, 1:3)   = grf_SI_unit(1:3)';             % right foot's GRF vector 
    dataGRF(iNode, 4:6)   = [CoP_r(1) CoP_r(2) CoP_r(3)];  % right foot's center of pressure 
    dataGRF(iNode, 7:9)   = grf_SI_unit(7:9)';             % left foot's GRF vector 
    dataGRF(iNode, 10:12) = [CoP_l(1) CoP_l(2) CoP_l(3)];  % left foot's center of pressure 
    dataGRF(iNode, 13:15) = grf_SI_unit(4:6)';            % right foot's moment vector
    dataGRF(iNode, 16:18) = grf_SI_unit(10:12)';          % left foot's moment vector
   
end 
     
% Write to file
filenameGRF = [filename '_kinetics_GRFs.mot'];
inDegrees = 0;
writeMotSto(times, dataGRF, names, filenameGRF, inDegrees);

%% Write GRFs to <filename>_kinetics_GRFs_subTrans.mot
for iNode = 1 : nNodes
    [grf, ~] = model.getGRF(xSubTrans(iNode,:)');
    [CoP_r, CoP_l] = model.getCoP(grf);   
    grf_SI_unit = convert_BW_to_N * grf;
    
    dataGRF(iNode, 1:3)   = grf_SI_unit(1:3)';             % right foot's GRF vector 
    dataGRF(iNode, 4:6)   = [CoP_r(1) CoP_r(2) CoP_r(3)];  % right foot's center of pressure 
    dataGRF(iNode, 7:9)   = grf_SI_unit(7:9)';             % left foot's GRF vector 
    dataGRF(iNode, 10:12) = [CoP_l(1) CoP_l(2) CoP_l(3)];  % left foot's center of pressure 
    dataGRF(iNode, 13:15) = grf_SI_unit(4:6)';            % right foot's moment vector
    dataGRF(iNode, 16:18) = grf_SI_unit(10:12)';          % left foot's moment vector
   
end 
     
% Write to file
filenameGRF = [filename '_kinetics_GRFs_subTrans.mot'];
inDegrees = 0;
writeMotSto(times, dataGRF, names, filenameGRF, inDegrees);

%% Write activation to <filename>_kinetics_activation.sto
% Extract data
a = x(:, model.extractState('a'));
names = model.muscles.Properties.RowNames; 

% Write to file
filenamea = [filename '_kinetics_activation.sto'];
inDegrees = 0;
writeMotSto(times, a, names, filenamea, inDegrees);


%% Write moments to <filename>_kinetics_moments.sto
% Extract data
names = cell(model.nDofs, 1);
for iDof = 1 : model.nDofs
    if ismember(model.dofs.Properties.RowNames{iDof}, pelvis_t_xyz)
        names{iDof} = sprintf('\t%s_force', model.dofs.Properties.RowNames{iDof});   
    else
        names{iDof} = sprintf('\t%s_moment', model.dofs.Properties.RowNames{iDof});
    end
end
M = nan(nNodes, model.nDofs);
for iNode = 1 : nNodes
    M(iNode, :) = model.getJointmoments(x(iNode,:)', u(iNode,:)');
end

% Write to file
filenameMom = [filename '_kinetics_moments.sto'];
inDegrees = 0;
writeMotSto(times, M, names, filenameMom, inDegrees);


%% Write muscle forces to <filename>_kinetics_force.sto
% Extract data
names = {model.muscles.Properties.RowNames{:}, ...                                                   % muscle forces
         'calcn_r_ExternalForce_1_Fx','calcn_r_ExternalForce_1_Fy','calcn_r_ExternalForce_1_Fz', ... % right foot's GRF vector 
         'calcn_r_ExternalForce_1_px','calcn_r_ExternalForce_1_py','calcn_r_ExternalForce_1_pz', ... % right foot's center of pressure
         'calcn_r_ExternalForce_1_Tx','calcn_r_ExternalForce_1_Ty','calcn_r_ExternalForce_1_Tz', ... % right foot's moment vector 
         'calcn_l_ExternalForce_2_Fx','calcn_l_ExternalForce_2_Fy','calcn_l_ExternalForce_2_Fz', ... % left foot's GRF vector 
         'calcn_l_ExternalForce_2_px','calcn_l_ExternalForce_2_py','calcn_l_ExternalForce_2_pz', ... % left foot's center of pressure
         'calcn_l_ExternalForce_2_Tx','calcn_l_ExternalForce_2_Ty','calcn_l_ExternalForce_2_Tz'};    % left foot's moment vector 
dataF = nan(nNodes, length(names));
for iNode = 1 : nNodes
    forces = model.getMuscleforces(x(iNode,:)');    % getting muscle forces
    [grf, ~] = model.getGRF(x(iNode,:)');           % getting GRF and ground reaction moment
    [CoP_r, CoP_l] = model.getCoP(grf);             % getting center of presusre
    grf_SI_unit = convert_BW_to_N * grf;            % convert GRFs to N
    
    dataF(iNode, 1:model.nMus)       = forces;                        % muscle forces
    dataF(iNode, model.nMus+(1:3))   = grf_SI_unit(1:3)';             % right foot's GRF vector 
    dataF(iNode, model.nMus+(4:6))   = [CoP_r(1) CoP_r(2) CoP_r(3)];  % right foot's center of pressure 
    dataF(iNode, model.nMus+(7:9))   = grf_SI_unit(7:9)';             % left foot's GRF vector 
    dataF(iNode, model.nMus+(10:12)) = [CoP_l(1) CoP_l(2) CoP_l(3)];  % left foot's center of pressure 
    dataF(iNode, model.nMus+(13:15)) = grf_SI_unit(4:6)';             % right foot's moment vector
    dataF(iNode, model.nMus+(16:18)) = grf_SI_unit(10:12)';           % left foot's moment vector
end 
     
% Write to file
filenameF = [filename '_kinetics_force.sto'];
inDegrees = 0;
writeMotSto(times, dataF, names, filenameF, inDegrees);


%% Write controls to <filename>_kinetics_controls.xml
% Extract data
a = u(:, model.extractControl('u'));
names = model.muscles.Properties.RowNames; 

% Write to file
filenameu = [filename '_kinetics_controls.xml'];
writeControlXML(times, a, names, filenameu);


%% Optional: Write marker positions to <filename>_markers.trc
if nargin > 3 && ~isempty(variableTable) && ~isempty(variableTable(strcmp(variableTable.type, 'marker'), :))
    % Get only the rows with marker data
    variableTableMarker = variableTable(strcmp(variableTable.type, 'marker'), :);
    markerNames = unique(variableTableMarker.name,'stable');

    % Simulate marker data
    q = x(:, model.extractState('q'));
    markerData = nan(nNodes, height(variableTableMarker));
    for iNode = 1 : nNodes
        markerData(iNode, :) = model.simuMarker(variableTableMarker, q(iNode, :)')*1000; % convert to mm
    end

    % Write to file
    filenameMarker = [filename '_markers.trc'];
    inMM = 1;
    writeMarkerTrc(times, markerData, markerNames, filenameMarker, inMM);

end


%% Optional: Write IMU data to <filename>_IMUs.csv
if nargin > 3 && ~isempty(variableTable) && ...
        (~isempty(variableTable(strcmp(variableTable.type, 'acc'), :)) || ~isempty(variableTable(strcmp(variableTable.type, 'gyro'), :)))
    % Get only the rows with IMU data
    variableTableIMU = variableTable(ismember(variableTable.type, {'acc', 'gyro'}), :);

    % Simulate IMU data
    iq = obj.model.extractState('q');
    iqd = obj.model.extractState('qdot');
    simVarAll = zeros(nNodes,height(variableTableIMU));
    for iNode = 1:nNodes
        % determine the q and qd and qdd
        q = x(iNode, iq)';
        qd = x(iNode, iqd)';
        if iNode < nNodes
            qdd = (x(iNode+1, iqd) - x(iNode, iqd) )' / h;
        else
            qdd = NaN(size(qd)); % Acc can't be computed for last node and will be nan
        end
        
        % simulate gyroscope signals
        simVarAll(iNode, :) = obj.model.simuAccGyro(variableTableIMU, q, qd, qdd);
    end
    
    % Add data to variables table
    variableTableIMU = removevars(variableTableIMU, {'mean', 'var'});
    variableTableIMU.mean = mat2cell(simVarAll, size(simVarAll, 1), ones(size(simVarAll, 2), 1))';

    % Write to file
    filenameIMU = [filename '_IMUs.csv'];
    writeIMUcsv(times, variableTableIMU, filenameIMU);

end

end



