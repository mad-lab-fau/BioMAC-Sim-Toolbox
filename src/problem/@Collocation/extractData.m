%======================================================================
%> @file Collocation/extractData.m
%> @brief Collocation function to extract simulated and tracked data of the problem
%> @details
%> Details: Collocation::extractData()
%>
%> @author Marlies Nitschke
%> @date July, 2021
%======================================================================

%======================================================================
%> @brief Function to extract simulated and tracked data of the problem
%>
%> @details
%> Using the fields of settings, you can specify which simulated data will be
%> extracted. For each type, you define the names. The fieldname corresponds
%> to the type of the variable. You do not have to use all of the fields!
%> Just use the ones you need.
%>
%> Example entries for all possible fields of settings:
%> @code
%> settings.translation   = {'pelvis_tx', 'pelvis_ty'};
%> settings.angle         = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
%> settings.qdot          = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
%> settings.moment        = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
%> settings.torque        = {'arm_flex_r', 'arm_add_r', 'arm_rot_r', 'elbow_flex_r', 'pro_sup_r'};
%> settings.GRF           = {'GRF_x_r', 'GRF_y_r'};
%> settings.CoP           = {'CoP_x_r'};
%> settings.acc           = accTable; % variables table specifying acc with at least the columns: type, name, unit, segment, position, direction
%> settings.gyro          = gyroTable; % variables table specifying gyro with at least the columns: type, name, unit, segment, position, direction
%> settings.marker        = markerTable; % variables table specifying markers with at least the columns: type, name, unit, segment, position, direction
%> settings.u             = {'Iliopsoas_r', 'Glutei_r', 'Hamstrings_r', 'Rectus_r', 'Vasti_r', 'Gastroc_r', 'Soleus_r', 'TibialisAnt_r'};
%> settings.a             = settings.u;
%> settings.s             = settings.u;
%> settings.sdot          = settings.u;
%> settings.LMTU          = settings.u;
%> settings.LCE           = settings.u;
%> settings.LSEE          = settings.u;
%> settings.LdotMTU       = settings.u;
%> settings.LdotCE        = settings.u;
%> settings.LdotSEEE      = settings.u;
%> settings.muscleForce   = settings.u;
%> settings.CEForce       = settings.u;
%> settings.musclePower   = settings.u;
%> settings.CEPower       = settings.u;
%> settings.SEEPower      = settings.u;
%> settings.muscleMetRate = settings.u;
%> settings.jointPower    = {'hip_flexion_r','knee_angle_r', 'ankle_angle_r'};
%> settings.CoM           = {'CoM_x', 'CoM_y'};
%> settings.footAngle     = {'angle_r', 'angle_l' };
%> settings.standing      = {'standing_r', 'standing_l'};
%> settings.duration      = 1;
%> @endcode
%>
%>
%> By default, the following variables will automatically be extracted if
%> you do not overwrite the respective fields in settings:
%> - translation: all
%> - angle: all joint angles of both sides (no global orientation)
%> - moment: all joint moments of both sides
%> - GRF: all of both sides
%> - a and u:
%>   - 2D: all of both sides
%>   - 3D: in totla 16 largest muscles of both sides
%>
%>
%> The function will further extract all tracked variables which are defined in
%> settings.
%>
%> Additionally, you can use the input variableTable to specify a variable table 
%> with additional (reference) data which was not tracked. You can for example 
%> use the variables table in a TrackingData object:
%> @code
%> variableTable = trackingdata.variables;
%> @endcode
%>
%>
%> @param  obj            Collocation class object
%> @param  X              Double matrix: State vector (i.e. result) of the problem
%> @param  settings       (optional) Struct: Describing which data should be extracted (see details above; use empty to skip)
%> @param  variableTable  (optional) Table: Additional reference data with the data beeing in the 
%>                        columns 'mean' and 'var'. (default: empty)
%> @param  getFullCycle   (optional) Boolean: Defines if a full gait cycle should be constructed
%>                        for symmetric simulations:
%>                        - 0: requested variable names for the half cycle
%>                        - 1: requested variable names for the entire cycle by using symmetry
%>
%> @retval simVarTable    Table: Summarizing all simulated, tracked and optionally additional data 
%>                        which was requested. The three different data kinds are saved in the following columns:
%>                        - simulated data: 'sim'
%>                        - tracked data: 'mean' and 'var'
%>                        - additional data: 'mean_extra' and 'var_extra'
%======================================================================
function simVarTable = extractData(obj, X, settings, variableTable, getFullCycle)

%% Read general information and do intialization
% Get properties
model = obj.model;
objTermsNames = {obj.objectiveTerms.name};

% Get default for settings
if isa(model, 'quad_11DOF') || isa(model, 'Quadruped')
    translationNames = {'trunk_tx', 'trunk_ty'};
    orientationNames = {'trunk_q'};
else
    translationNames = {'pelvis_tx', 'pelvis_ty', 'pelvis_tz'};
    orientationNames = {'pelvis_tilt', 'pelvis_list','pelvis_rotation','pelvis_obliquity'};
end
if nargin < 3 || isempty(settings)
    settings = struct();
end
if ~isfield(settings, 'translation')
    % default: all translations
    dofNames = model.dofs.Properties.RowNames;
    settings.translation = dofNames(ismember(dofNames, translationNames));
end
if ~isfield(settings, 'angle')
    % default: all joint angles of both sides (no global orientation)
    dofNames = model.dofs.Properties.RowNames;
    settings.angle = dofNames(~ismember(dofNames, [translationNames, orientationNames]));
end
if ~isfield(settings, 'moment')
    % default: all joint moments of both sides
    dofNames = model.dofs.Properties.RowNames;
    settings.moment = dofNames(~ismember(dofNames, [translationNames, orientationNames]));
end
if ~isfield(settings, 'GRF')
    % default: all GRFs of both sides
    if isa(model, 'Gait2dc') || isa(model, 'Gait2d_osim')
        settings.GRF = {'GRF_x_r', 'GRF_y_r', 'GRF_x_l', 'GRF_y_l'};
    elseif isa(model, 'Gait3d')
        settings.GRF = {'GRF_x_r', 'GRF_y_r', 'GRF_z_r', 'GRF_x_l', 'GRF_y_l', 'GRF_z_l'};
    elseif isa(model, 'quad_11DOF') || isa(model,'Quadruped')
        settings.GRF = {'GRF_x_r_f', 'GRF_y_r_f', 'GRF_x_l_f', 'GRF_y_l_f', 'GRF_x_r_b', 'GRF_y_r_b', 'GRF_x_l_b', 'GRF_y_l_b'};
    else
        error('Model is unknown');
    end
end
typesMusc = {'a', 'u'};
for iType = 1 : numel(typesMusc)
    type = typesMusc{iType};
    if ~isfield(settings, type)
        % default: all for 2D and 16 largest for 3D of both sides
        if isa(model, 'Gait2dc') || isa(model, 'Gait2d_osim')
            settings.(type) = model.muscles.Properties.RowNames;
        elseif isa(model, 'Gait3d')
            muscleNames = model.muscles.Properties.RowNames;
            nMuscles = 16;
            [~, idxSorting] = sort(model.muscles.weight, 'descend');
            settings.(type) = muscleNames(idxSorting(1:nMuscles));
        elseif isa(model, 'quad_11DOF') || isa(model,'Quadruped')
            settings.(type) = [];
        else
            error('Model is unknown');
        end
    end
end

% Get other defaults
if nargin < 4
   variableTable = []; 
end

if obj.isSymmetric && nargin > 4 && getFullCycle
    doSymmetry = 1;
else
    doSymmetry = 0;
end
        
% Initialize a table containing all data
simVarTable = table();

%% Get translations
typeStr = 'translation';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    unitStr = 'm';

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxTrans = model.extractState('q', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxTrans));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxTrans));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
        % Unit displacement for translation
        unitdisplacement = zeros(1, numel(idxTrans));
        unitdisplacement(ismember(idxTrans, intersect(model.idxForward, idxTrans))) = 1;
        displacement = unitdisplacement * X(obj.idx.speed) * X(obj.idx.dur);
        displacement = [zeros(1, numel(idxTrans)); displacement];
        displacement = num2cell(displacement(:)); % concat it alternating
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({unitStr}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked
    simVarTableType.mean = cell(height(simVarTableType), 1);  
    simVarTableType.var = cell(height(simVarTableType), 1);  
    idxObjTermTrans = find(strcmp(objTermsNames, 'trackTranslations')); % tracked in separate term
    idxObjTermTransAng = find(strcmp(objTermsNames, 'trackTranslationsAndAngles')); % tracked with angles
    idxObjTerm = [idxObjTermTrans, idxObjTermTransAng]; % If not empty, translations were tracked
    if length(idxObjTerm) == 1
        variables = obj.objectiveTerms(idxObjTerm).varargin{:}.variables;
        variables = convertUnit(variables, typeStr, 'mm', unitStr, 0.001); % convert from mm to m if needed
        for iTrans = 1 : height(simVarTableType)
           idxVar = find(strcmp(variables.type, typeStr) & strcmp(variables.name, simVarTableType.name{iTrans}));
           if ~isempty(idxVar)
               assert(strcmp(variables.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variables.unit{idxVar}, typeStr);
               simVarTableType.mean{iTrans} = variables.mean{idxVar};
               simVarTableType.var{iTrans}  = variables.var{idxVar};
           end
        end
    elseif length(idxObjTerm) > 1
        error('Translations were tracked in multiple objective terms. The function extractData() can not deal with this.')        
    end 
    
    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra = cell(height(simVarTableType), 1);  
    if ~isempty(variableTable)
        variableTable = convertUnit(variableTable, typeStr, 'mm', unitStr, 0.001); % convert from mm to m if needed
        for iAng = 1 : height(simVarTableType)
            idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iAng}));
            if ~isempty(idxVar)
                assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
                simVarTableType.mean_extra{iAng} = variableTable.mean{idxVar};
                simVarTableType.var_extra{iAng}  = variableTable.var{idxVar};
            end
        end
    end
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym, displacement);
    end

    % Add simVarTableType to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get angles
typeStr = 'angle';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    unitStr = 'deg';

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxAngle = model.extractState('q', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxAngle));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left
        signSym1 = num2cell(model.idxSymmetry.xsign(idxAngle));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({unitStr}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked
    simVarTableType.mean = cell(height(simVarTableType), 1);  
    simVarTableType.var = cell(height(simVarTableType), 1);   
    idxObjTermAngles = find(strcmp(objTermsNames, 'trackAngles')); % tracked in separate term
    idxObjTermTransAng = find(strcmp(objTermsNames, 'trackTranslationsAndAngles')); % tracked with translations
    idxObjTerm = [idxObjTermAngles, idxObjTermTransAng]; % If not empty, angles were tracked
    if length(idxObjTerm) == 1
        variables = obj.objectiveTerms(idxObjTerm).varargin{:}.variables;
        variables = convertUnit(variables, typeStr, 'rad', unitStr, 180/pi); % convert from rad to deg if needed
        for iAng = 1 : height(simVarTableType)
           idxVar = find(strcmp(variables.type, typeStr) & strcmp(variables.name, simVarTableType.name{iAng}));
           if ~isempty(idxVar)
               assert(strcmp(variables.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variables.unit{idxVar}, typeStr);
               simVarTableType.mean{iAng} = variables.mean{idxVar};
               simVarTableType.var{iAng}  = variables.var{idxVar};
           end
        end
    elseif length(idxObjTerm) > 1
        error('Angles were tracked in multiple objective terms. The function extractData() can not deal with this.')   
    end
    
    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra = cell(height(simVarTableType), 1); 
    if ~isempty(variableTable)
        variableTable = convertUnit(variableTable, typeStr, 'rad', unitStr, 180/pi); % convert from rad to deg if needed
        for iAng = 1 : height(simVarTableType)
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iAng}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{iAng} = variableTable.mean{idxVar};
               simVarTableType.var_extra{iAng}  = variableTable.var{idxVar};
           end
        end
    end 
     
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get qdots
typeStr = 'qdot';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxAngle = model.extractState('qdot', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxAngle));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left
        signSym1 = num2cell(model.idxSymmetry.xsign(idxAngle));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'deg/s'}, numel(name), 1);
    unit(ismember(name, translationNames)) = {'m/s'}; % Translation has different unit
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for qdot!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for qdot
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get moments
typeStr = 'moment';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    unitStr = 'Nm';

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxDof = model.extractState('q', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxDof));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxDof));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({unitStr}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for moments!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);

    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        for iMom = 1 : height(simVarTableType)
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iMom}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{iMom} =  variableTable.mean{idxVar};
               simVarTableType.var_extra{iMom}  =  variableTable.var{idxVar};
           end
        end
    end
         
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get torques
typeStr = 'torque';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    unitStr = 'Nm';

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxTor = model.extractControl('torque', settings.(typeStr));
        namesSym = model.controls.name(model.idxSymmetry.uindex(idxTor));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.usign(idxTor));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'Nm'}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for torques!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);

    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        for iMom = 1 : height(simVarTableType)
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iMom}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{iMom} =  variableTable.mean{idxVar};
               simVarTableType.var_extra{iMom}  =  variableTable.var{idxVar};
           end
        end
    end

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get GRFs
typeStr = 'GRF';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    unitStr = 'BW';
    bodyweight = model.bodymass*norm(model.gravity);

    % Create table to get simulated data
    if doSymmetry
        % Hard code the indices to get names which are symmetric to the
        % requested ones
        namesAll = {'GRF_x_r', 'GRF_y_r', 'GRF_z_r', 'GRF_x_l', 'GRF_y_l', 'GRF_z_l'}; % See todo in header
        isRequested = ismember(namesAll, settings.(typeStr));
        idxSymGRF = [4, 5, 6, 1, 2, 3];
        namesSym = namesAll(idxSymGRF(isRequested));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSymGRF = [1, 1, -1, 1, 1, -1];
        signSym1 = num2cell(signSymGRF(isRequested));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({unitStr}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var = cell(height(simVarTableType), 1);
    idxObjTerm = find(strcmp(objTermsNames, 'trackGRF')); 
    if length(idxObjTerm) == 1
        variables = obj.objectiveTerms(idxObjTerm).varargin{:}.variables;
        variables = convertUnit(variables, typeStr, 'N', unitStr, 1/bodyweight); % convert from N to BW if needed
        variables = convertUnit(variables, typeStr, 'BW%', unitStr, 0.01); % convert from BW% to BW if needed
        for iGRF = 1 : height(simVarTableType)
           idxVar = find(strcmp(variables.type, typeStr) & strcmp(variables.name, simVarTableType.name{iGRF}));
           if ~isempty(idxVar)
               assert(strcmp(variables.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variables.unit{idxVar}, typeStr);
               simVarTableType.mean{iGRF} = variables.mean{idxVar};
               simVarTableType.var{iGRF}  = variables.var{idxVar};
           end
        end    
    end  
    
    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        variableTable = convertUnit(variableTable, typeStr, 'N', unitStr, 1/bodyweight); % convert from N to BW if needed
        variableTable = convertUnit(variableTable, typeStr, 'BW%', unitStr, 0.01); % convert from BW% to BW if needed
        for iGRF = 1 : height(simVarTableType)
            idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iGRF}));
            if ~isempty(idxVar)
                assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
                simVarTableType.mean_extra{iGRF} =  variableTable.mean{idxVar};
                simVarTableType.var_extra{iGRF}  =  variableTable.var{idxVar};
            end
        end
    end

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get CoPs
typeStr = 'CoP';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    unitStr = 'm';

    % Create table to get simulated data
    if doSymmetry
        % Hard code the indices to get names which are symmetric to the
        % requested ones
        namesAll = {'CoP_x_r', 'CoP_y_r', 'CoP_z_r', 'CoP_x_l', 'CoP_y_l', 'CoP_z_l'}; % See todo in header
        isRequested = ismember(namesAll, settings.(typeStr));
        idxSymCoP = [4, 5, 6, 1, 2, 3];
        namesSym = namesAll(idxSymCoP(isRequested));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSymCoP = [1, 1, -1, 1, 1, -1];
        signSym1 = num2cell(signSymCoP(isRequested));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({unitStr}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for center of pressure!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);

    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        variableTable = convertUnit(variableTable, typeStr, 'mm', unitStr, 0.001); % convert from mm to m if needed
        for iCoP = 1 : height(simVarTableType)
            idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iCoP}));
            if ~isempty(idxVar)
                assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
                simVarTableType.mean_extra{iCoP} =  variableTable.mean{idxVar};
                simVarTableType.var_extra{iCoP}  =  variableTable.var{idxVar};
            end
        end
    end

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get accleration data
typeStr = 'acc';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr))

    % Get table to get simulated data and remove other data in the table
    idxCol = ismember(settings.(typeStr).Properties.VariableNames, {'type', 'name', 'unit', 'segment', 'position', 'direction'});
    simVarTableType = settings.(typeStr)(:, idxCol);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    assert(numel(unique(simVarTableType.unit)) == 1, 'Acc table should only contain one unit');
    unitStr = simVarTableType.unit{1}; % Take first

    
    % Get tracking data (mean and var) if it was tracked
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    idxObjTerm = find(strcmp(objTermsNames, 'trackAcc'));
    if length(idxObjTerm) == 1
        variables = obj.objectiveTerms(idxObjTerm).varargin{1}.variables;
        if strcmp(unitStr, 'm/s^2')
            variables = convertUnit(variables, typeStr, 'mm/s^2', unitStr, 0.001); % convert from mm/s^2 to m/s^2 if needed
        elseif strcmp(unitStr, 'mm/^2')
            variables = convertUnit(variables, typeStr, 'm/^2', unitStr, 1000); % convert from m/s^2 to mm/s^2 if needed
        end
        for iAcc = 1 : height(simVarTableType)
            nDim = length(simVarTableType.direction(iAcc,:));
            idxVar = find(strcmp(variables.type, typeStr) & ...
                          strcmp(variables.name, simVarTableType.name{iAcc}) & ...
                          strcmp(variables.segment, simVarTableType.segment{iAcc}) & ...
                          all(variables.direction(:, 1:nDim)==simVarTableType.direction(iAcc,:), 2) & ...
                          all(variables.position(:, 1:nDim)==simVarTableType.position(iAcc,:), 2));
            if ~isempty(idxVar)
                assert(strcmp(variables.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variables.unit{idxVar}, typeStr);
                simVarTableType.mean{iAcc} = variables.mean{idxVar};
                simVarTableType.var{iAcc}  = variables.var{idxVar};
            end
        end
    end

    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        if strcmp(unitStr, 'm/s^2')
            variableTable = convertUnit(variableTable, typeStr, 'mm/s^2', unitStr, 0.001); % convert from mm/s^2 to m/s^2 if needed
        elseif strcmp(unitStr, 'mm/^2')
            variableTable = convertUnit(variableTable, typeStr, 'm/^2', unitStr, 1000); % convert from m/s^2 to mm/s^2 if needed
        end
        for iAcc = 1 : height(simVarTableType)
            nDim = length(simVarTableType.direction(iAcc,:));
            idxVar = find(strcmp(variableTable.type, typeStr) & ...
                          strcmp(variableTable.name, simVarTableType.name{iAcc}) & ...
                          strcmp(variableTable.segment, simVarTableType.segment{iAcc}) & ...
                          all(variableTable.direction(:, 1:nDim)==simVarTableType.direction(iAcc,:), 2) & ...
                          all(variableTable.position(:, 1:nDim)==simVarTableType.position(iAcc,:), 2));
            if ~isempty(idxVar)
                assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
                simVarTableType.mean_extra{iAcc} =  variableTable.mean{idxVar};
                simVarTableType.var_extra{iAcc}  =  variableTable.var{idxVar};
            end
        end
    end

    % Get NaNs for second half of cycle since we can not assume symmetric placement
    if doSymmetry
       simVarTableType = addNaNsForSecondHalf(simVarTableType);
    end

    % Add simVarTable to table containing all data
    simVarTable = [addDirPosSegCol(simVarTableType, simVarTable); simVarTableType];

end

%% Get gyroscope data
typeStr = 'gyro';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr))

    % Get table to get simulated data and remove other data in the table
    idxCol = ismember(settings.(typeStr).Properties.VariableNames, {'type', 'name', 'unit', 'segment', 'position', 'direction'});
    simVarTableType = settings.(typeStr)(:, idxCol);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    assert(numel(unique(simVarTableType.unit)) == 1, 'Gyro table should only contain one unit');
    unitStr = simVarTableType.unit{1}; % Take first
    
    % Get tracking data (mean and var) if it was tracked
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    idxObjTerm = find(strcmp(objTermsNames, 'trackGyro'));
    if length(idxObjTerm) == 1
        variables = obj.objectiveTerms(idxObjTerm).varargin{1}.variables;
        if strcmp(unitStr, 'rad/s')
            variables = convertUnit(variables, typeStr, 'deg/s', unitStr, pi/180); % convert from deg/s to rad/s if needed
        elseif strcmp(unitStr, 'deg/s')
            variables = convertUnit(variables, typeStr, 'rad/s', unitStr, 180/pi); % convert from rad/s to deg/s if needed
        end
        for iGyro = 1 : height(simVarTableType)
            nDim = length(simVarTableType.direction(iGyro,:));
            idxVar = find(strcmp(variables.type, typeStr) & ...
                          strcmp(variables.name, simVarTableType.name{iGyro}) & ...
                          strcmp(variables.segment, simVarTableType.segment{iGyro}) & ...
                          all(variables.direction(:, 1:nDim)==simVarTableType.direction(iGyro,:), 2) & ...
                          all(variables.position(:, 1:nDim)==simVarTableType.position(iGyro,:), 2));
            if ~isempty(idxVar)
                assert(strcmp(variables.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variables.unit{idxVar}, typeStr);
                simVarTableType.mean{iGyro} = variables.mean{idxVar};
                simVarTableType.var{iGyro}  = variables.var{idxVar};
            end
        end
    end

    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        if strcmp(unitStr, 'rad/s')
            variableTable = convertUnit(variableTable, typeStr, 'deg/s', unitStr, pi/180); % convert from deg/s to rad/s if needed
        elseif strcmp(unitStr, 'deg/s')
            variableTable = convertUnit(variableTable, typeStr, 'rad/s', unitStr, 180/pi); % convert from rad/s to deg/s if needed
        end
        for iGyro = 1 : height(simVarTableType)
            nDim = length(simVarTableType.direction(iGyro,:));
            idxVar = find(strcmp(variableTable.type, typeStr) & ...
                          strcmp(variableTable.name, simVarTableType.name{iGyro}) & ...
                          strcmp(variableTable.segment, simVarTableType.segment{iGyro}) & ...
                          all(variableTable.direction(:, 1:nDim)==simVarTableType.direction(iGyro,:), 2) & ...
                          all(variableTable.position(:, 1:nDim)==simVarTableType.position(iGyro,:), 2));
            if ~isempty(idxVar)
                assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
                simVarTableType.mean_extra{iGyro} =  variableTable.mean{idxVar};
                simVarTableType.var_extra{iGyro}  =  variableTable.var{idxVar};
            end
        end
    end

    % Get NaNs for second half of cycle since we can not assume symmetric placement
    if doSymmetry
       simVarTableType = addNaNsForSecondHalf(simVarTableType);
    end

    % NaNs have to be added for direction and position
    if ~(isfield(settings, 'acc') && ~isempty(settings.('acc'))) % Acc was not plotted
        simVarTable = addDirPosSegCol(simVarTableType, simVarTable);
    end
    
    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];
    
end

%% Get marker data (only working for tracked data)
typeStr = 'marker';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr))

    % Get table to get simulated data and remove other data in the table
    idxCol = ismember(settings.(typeStr).Properties.VariableNames, {'type', 'name', 'unit', 'segment', 'position', 'direction'});
    simVarTableType = settings.(typeStr)(:, idxCol);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    assert(numel(unique(simVarTableType.unit)) == 1, 'Marker table should only contain one unit');
    unitStr = simVarTableType.unit{1}; % Take first

    % Get tracking data if it was tracked
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    idxObjTerm = find(strcmp(objTermsNames, 'trackMarker'));
    if length(idxObjTerm) == 1
        variables = obj.objectiveTerms(idxObjTerm).varargin{:}.variables;
        if strcmp(unitStr, 'm')
            variables = convertUnit(variables, typeStr, 'mm', unitStr, 0.001); % convert from mm to m if needed
        elseif strcmp(unitStr, 'mm')
            variables = convertUnit(variables, typeStr, 'm', unitStr, 1000); % convert from m to mm if needed
        end
        for iMarker = 1 : height(simVarTableType)
            nDim = length(simVarTableType.direction(iMarker,:));
            idxVar = find(strcmp(variables.type, typeStr) & ...
                          strcmp(variables.name, simVarTableType.name{iMarker}) & ...
                          strcmp(variables.segment, simVarTableType.segment{iMarker}) & ...
                          all(variables.direction(:, 1:nDim)==simVarTableType.direction(iMarker,:), 2) & ...
                          all(variables.position(:, 1:nDim)==simVarTableType.position(iMarker,:), 2));
            if ~isempty(idxVar)
                assert(strcmp(variables.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variables.unit{idxVar}, typeStr);
                simVarTableType.mean{iMarker} = variables.mean{idxVar};
                simVarTableType.var{iMarker}  = variables.var{idxVar};
            end
        end
    end
    
    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        if strcmp(unitStr, 'm')
            variableTable = convertUnit(variableTable, typeStr, 'mm', unitStr, 0.001); % convert from mm to m if needed
        elseif strcmp(unitStr, 'mm')
            variableTable = convertUnit(variableTable, typeStr, 'm', unitStr, 1000); % convert from m to mm if needed
        end
        for iMarker = 1 : height(simVarTableType)
            nDim = length(simVarTableType.direction(iMarker,:));
            idxVar = find(strcmp(variableTable.type, typeStr) & ...
                          strcmp(variableTable.name, simVarTableType.name{iMarker}) & ...
                          strcmp(variableTable.segment, simVarTableType.segment{iMarker}) & ...
                          all(variableTable.direction(:, 1:nDim)==simVarTableType.direction(iMarker,:), 2) & ...
                          all(variableTable.position(:, 1:nDim)==simVarTableType.position(iMarker,:), 2));
            if ~isempty(idxVar)
                assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
                simVarTableType.mean_extra{iMarker} =  variableTable.mean{idxVar};
                simVarTableType.var_extra{iMarker}  =  variableTable.var{idxVar};
            end
        end
    end

    % Get NaNs for second half of cycle since we can not assume symmetric placement
    if doSymmetry
       simVarTableType = addNaNsForSecondHalf(simVarTableType);
    end

    % NaNs have to be added for direction and position
    if (~(isfield(settings, 'acc') && ~isempty(settings.('acc')))) || ...
       (~(isfield(settings, 'gyro') && ~isempty(settings.('gyro')))) % Acc and gyro were not requested
        simVarTable = addDirPosSegCol(simVarTableType, simVarTable);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; simVarTableType];

end

%% Get u
typeStr = 'u';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxu = model.extractControl('u', settings.(typeStr));
        namesSym = model.controls.name(model.idxSymmetry.uindex(idxu));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.usign(idxu));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({''}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for u!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        for iu = 1 : height(simVarTableType)
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iu}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, ''), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{iu} =  variableTable.mean{idxVar};
               simVarTableType.var_extra{iu}  =  variableTable.var{idxVar};
           end
        end
    end
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get a
typeStr = 'a';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({''}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for a!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        for ia = 1 : height(simVarTableType)
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{ia}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, ''), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{ia} =  variableTable.mean{idxVar};
               simVarTableType.var_extra{ia}  =  variableTable.var{idxVar};
           end
        end
    end
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get s (length state of contractile element)
typeStr = 's';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('s', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({''}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for s!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for s
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end


%% Get sdot (velocity state of contractile element)
typeStr = 'sdot';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('s', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'1/s'}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for sdot!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for sdot
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get LMTU (length of muscle-tendon-unit)
typeStr = 'LMTU';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'m'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for LMTU!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for LMTU
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get LCE (length of contractile element)
typeStr = 'LCE';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'m'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for LCE!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for LCE
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get LSEE (length of serial elastic element)
typeStr = 'LSEE';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'m'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for LSEE!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for LSEE
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get LdotMTU (velocity of muscle-tendon-unit)
typeStr = 'LdotMTU';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'m/s'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for LdotMTU!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for qdot
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get LdotCE (velocity of contractil element)
typeStr = 'LdotCE';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'m/s'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for LdotCE!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for LdotCE
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get LdotSEE (velocity of serial elastic element)
typeStr = 'LdotSEE';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'m/s'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for LdotSEE!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for LdotSEE
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get muscle forces
typeStr = 'muscleForce';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'N'}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for muscleForce!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        for iMus = 1 : height(simVarTableType)
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iMus}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, 'N'), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{iMus} =  variableTable.mean{idxVar};
               simVarTableType.var_extra{iMus}  =  variableTable.var{idxVar};
           end
        end
    end
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get CE forces
typeStr = 'CEForce';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'N'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for CEForce!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for CEForce
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get powers of muscle-tendon-unit
typeStr = 'musclePower';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'W'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for musclePower!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for musclePower
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get powers of contractile elements
typeStr = 'CEPower';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'W'}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for CEPower!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for CEPower
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get powers of seriel elastic elements
typeStr = 'SEEPower';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'W'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for SEEPower!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for SEEPower
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get metabolic rate of muscles
typeStr = 'muscleMetRate';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxa = model.extractState('a', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxa));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxa));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'W/kg'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType, 'umberger');

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for muscleMetRate!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for muscleMetRate
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end


%% Get joint power
typeStr = 'jointPower';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    unitStr = 'W/BW';
    bodyweight = model.bodymass*norm(model.gravity);

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        idxDof = model.extractState('q', settings.(typeStr));
        namesSym = model.states.name(model.idxSymmetry.xindex(idxDof));
        name = {settings.(typeStr){:}; namesSym{:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym1 = num2cell(model.idxSymmetry.xsign(idxDof));
        signSym2 = num2cell(ones(size(signSym1))); % use 1 as sign for this side which was already there
        signSym = {signSym2{:}; signSym1{:}}; % concat it alternating
        signSym = signSym(:);
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({unitStr}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for jointPower!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);

    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        variableTable = convertUnit(variableTable, typeStr, 'W', unitStr, 1/bodyweight); % convert from W to W/BW if needed
        for iMom = 1 : height(simVarTableType)
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iMom}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{iMom} =  variableTable.mean{idxVar};
               simVarTableType.var_extra{iMom}  =  variableTable.var{idxVar};
           end
        end
    end

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end


%% Get center of mass
typeStr = 'CoM';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))

    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones (simply the same names)
        name = {settings.(typeStr){:}; settings.(typeStr){:}}; % concat it alternating
        name = name(:);
        % Sign which has to be used to concatenate right and left (first and second half of gait cycle)
        signSym = num2cell(ones(size(name))); % use one for all
        % Unit displacement for translation
        supportedNames = {'CoM_x', 'CoM_y', 'CoM_z'};
        translationNames = {'pelvis_tx', 'pelvis_ty', 'pelvis_tz'};
        idxSimVar = zeros(size(settings.(typeStr)));
        for iName = 1 : numel(settings.(typeStr))
            idxSimVar(iName) = find(strcmp(supportedNames, settings.(typeStr){iName}));
        end
        idxTrans = model.extractState('q', translationNames(idxSimVar));
        unitdisplacement = zeros(1, numel(idxTrans));
        unitdisplacement(ismember(idxTrans, intersect(model.idxForward, idxTrans))) = 1;
        displacement = unitdisplacement * X(obj.idx.speed) * X(obj.idx.dur);
        displacement = [zeros(1, numel(idxTrans)); displacement];
        displacement = num2cell(displacement(:)); % concat it alternating
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'m'}, numel(name), 1);
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for CoM!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);

    % Get additional data
    % => We currently assume that there is no additional data for CoM
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);

    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym, displacement);
    end

    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get stance phase
typeStr = 'standing';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    if doSymmetry
        % Get names which are symmtric to the requested ones
        requestedNames = {'standing_r', 'standing_l'};
        % Do not need a sign change => Multiply by 1
        signSym = {1; 1};
        % Do always just return the right one
        settings.(typeStr) = {'standing_r'};
    else
        requestedNames = settings.(typeStr);
    end

    % Create table to get simulated GRF (table will not be saved in simVarTable)
    possibleNames = {'standing_r', 'standing_l'};
    requiredGRFNames = {'GRF_y_r'; 'GRF_y_l'};
    name = requiredGRFNames(ismember(possibleNames, requestedNames));
    type = repmat({'GRF'}, numel(name), 1);
    unit = repmat({'BW'}, numel(name), 1);
    simVarTableGRF = table(type, name, unit);
    
    % Get simulated GRF
    simVarTableGRF = obj.getSimData(X, simVarTableGRF);

    % Get full cycle by conctenating both sides
    if doSymmetry
        simVarTableGRF.mean = cell(height(simVarTableGRF), 1);
        simVarTableGRF.var  = cell(height(simVarTableGRF), 1);
        simVarTableGRF.mean_extra = cell(height(simVarTableGRF), 1);
        simVarTableGRF.var_extra  = cell(height(simVarTableGRF), 1);
        simVarTableGRF = getFullCycleFromHalf(simVarTableGRF, signSym);
    end
        
    % Get standing and HS and TO events
    nNames = length(settings.(typeStr));
    nameSt = settings.(typeStr)(:);
    sim = cell(nNames*3, 1);
    nameIdx = cell(nNames*2, 1);
    for iName = 1 : nNames
        verGRF = simVarTableGRF.sim{iName} * model.bodymass * norm(model.gravity); % vertical GRF in N
        [standing, idxHSs, idxTOs] = getStancePhase(verGRF);
        sim{iName} = standing;
        sim{iName*nNames+1} = idxHSs;
        sim{iName*nNames+2} = idxTOs;
        nameIdx{1+(iName-1)*2} = ['HS_' nameSt{iName}(end)];
        nameIdx{2+(iName-1)*2} = ['TO_' nameSt{iName}(end)];
    end
    
    % Create simVarTable to for standing, idxHSs, and idxTOs
    % Standing
    typeSt = repmat({'standing'}, nNames, 1);
    unitSt = repmat({'boolean'}, nNames, 1);
    % idxHSs and idxTOs
    typeIdx = repmat({'standingEvent'}, nNames*2, 1);
    unitIdx = repmat({'index'}, nNames*2, 1);
    % Combine the vectors
    name = [nameSt; nameIdx];
    type = [typeSt; typeIdx];
    unit = [unitSt; unitIdx];
    simVarTableType = table(type, name, unit, sim);
    
    % Get tracking data (mean and var) if it was tracked
    % => There is currently no tracking method for standing!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for standing
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    
    % Add simVarTableType to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];
    
end

%% Get foot angles
typeStr = 'footAngle';
if isfield(settings, typeStr) && ~isempty(settings.(typeStr)) && iscellstr(settings.(typeStr))
    
    % Create table to get simulated data
    if doSymmetry
        % Get names which are symmtric to the requested ones
        name = {'angle_r'; 'angle_l'};
        % Do not need a sign change => Multiply by 1
        signSym = {1; 1};
    else
        name = settings.(typeStr)(:);
    end
    type = repmat({typeStr}, numel(name), 1);
    unit = repmat({'deg'}, numel(name), 1);
    simVarTableType = table(type, name, unit);
    
    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);
    
    % Get tracking data (mean and var) if it was tracked 
    % => There is currently no tracking method for footAngle!
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var  = cell(height(simVarTableType), 1);
    
    % Get additional data
    % => We currently assume that there is no additional data for footAngle
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra  = cell(height(simVarTableType), 1);
    
    % Get full cycle by conctenating both sides
    if doSymmetry
       simVarTableType = getFullCycleFromHalf(simVarTableType, signSym);
    end

    % If standing is extracted get angle at HS and TO
    typeStrSt = 'standing';
    if isfield(settings, typeStrSt) && ~isempty(settings.(typeStrSt)) && iscellstr(settings.(typeStrSt))
        
        % Create table including simulated data
        doRight = ismember('angle_r', settings.(typeStr)) && ismember('standing_r', settings.(typeStrSt));
        doLeft = ismember('angle_l', settings.(typeStr)) && ismember('standing_l', settings.(typeStrSt));
        name = cell((doRight+doLeft)*2, 1);
        sim = cell((doRight+doLeft)*2, 1);
        % right
        if doRight
            name{1} = insertAfter('angle_r', '_', 'HS_'); 
            name{2} = insertAfter('angle_r', '_', 'TO_');
            iHS = simVarTable.sim{strcmp(simVarTable.name, 'HS_r')};
            iTO = simVarTable.sim{strcmp(simVarTable.name, 'TO_r')};
            sim{1} = simVarTableType.sim{strcmp(simVarTableType.name, 'angle_r')}(iHS);
            sim{2} = simVarTableType.sim{strcmp(simVarTableType.name, 'angle_r')}(iTO);
        end
        % left
        if doLeft
            name{doRight*2+1} = insertAfter('angle_l', '_', 'HS_');
            name{doRight*2+2} = insertAfter('angle_l', '_', 'TO_');
            iHS = simVarTable.sim{strcmp(simVarTable.name, 'HS_l')};
            iTO = simVarTable.sim{strcmp(simVarTable.name, 'TO_l')};
            sim{doRight*2+1} = simVarTableType.sim{strcmp(simVarTableType.name, 'angle_l')}(iHS);
            sim{doRight*2+2} = simVarTableType.sim{strcmp(simVarTableType.name, 'angle_l')}(iTO);
        end
        type = repmat({'footAngleEvent'}, numel(name), 1);
        unit = repmat({'deg'}, numel(name), 1);
        simVarTableHSTO = table(type, name, unit, sim);
        
        % Get tracking data (mean and var) if it was tracked
        % => There is currently no tracking method for footAngle!
        simVarTableHSTO.mean = cell(height(simVarTableHSTO), 1);
        simVarTableHSTO.var  = cell(height(simVarTableHSTO), 1);
        
        % Get additional data
        % => We currently assume that there is no additional data for footAngle
        simVarTableHSTO.mean_extra = cell(height(simVarTableHSTO), 1);
        simVarTableHSTO.var_extra  = cell(height(simVarTableHSTO), 1);
        
        % Add table to simVarTable
        simVarTableType = [simVarTableType; simVarTableHSTO];
        
    end
    
    % Add simVarTable to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end

%% Get duration
typeStr = 'duration';
if isfield(settings, typeStr) && settings.(typeStr)

    % Create table to get simulated data
    name = {typeStr};
    type = {typeStr};
    unit = {'s'};
    simVarTableType = table(type, name, unit);

    % Get simulated data
    simVarTableType = obj.getSimData(X, simVarTableType);

    % Get tracking data (mean and var) if it was tracked
    simVarTableType.mean = cell(height(simVarTableType), 1);
    simVarTableType.var = cell(height(simVarTableType), 1);
    idxObjTerm = find(strcmp(objTermsNames, 'trackDuration'));
    if length(idxObjTerm) == 1
        variables = obj.objectiveTerms(idxObjTerm).varargin{:}.variables;
        for iDur = 1 : height(simVarTableType) % There should be only one duration entry, but this does not hurt...
           idxVar = find(strcmp(variables.type, typeStr) & strcmp(variables.name, simVarTableType.name{iDur}));
           if ~isempty(idxVar)
               assert(strcmp(variables.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variables.unit{idxVar}, typeStr);
               simVarTableType.mean{iDur} = variables.mean{idxVar};
               simVarTableType.var{iDur}  = variables.var{idxVar};
           end
        end
    end

    % Get additional data
    simVarTableType.mean_extra = cell(height(simVarTableType), 1);
    simVarTableType.var_extra = cell(height(simVarTableType), 1);
    if ~isempty(variableTable)
        for iDur = 1 : height(simVarTableType) % There should be only one duration entry, but this does not hurt...
           idxVar = find(strcmp(variableTable.type, typeStr) & strcmp(variableTable.name, simVarTableType.name{iDur}));
           if ~isempty(idxVar)
               assert(strcmp(variableTable.unit{idxVar}, unitStr), 'Unit %s is not supported for type %s.', variableTable.unit{idxVar}, typeStr);
               simVarTableType.mean_extra{iDur} = variableTable.mean{idxVar};
               simVarTableType.var_extra{iDur}  = variableTable.var{idxVar};
           end
        end
    end

    % Get full cycle by multiplying the duration by two
    if doSymmetry
       simVarTableType.sim{1} = simVarTableType.sim{1} *2;
    end

    % Add simVarTableType to table containing all data
    simVarTable = [simVarTable; addDirPosSegCol(simVarTable, simVarTableType)];

end


end


%> @cond DO_NOT_DOCUMENT
%======================================================================
%> @brief Function to concatenate both sides to get full cycle
%>
%> @details
%> It assumes that two consecutive rows belong to each other.
%>
%> @param   simVarTableHalf  Table: Table containing rows for both side with the half of the cycle
%> @param   signSym          Cell: Signs for every row to concatenate the rows
%> @param   displacement     (optional) Cell: Displacement of simulated data which has to be added when concatenating
%> @retval  simVarTableFull  Table: Table containing rows for one side with the full cycle
%======================================================================
function simVarTableFull = getFullCycleFromHalf(simVarTableHalf, signSym, displacement)

% Get basic columns
type = simVarTableHalf.type(1:2:end);
name = simVarTableHalf.name(1:2:end);
unit = simVarTableHalf.unit(1:2:end);

% Concatenate data
sim        = cell(numel(type), 1);
mean       = cell(numel(type), 1);
var        = cell(numel(type), 1);
mean_extra = cell(numel(type), 1);
var_extra   = cell(numel(type), 1);
for iVar = 1 : numel(type)

    % Simulated data
    if nargin < 3
        % Concatenate with without displacement
        sim{iVar} = [simVarTableHalf.sim{iVar*2-1}; simVarTableHalf.sim{iVar*2} * signSym{iVar*2}];
    else
        % Concatenate with with displacement
        sim{iVar} = [simVarTableHalf.sim{iVar*2-1}; simVarTableHalf.sim{iVar*2} * signSym{iVar*2} + displacement{iVar*2}];
    end

    % Add NaNs if data of only one side was tracked
    if length(simVarTableHalf.mean{iVar*2-1}) ~= length(simVarTableHalf.mean{iVar*2})
       if isempty(simVarTableHalf.mean{iVar*2-1})
          simVarTableHalf.mean{iVar*2-1} = nan(size(simVarTableHalf.mean{iVar*2}));
          simVarTableHalf.var{iVar*2-1}  = nan(size(simVarTableHalf.var{iVar*2}));
       end
       if isempty(simVarTableHalf.mean{iVar*2})
          simVarTableHalf.mean{iVar*2} = nan(size(simVarTableHalf.mean{iVar*2-1}));
          simVarTableHalf.var{iVar*2}  = nan(size(simVarTableHalf.var{iVar*2-1}));
       end
    end

    % Add NaNs if data of only one side was given in the additional data
    if length(simVarTableHalf.mean_extra{iVar*2-1}) ~= length(simVarTableHalf.mean_extra{iVar*2})
       if isempty(simVarTableHalf.mean_extra{iVar*2-1})
          simVarTableHalf.mean_extra{iVar*2-1} = nan(size(simVarTableHalf.mean_extra{iVar*2}));
          simVarTableHalf.var_extra{iVar*2-1}  = nan(size(simVarTableHalf.var_extra{iVar*2}));
       end
       if isempty(simVarTableHalf.mean_extra{iVar*2})
          simVarTableHalf.mean_extra{iVar*2} = nan(size(simVarTableHalf.mean_extra{iVar*2-1}));
          simVarTableHalf.var_extra{iVar*2}  = nan(size(simVarTableHalf.var_extra{iVar*2-1}));
       end
    end

    % Tracked and additional data
    if nargin < 3
        % Concatenate with without displacement
        mean{iVar}       = [simVarTableHalf.mean{iVar*2-1};       simVarTableHalf.mean{iVar*2}       * signSym{iVar*2}];
        var{iVar}        = [simVarTableHalf.var{iVar*2-1};        simVarTableHalf.var{iVar*2}        * signSym{iVar*2}];
        mean_extra{iVar} = [simVarTableHalf.mean_extra{iVar*2-1}; simVarTableHalf.mean_extra{iVar*2} * signSym{iVar*2}];
        var_extra{iVar}  = [simVarTableHalf.var_extra{iVar*2-1};  simVarTableHalf.var_extra{iVar*2}  * signSym{iVar*2}];
    else
        % Do not concatenate but add NaNs do not know the displacement
        mean{iVar}       = [simVarTableHalf.mean{iVar*2-1};       nan(size(simVarTableHalf.mean{iVar*2-1}))];
        var{iVar}        = [simVarTableHalf.var{iVar*2-1};        nan(size(simVarTableHalf.var{iVar*2-1}))];
        mean_extra{iVar} = [simVarTableHalf.mean_extra{iVar*2-1}; nan(size(simVarTableHalf.mean_extra{iVar*2-1}))];
        var_extra{iVar}  = [simVarTableHalf.var_extra{iVar*2-1};  nan(size(simVarTableHalf.var_extra{iVar*2-1}))];
    end

end

% Make table
simVarTableFull = table(type, name, unit, sim, mean, var, mean_extra, var_extra);

end

%======================================================================
%> @brief Function to add nans to symmetric data for seconds half
%>
%> @details
%> For acc, gyro, and marker, it adds NaNs for the second half of the gait cycle
%> since we can not assume symmetric sensor/marker positions.
%>
%> @param   simVarTableHalf  Table: Table containing rows for one side with the half of the cycle
%> @retval  simVarTableFull  Table: Table containing rows for one side with the full cycle (NaNs for second half)
%======================================================================
function simVarTableFull = addNaNsForSecondHalf(simVarTableHalf)

% Copy the table
simVarTableFull = simVarTableHalf;

% Add NaNs
for iVar = 1 : height(simVarTableHalf)
   simVarTableFull.sim{iVar}        = [simVarTableHalf.sim{iVar};  nan(size(simVarTableHalf.sim{iVar}))];
   simVarTableFull.mean{iVar}       = [simVarTableHalf.mean{iVar}; nan(size(simVarTableHalf.mean{iVar}))];
   simVarTableFull.var{iVar}        = [simVarTableHalf.var{iVar};  nan(size(simVarTableHalf.var{iVar}))];
   simVarTableFull.mean_extra{iVar} = [simVarTableHalf.mean_extra{iVar};  nan(size(simVarTableHalf.mean_extra{iVar}))];
   simVarTableFull.var_extra{iVar}  = [simVarTableHalf.var_extra{iVar};   nan(size(simVarTableHalf.var_extra{iVar}))];
end

end


%======================================================================
%> @brief Function to add direction, position, and segment columns to simVarTable
%>
%> @details
%> It adds columns with NaNs for direction and position and empty chars for segment
%> if they these columns in simVarTableAll.
%>
%> @param   simVarTableAll  Table: Table containing all data
%> @param   simVarTable     Table: New table which should be added to simVarTableAll
%> @param   simVarTable     Table: Adapted table which should be added to simVarTableAll
%======================================================================
function simVarTable = addDirPosSegCol(simVarTableAll, simVarTable)

if isempty(simVarTable)
    % Do not add the columns since concatenation would not work afterwards
    return;
end

colHeaders = {'direction', 'position'};

for iCol = 1 : numel(colHeaders)
    
    doColumnExist = find(strcmp(colHeaders{iCol},simVarTableAll.Properties.VariableNames));
    if (doColumnExist)
        nDim = size(simVarTableAll.(colHeaders{iCol}), 2);
        simVarTable.(colHeaders{iCol}) = nan(size(simVarTable,1), nDim);
    end
    
end

colHeader = 'segment';
doColumnExist = find(strcmp(colHeader,simVarTableAll.Properties.VariableNames));

if (doColumnExist)
    simVarTable.(colHeader) = repmat({''}, height(simVarTable), 1);
end

end
%> @endcond
