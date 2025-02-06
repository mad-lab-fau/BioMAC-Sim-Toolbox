% ======================================================================
%> @file @TrackingData/writeMotionToOsim.m
%> @brief TrackingData function to write motion from tracking data object to OpenSim files
%> @details
%> Details: TrackingData::writeMotionToOsim()
%>
%> @author Marlies Nitschke
%> @date July, 2019
% ======================================================================

% ======================================================================
%> @brief Function to write motion from tracking data object to OpenSim files
%>
%> @details
%> Saves the motion to the following files:
%> 1. <filename>_kinematics.mot containing the joint angles
%> 2. <filename>_kinetics_GRFs.mot containing GRFs
%> 3. <filename>_kinetics_moments.sto containing moments
%>
%> No resampling is applied.
%>
%> @param  obj      TrackingData class object which should saved to OpenSim files.
%>                  It has to contain the following fields:
%>                   | Type        | Name                                                             | 
%>                   |-------------|------------------------------------------------------------------|
%>                   | duration    |                                                                  |
%>                   | angle       | all joint angles of the model                                    |
%>                   | translation | all translations of the model                                    |
%>                   | GRF         | 'GRF_x_r', 'GRF_y_r', 'GRF_z_r', 'GRF_x_l', 'GRF_y_l', 'GRF_z_l' |
%>                   | CoP         | 'CoP_x_r', 'CoP_y_r', 'CoP_z_r', 'CoP_x_l', 'CoP_y_l', 'CoP_z_l' |
%>                   | GRM         | 'GRM_x_r', 'GRM_y_r', 'GRM_z_r', 'GRM_x_l', 'GRM_y_l', 'GRM_z_l' |
%>                   | moment      | all DOFs of the model                                            |
%>
%> @param  model    Model: Model associated with the tracking data. This is
%>                  required to obtain names of DOFs, bodymass, and gravity
%> @param  filename String: Filename to save the results. See details for the used postfixes
% ======================================================================
function writeMotionToOsim(obj, model, filename)

% Get time vector
nSamples = obj.nSamples;
duration = obj.variables.mean{strcmp(obj.variables.type,'duration')};
times = (0:nSamples-1)/nSamples*duration;

% Get kinematics
idxq = find(strcmp(obj.variables.type,'angle') | strcmp(obj.variables.type,'translation'));
assert(numel(idxq)==model.nDofs, 'Kinematics of all DOFs have to be in the TrackingData object.');
q = [obj.variables.mean{idxq}];
names = {obj.variables.name{idxq}};
pelvis_t_xyz = {'pelvis_tx', 'pelvis_ty', 'pelvis_tz'};
idxConvert = ~ismember(names, pelvis_t_xyz); % all which are not translation
q(:, idxConvert) = q(:, idxConvert) / pi * 180;

% Save kinematics
filenameKinem = [filename '_kinematics.mot'];
inDegrees = 1;
writeMotSto(times, q, names, filenameKinem, inDegrees);

% Get GRF data including CoP and GRM
type = 'GRF';
namesGiven = {'GRF_x_r', 'GRF_y_r', 'GRF_z_r', 'GRF_x_l', 'GRF_y_l', 'GRF_z_l'};
grf = getGRFData(obj, type, namesGiven);
type = 'CoP';
namesGiven = {'CoP_x_r', 'CoP_y_r', 'CoP_z_r', 'CoP_x_l', 'CoP_y_l', 'CoP_z_l'};
CoP = getGRFData(obj, type, namesGiven);
type = 'GRM';
namesGiven = {'GRM_x_r', 'GRM_y_r', 'GRM_z_r', 'GRM_x_l', 'GRM_y_l', 'GRM_z_l'};
grm = getGRFData(obj, type, namesGiven);

names = {'ground_force_vx','ground_force_vy','ground_force_vz', ...  % right foot's GRF vector
    'ground_force_px','ground_force_py','ground_force_pz', ...       % right foot's center of pressure
    'l_ground_force_vx','l_ground_force_vy','l_ground_force_vz', ... % left  foot's GRF vector
    'l_ground_force_px','l_ground_force_py','l_ground_force_pz', ... % left  foot's center of pressure
    'ground_torque_x','ground_torque_y','ground_torque_z', ...       % right foot's moment vector
    'l_ground_torque_x','l_ground_torque_y','l_ground_torque_z'};    % left  foot's moment vector

convert_BW_to_N = norm(model.gravity) * model.bodymass;  % convert BW to N
dataGRF(:, 1:3)   = grf(:, 1:3) * convert_BW_to_N; % right foot's GRF vector
dataGRF(:, 4:6)   = CoP(:, 1:3);                   % right foot's center of pressure
dataGRF(:, 7:9)   = grf(:, 4:6) * convert_BW_to_N; % left foot's GRF vector
dataGRF(:, 10:12) = CoP(:, 4:6);                   % left foot's center of pressure
dataGRF(:, 13:15) = grm(:, 1:3) * convert_BW_to_N; % right foot's GRF vector
dataGRF(:, 16:18) = grm(:, 4:6) * convert_BW_to_N; % left foot's GRF vector

% Save GRF data
filenameGRF = [filename '_kinetics_GRFs.mot'];
inDegrees = 0;
writeMotSto(times, dataGRF, names, filenameGRF, inDegrees);

% Get moments
idxMom = find(strcmp(obj.variables.type,'moment'));
assert(numel(idxMom)==model.nDofs, 'Moments of all DOFs have to be in the TrackingData object.');
M = [obj.variables.mean{idxMom}];
namesDOF = {obj.variables.name{idxMom}};
names = cell(model.nDofs, 1);
for iDof = 1 : model.nDofs
    if ismember(model.dofs.Properties.RowNames{iDof}, pelvis_t_xyz)
        names{iDof} = sprintf('\t%s_force', namesDOF{iDof});
    else
        names{iDof} = sprintf('\t%s_moment', namesDOF{iDof});
    end
end

% Write to file
filenameMom = [filename '_kinetics_moments.sto'];
inDegrees = 0;
writeMotSto(times, M, names, filenameMom, inDegrees);


end


%> @cond DO_NOT_DOCUMENT
%======================================================================
%> @brief Function to extract and sort GRF, CoP and GRM data from the variables table
%======================================================================
function data = getGRFData(trackingData, type, namesGiven)

% Get indices in var table
idx = find(strcmp(trackingData.variables.type,type));
assert(numel(idx)==6, ['All 6 ' type ' should be given in the TrackingData object.']);

% Get data and names
data  = [trackingData.variables.mean{idx}];
names = {trackingData.variables.name{idx}};

% Sort data
[isMem, idxSort] = ismember(names,namesGiven);
assert(sum(isMem)==6, ['All 6 ' type ' should be given in the TrackingData object with the correct names.']);
data(:, 1:end) = data(:, idxSort); % Sort correctly

end
%> @endcond

