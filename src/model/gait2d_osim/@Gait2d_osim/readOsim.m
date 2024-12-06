%======================================================================
%> @file @Gait2d_osim/readOsim.m
%> @brief Gait2d_osim function to load OpenSim model in Matlab as structure
%> @details
%> Details: Gait2d_osim::readOsim()
%> The code is copied from gait3d entirely, with adjustments for gait2d_osim
%> 
%> @author Eva Dorschky, Ton van den Bogert, Marlies Nitschke
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Function to load OpenSim model in Matlab as structure
%> @protected 
%> @static
%>
%> @details
%> Saves the model struct using strrep(osimfile, '.osim', '.mat').
%>
%> @param  osimfile     String: Filename of opensim file including path and extension
%> @retval model        Struct: Parameters of OpenSim model
%======================================================================
function model = readOsim(osimfile)
%GAI3D_READOSIM Initializes a MATLAB struct from an opensim file.
%   osimfile: The OpenSIM model file from which to initialize
try
    fileID = fopen(osimfile);
    fileID_read =fread(fileID);
catch
    error('Opensim file not found. Please provide the full path to the Opensim file or change to the according folder.')
end

md = java.security.MessageDigest.getInstance('SHA-256');
osim_sha256 = typecast(md.digest(uint8(fileID_read))', 'uint8');

warning('Opensim library let Matlab crash if IPOPT exits with an exception. Open a new matlab session if you want to run optimizations with IPOPT.');

% import OpenSim namespace
import org.opensim.modeling.*;

% load the model - For HPC: dont use 'which' when the given filepath is absolute
if osimfile(1) == filesep || osimfile(1) == '~' || osimfile(1) == '$'
    osim_path = osimfile;
else
    osim_path = which(osimfile);
end

Mod = Model(osim_path);
% initialize the system and get the initial state
state = Mod.initSystem();

% equilibrate all muscles
Mod.equilibrateMuscles(state);

% get version of OpenSim
opensimVersion = char(org.opensim.modeling.opensimCommon.GetVersion());
opensimVersion = opensimVersion(1:3); % Only first chars ar important
opensimVersion33 = '3.3';
opensimVersion40 = '4.0';
opensimVersion41 = '4.1';
opensimVersion43 = '4.3';
opensimVersion45 = '4.5';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% basic model parameters
model.osim.name = char(Mod.getName());
save_path = osim_path(strfind(char(Mod.getInputFileName()), 'model')-1:end);
model.osim.file = save_path;
file = dir(char(Mod.getInputFileName()));
model.osim.modified = file.date;

% gravity
tmpVec3f = Mod.getGravity();
model.gravity = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of joints
JointSet = Mod.getJointSet();
nJoints = JointSet.getSize();

% loop through joints
for ijoints = 1:nJoints
    
    currentJoint = JointSet.get(ijoints-1);
    
    % basics
    joints(ijoints,1).name = char(currentJoint.getName());
    
    % location in parent segment
    tmpVec3f = Vec3();
    switch opensimVersion
        case opensimVersion33
            currentJoint.getLocationInParent(tmpVec3f);
            joints(ijoints,1).location = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
        case opensimVersion40
            % I couldn't find an according function in opensim 4
            A = fileread(osimfile);
            b = strfind(A,'<location_in_parent>');
            c = strfind(A,'</location_in_parent>');
            l = A(b(ijoints):c(ijoints));
            l = l(21:end-1);
            l = str2num(l);
            joints(ijoints,1).location = l; 
        case {opensimVersion41, opensimVersion43, opensimVersion45}
            % OpenSim 4.1. Asuming that the first PhysicalOffsetFrame is
            % the parent
            parentFrame = currentJoint.getParentFrame;
            potentialParentOffsetFrame = currentJoint.get_frames(0);
            assert(strcmp(parentFrame.getName(), potentialParentOffsetFrame.getName()), ...
                'First PhysicalOffsetFrame is not parent.');
            tmpVec3f = potentialParentOffsetFrame.get_translation();
            joints(ijoints,1).location = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    end
    
    % orientation in parent segment (seldom used)
    tmpVec3f = Vec3();
    switch opensimVersion
        case opensimVersion33
            currentJoint.getOrientationInParent(tmpVec3f);
            joints(ijoints,1).orientation = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
        case opensimVersion40
            % I couldn't find an according function in opensim 4
            A = fileread(osimfile);
            b = strfind(A,'<orientation_in_parent>');
            c = strfind(A,'</orientation_in_parent>');
            l = A(b(ijoints):c(ijoints));
            l = l(24:end-1);
            l = str2num(l);
            
            joints(ijoints,1).orientation= l;
        case {opensimVersion41, opensimVersion43, opensimVersion45}
            tmpVec3f = potentialParentOffsetFrame.get_orientation();
            joints(ijoints,1).orientation = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    end
    
    % spatial transforms
    switch char(currentJoint.getConcreteClassName)
        case 'CustomJoint'
            % get transforms
            currentJoint = CustomJoint.safeDownCast(currentJoint);
            transforms = currentJoint.getSpatialTransform();
            
            % rotations
            r1 = transforms.getTransformAxis(0);
            tmpVec3f = r1.getAxis();
            joints(ijoints,1).r1_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
            joints(ijoints,1).r1_coefs = osimfunction2coefs(r1.getFunction());
            
            r2 = transforms.getTransformAxis(1);
            tmpVec3f = r2.getAxis();
            joints(ijoints,1).r2_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
            joints(ijoints,1).r2_coefs = osimfunction2coefs(r2.getFunction());
            
            r3 = transforms.getTransformAxis(2);
            tmpVec3f = r3.getAxis();
            joints(ijoints,1).r3_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
            joints(ijoints,1).r3_coefs = osimfunction2coefs(r3.getFunction());
            
            % translations
            t1 = transforms.getTransformAxis(3);
            tmpVec3f = t1.getAxis();
            joints(ijoints,1).t1_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
            joints(ijoints,1).t1_coefs = osimfunction2coefs(t1.getFunction());
            
            t2 = transforms.getTransformAxis(4);
            tmpVec3f = t2.getAxis();
            joints(ijoints,1).t2_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
            joints(ijoints,1).t2_coefs = osimfunction2coefs(t2.getFunction());
            
            t3 = transforms.getTransformAxis(5);
            tmpVec3f = t3.getAxis();
            joints(ijoints,1).t3_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
            joints(ijoints,1).t3_coefs = osimfunction2coefs(t3.getFunction());
            
        case 'PinJoint' % rotation around Z
            % rotations
            joints(ijoints,1).r1_axis = [0, 0, 1];
            joints(ijoints,1).r1_coefs = [0, 0, 0, 0, 0];
            
            joints(ijoints,1).r2_axis = [0, 1, 0];
            joints(ijoints,1).r2_coefs = [0, 0, 0, 0, 0];
            
            joints(ijoints,1).r3_axis = [1, 0, 0];
            joints(ijoints,1).r3_coefs = [0, 0, 0, 0, 0];
            
            % translations
            joints(ijoints,1).t1_axis = [1, 0, 0];
            joints(ijoints,1).t1_coefs = [0, 0, 0, 0, 0];
            
            joints(ijoints,1).t2_axis = [0, 1, 0];
            joints(ijoints,1).t2_coefs = [0, 0, 0, 0, 0];
            
            joints(ijoints,1).t3_axis = [0, 0, 1];
            joints(ijoints,1).t3_coefs = [0, 0, 0, 0, 0];
            
            
        case 'WeldJoint' % no movement
            % rotations
            joints(ijoints,1).r1_axis = nan(1, 3);
            joints(ijoints,1).r1_coefs = nan(1, 5);
            
            joints(ijoints,1).r2_axis = nan(1, 3);
            joints(ijoints,1).r2_coefs = nan(1, 5);
            
            joints(ijoints,1).r3_axis = nan(1, 3);
            joints(ijoints,1).r3_coefs = nan(1, 5);
            
            % translations
            joints(ijoints,1).t1_axis = nan(1, 3);
            joints(ijoints,1).t1_coefs = nan(1, 5);
            
            joints(ijoints,1).t2_axis = nan(1, 3);
            joints(ijoints,1).t2_coefs = nan(1, 5);
            
            joints(ijoints,1).t3_axis = nan(1, 3);
            joints(ijoints,1).t3_coefs = nan(1, 5);
            
        otherwise
            error('Joint type ''%s'' not known', currentJoint.getConcreteClassName);
    end
   
    % name of parent segment this joint is attached to
    switch opensimVersion
        case opensimVersion33
            joints(ijoints,1).parent_segment=char(currentJoint.getParentBody().getName());
        case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
            parentSegment = char(currentJoint.getParentFrame);
            parentSegment = erase(parentSegment,'_offset');
            joints(ijoints,1).parent_segment=parentSegment;
    end
    
    % name of segment this joint belongs to
    switch opensimVersion
        case opensimVersion33
            joints(ijoints,1).segment=char(currentJoint.getBody().getName());
        case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
            joints(ijoints,1).segment= erase(char(currentJoint.getChildFrame),'_offset');
    end
end

% Make table
joints = struct2table(joints);
% Set name as row name
joints.Properties.RowNames = joints.name;
joints.name = [];
% Set property
model.joints = joints;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEGMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the set of segments
BodySet = Mod.getBodySet();
nSegments = BodySet.getSize();
startSegment = 1;
if ~strcmp(char(BodySet.get(0).getName()),'ground')
    nSegments = nSegments+1;
    startSegment = 2;
    segments(1,1).name ='ground';
    segments(1,1).mass = 0;
    segments(1,1).mass_center = [0 0 0];
    segments(1,1).inertia = zeros(3);
    segments(1,1).parent_joint= '';
end

% loop through segments
for isegment = startSegment:nSegments
    currentSegment = BodySet.get(isegment-startSegment+1-1);
   
    % name and mass
    segments(isegment,1).name = char(currentSegment.getName());
    segments(isegment,1).mass = currentSegment.getMass();
    
    % mass center
    tmpVec3 = ArrayDouble.createVec3([1 0 0]);
    switch opensimVersion
        case opensimVersion33
            currentSegment.getMassCenter(tmpVec3);
        case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
            tmpVec3 = currentSegment.get_mass_center;
    end
    tmpAD = ArrayDouble.getValuesFromVec3(tmpVec3);
    segments(isegment,1).mass_center = [tmpAD.getitem(0) tmpAD.getitem(1) tmpAD.getitem(2)];
    
    % inertia matrix
    tmpMat33=Mat33(0,0,0,0,0,0,0,0,0);
    switch opensimVersion
        case opensimVersion33
            currentSegment.getInertia(tmpMat33);
            segments(isegment,1).inertia = zeros(3,3);
            for x=1:3
                for y=1:3
                    segments(isegment,1).inertia(x,y)=tmpMat33.get(x-1,y-1);
                end
            end
        case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
            %"The elements of the inertia tensor (Vec6) as [Ixx Iyy Izz Ixy Ixz Iyz] " "measured about the mass_center and not the body origin."
            tmpVec6 = currentSegment.get_inertia;
            segments(isegment,1).inertia = zeros(3,3);
            segments(isegment,1).inertia(1,1)=tmpVec6.get(0);
            segments(isegment,1).inertia(2,2)=tmpVec6.get(1);
            segments(isegment,1).inertia(3,3)=tmpVec6.get(2);
            segments(isegment,1).inertia(1,2)=tmpVec6.get(3);
            segments(isegment,1).inertia(1,3)=tmpVec6.get(4);
            segments(isegment,1).inertia(2,3)=tmpVec6.get(5);
            segments(isegment,1).inertia(2,1)=tmpVec6.get(3);
            segments(isegment,1).inertia(3,1)=tmpVec6.get(4);
            segments(isegment,1).inertia(3,2)=tmpVec6.get(5);
    end
    
    
    % name of parent joint this segment is attached to
    segments(isegment,1).parent_joint='';
    switch opensimVersion
        case opensimVersion33
            if currentSegment.hasJoint()
                segments(isegment,1).parent_joint=char(currentSegment.getJoint().getName());
            end
        case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
            % no according function for OpenSim 4 was found
            for ijoints = 1:nJoints
                currentJointChildName = erase(char(JointSet.get(ijoints-1).getChildFrame().getName()),'_offset');
                if strcmp(currentJointChildName, currentSegment.getName())
                    segments(isegment,1).parent_joint = char(JointSet.get(ijoints-1).getName());
                    break;
                end
            end
    end
    
end

% Make table
segments = struct2table(segments);
% Set name as row name
segments.Properties.RowNames = segments.name;
segments.name = [];
% Set property
model.segments = segments;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of DOFs
CoordSet = Mod.getCoordinateSet();
nDofs = CoordSet.getSize();

counter=0;
dofs_segment_list = zeros(nSegments,nDofs);
% loop through dofs
for idof = 1:nDofs
    
    currentDof = CoordSet.get(idof-1);
   
    % remove all locked DoFs from the model
    % CAREFUL, this changes indices
    
    if ~currentDof.getLocked(state)
        
        counter=counter+1;
        switch opensimVersion
            case opensimVersion33
                currentBody = char(currentDof.getJoint.getBody().getName());
            case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
                currentBody = char(currentDof.getJoint.getChildFrame);
        end
        
        if contains(currentBody, '_offset')
            currentBody = erase(currentBody,'_offset');
        end
        iBody = find(strcmp(segments.Properties.RowNames,currentBody));
        dofs_segment_list(iBody,counter)= 1; %matrix showing the dofs for every segment
        
        % basics
        dofs(counter,1).name = char(currentDof.getName());
        
        % range of motion in radians
        dofs(counter,1).range_osim(1) = currentDof.getRangeMin();
        dofs(counter,1).range_osim(2) = currentDof.getRangeMax();
        
        % get default value (= neutral position)
        dofs(counter, 1).neutral_position = currentDof.get_default_value();
        
        % joint
        dofs(counter,1).joint    = char(currentDof.getJoint);
        
        % name of joint this DOF belongs to
        % currently commented out, because currentSegment was not correct, it was the very last segment always
        % we don't need this field now, maybe later if we automatically generate Autolev code from the Opensim model
        % dofs(counter,1).joint=char(currentSegment.getJoint().getName());
    end
end
nDofs = counter;

% Make table
dofs = struct2table(dofs);
% Set name as row name
dofs.Properties.RowNames = dofs.name;
dofs.name = [];

% Set property
model.dofs = dofs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUSCLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of muscles
MuscleSet = Mod.getMuscles();
nMus = MuscleSet.getSize();
% loop through muscles
for imus = 1:nMus
    currentMuscle = MuscleSet.get(imus-1);
    
    % basic muscle properties
    muscles(imus,1).name = char(currentMuscle.getName());
    muscles(imus,1).fmax = currentMuscle.getMaxIsometricForce();
    muscles(imus,1).lceopt = currentMuscle.getOptimalFiberLength();
    muscles(imus,1).lslack = currentMuscle.getTendonSlackLength();
    muscles(imus,1).pennatopt = currentMuscle.getPennationAngleAtOptimalFiberLength();
    muscles(imus,1).vmax = currentMuscle.getMaxContractionVelocity();
    
    musType = char(currentMuscle.getConcreteClassName);
    switch musType
        case 'Thelen2003Muscle'
            currentMuscle = Thelen2003Muscle.safeDownCast(currentMuscle);
            switch opensimVersion
                case opensimVersion33
                    muscles(imus,1).kactive = currentMuscle.getKshapeActive();
                    muscles(imus,1).flen = currentMuscle.getFlen();
                case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
                    muscles(imus,1).kactive = currentMuscle.get_KshapeActive();
                    muscles(imus,1).flen = currentMuscle.get_Flen();
            end
        case 'Millard2012EquilibriumMuscle'
            currentMuscle = Millard2012EquilibriumMuscle.safeDownCast(currentMuscle);
            muscles(imus,1).kactive = sqrt(currentMuscle.getActiveForceLengthCurve().getShallowAscendingSlope);
            muscles(imus,1).flen =  currentMuscle.getForceVelocityCurve.getMaxEccentricVelocityForceMultiplier;
    end
    % Thelen 2003 muscle properties
    muscles(imus,1).tact = currentMuscle.getActivationTimeConstant();
    % muscles(imus,1).minact = currentMuscle.getMinimumActivation();
    muscles(imus,1).tdeact = currentMuscle.getDeactivationTimeConstant();
    % muscles(imus,1).fmaxts = currentMuscle.getFmaxTendonStrain();
    % muscles(imus,1).fmaxms = currentMuscle.getFmaxMuscleStrain();
    
    % muscles(imus,1).kpassive = currentMuscle.getKshapePassive();
    % muscles(imus,1).af = currentMuscle.getAf();
        
            
         
    % muscles(imus,1).velthres = currentMuscle.getForceVelocityExtrapolationThreshold();
    % muscles(imus,1).minflen = currentMuscle.getMinimumFiberLength();
    % muscles(imus,1).maxpenn = currentMuscle.getMaximumPennationAngle();
    
    % get muscle path and DOFs this muscle crosses
    muspath = currentMuscle.getGeometryPath();
    PtSet = muspath.getPathPointSet();
    nPts = PtSet.getSize();
    
    %disp(muscles(imus,1).name);
    
    origin_segment = PtSet.get(0).getBody();
    insertion_segment = PtSet.get(nPts-1).getBody();
    
    current_segment = insertion_segment;
    %disp(['  current_segment: ', char(current_segment.getName())]);
    
    dof_count = 0;
    dof_list = {};
    seg_list = segments.Properties.RowNames;
    isegment = find(strcmp(seg_list,char(current_segment.getName)));
    osegment = find(strcmp(seg_list,char(origin_segment.getName)));
    dof_indexes = [];
    
    while (isegment ~= osegment)
        switch opensimVersion
            case opensimVersion33
                if ~current_segment.hasJoint()
                    error(['No muscle path found for muscle ',muscles(imus,1).name]);
                end
            case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
                %no accoring function for OpenSim 4 was found
        end
        
        % accumulate dof_indexes
        dof_indexes  = [dof_indexes,find(dofs_segment_list(isegment,:)==1)];
        
        switch opensimVersion
            case opensimVersion33
                current_joint = current_segment.getJoint();
                current_segment = current_joint.getParentBody();
                isegment = find(strcmp(seg_list,char(current_segment.getName)));
            case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
                %this is a dirty temporary hack, but no according function
                %was found for OpenSim 4
                isegment = find(strcmp(seg_list,joints.parent_segment(isegment-1)));
        end
    end
    
    dof_count = length(dof_indexes);
    dof_list = dofs.Properties.RowNames(dof_indexes); 
    muscles(imus,1).dof_count = dof_count;
    muscles(imus,1).dof_indexes = dof_indexes;
    muscles(imus,1).dof_names = dof_list';
    muscles(imus,1).weight =  muscles(imus,1).fmax * muscles(imus,1).lceopt;
end
model.muscles = struct2table(muscles); % convert structure to table

% Make table
muscles = struct2table(muscles);
% Set name as row name
muscles.Properties.RowNames = muscles.name;
muscles.name = [];
% Set property
model.muscles = muscles;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARKERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of markers
MarkerSet = Mod.getMarkerSet();
nMarkers = MarkerSet.getSize();

% loop through all markers
for imarkers = 1:nMarkers
    currentMarker = MarkerSet.get(imarkers-1);
    
    % basic properties
    markers(imarkers,1).name = char(currentMarker.getName());
    switch opensimVersion
        case opensimVersion33
            markers(imarkers,1).segment = char(currentMarker.getBodyName());
            markers(imarkers,1).fixed = currentMarker.getFixed();
            tmpVec3 = currentMarker.getOffset();
        case {opensimVersion40, opensimVersion41, opensimVersion43, opensimVersion45}
            markers(imarkers,1).segment = char(currentMarker.getParentFrame());
            markers(imarkers,1).fixed = currentMarker.get_fixed();
             tmpVec3 = currentMarker.get_location();
    end
    markers(imarkers,1).position = [tmpVec3.get(0) tmpVec3.get(1) tmpVec3.get(2)];
    
    % find index
    segment_names = segments.Properties.RowNames;
    for ibodys = 1:nSegments
        if strcmp(segment_names(ibodys), markers(imarkers,1).segment)
            markers(imarkers,1).segmentindex = ibodys;
        end
    end
end

% The default gait10dof18musc does not have a markerset, therefore:
if nMarkers
    % Make table
    markers = struct2table(markers);
    % Set name as row name
    markers.Properties.RowNames = markers.name;
    markers.name = [];
else
    markers = [];
end

% Set property
model.markers = markers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTACT POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get contact set
ContactSet = Mod.getContactGeometrySet();
nContacts = ContactSet.getSize();
warning('Gait2d_osim:ReadOsim', 'Default Values for the contact model are used. Default parameters borrowed from gait3d');

for icontact = 1:nContacts
    current_contact = ContactSet.get(icontact-1);
    % No switching between opensim versions here:
    % Contact spheres are only allowed with 4.x
    if strcmp(char(current_contact.getName()),'floor')
        continue;
    end
    contacts(icontact,1).name = char(current_contact.getName());

    % ToDo: this only allows for contacts at the feet, since getFrame() somhow doesnt work. 
    if strcmp(contacts(icontact,1).name(end),'r')
        contacts(icontact,1).segment = 'foot_r';
        contacts(icontact,1).segmentindex = 0;
    elseif strcmp(contacts(icontact,1).name(end),'l')
        contacts(icontact,1).segment = 'foot_l';
        contacts(icontact,1).segmentindex = 1;
    end
    % Defining the exact ground contact point where the sphere hits the ground
    location = current_contact.getLocation();
    radius = ContactSphere.safeDownCast(current_contact);
    radius = radius.getRadius();
    contacts(icontact,1).position = [location.get(0) location.get(1)-radius location.get(2)];
    
    % contact variables are Fx,Fy, xc,yc,xf,yf (in m)
    contacts(icontact,1).range = [-2 -0.5 -1 -0.5; ...
                         2   3   5  1.0];
    


end
%set Table
if nContacts
    contacts = struct2table(contacts);
    contacts.Properties.RowNames = contacts.name;
    contacts.name = [];
    model.CPs = contacts;
end

%%%% Save the .mat file for the model
save(strrep(osim_path, '.osim', '.mat'),'model','osim_sha256');

    %> @cond DO_NOT_DOCUMENT
    %======================================================================
    %> @brief Function to get coefficients for a specific fuction
    %>
    %> @details
    %> Different functions can be used for translation and rotation of a
    %> joint location. This function reads the translation and rotation
    %> polynomial coefficients from the function object f. The polynomial
    %> is of order 4. Optional, there can be a scaling of the function due to
    %> model scaling.
    %>
    %> @param   f      org.opensim.modeling.function: function describing
    %>                 the translation or rotation
    %> @param   scale  (optional) Double: Scaling factor
    %> @retval  coefs  Double vector: (1x5)
    %======================================================================
    function [coefs] = osimfunction2coefs(f, scale)
        
        if nargin < 2
           scale = 1; 
        end
        
        import org.opensim.modeling.*;
        coefs=[];
        
        % Constant factor
        fx=Constant.safeDownCast(f);
        if ~isempty(fx)
            coefs = [scale*fx.getValue() 0 0 0 0];
        end
        
        % Linear function
        fx=LinearFunction.safeDownCast(f);
        if ~isempty(fx)
            coefs = [scale*fx.getIntercept() scale*fx.getSlope() 0 0 0];
        end
        
        % Spline function
        fx=SimmSpline.safeDownCast(f);
        if ~isempty(fx)
            s = fx.getSize();
            
            knots=zeros(2,s);
            for i=1:s
                knots(1,i) = fx.getX(i-1);
                knots(2,i) = fx.getY(i-1);
            end
            
            coefs = polyfit(knots(1,:),scale*knots(2,:),4);
        end
        
        % Multiplier function for scaling 
        fx = MultiplierFunction.safeDownCast(f);
        if ~isempty(fx)
            scale = fx.getScale();
            % Call this function recursivly using the scale factor
            coefs = osimfunction2coefs(fx.getFunction(),scale);
        end
        
        if isempty(coefs)
            error('Function Type unknown');
        end
    end
    %> @endcond

end



