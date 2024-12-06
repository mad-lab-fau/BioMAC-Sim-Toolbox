%======================================================================
%> @file @Gait2dc/readParameter.m
%> @brief Gait2dc function to read parameter matrix into tables
%> @details
%> Details: Gait2dc::readParameter()
%> 
%> @author Ton, Eva, Iris, Marlies
%> @date July, 2017
%======================================================================

%======================================================================
%> @brief Function to read parameter matrix into tables
%> @protected
%>
%> @details 
%> Function reads the parameters from the matrix Gait2dc.parameter and sets
%> Gait2dc.segments, Gait2dc.joints, Gait2dc.dofs, Gait2dc.muscles and Gait2dc.CPs. 
%> If you change sth here, scaleParameters() should be also checked!!!
%> Sets also the default range of Gait2dc.CPs.
%> Gait2dc.bodymass and Gait2dc.bodyheight are set according to the input of the
%> constructor (value in excelfile is ignored).
%>
%> @param   obj           Gait2dc class object
%======================================================================
function readParameter(obj)

% Get indices where a new section is starting
idxParamSections = obj.idxParamSections;
% 2 lines of header before first numeric value occurs
headerOffset = 2; 


% Set gravity, drag_coefficient, wind_speed and slope
obj.gravity          = obj.parameter(idxParamSections(1) +1, 2);
obj.drag_coefficient = obj.parameter(idxParamSections(1) +2, 2);
obj.wind_speed       = obj.parameter(idxParamSections(1) +3, 2);
obj.slope            = obj.parameter(idxParamSections(1) +4, 2);


% Bodymass and bodyheight is set according to the input of the
% constructor.

% Set segments:
% Hardcode name, index and name of parent joint
segments(1, 1).name         = 'pelvis';
segments(1, 1).iSection     = 1;
segments(1, 1).iRowInSection= 11;
segments(1, 1).parent_joint = '';

segments(2, 1).name = 'femur_r';
segments(2, 1).iSection     = 1;
segments(2, 1).iRowInSection= 12;
segments(2, 1).parent_joint = 'hip_r';

segments(3, 1).name = 'tibia_r';
segments(3, 1).iSection     = 1;
segments(3, 1).iRowInSection= 13;
segments(3, 1).parent_joint = 'knee_r';

segments(4, 1).name = 'foot_r';
segments(4, 1).iSection     = 1;
segments(4, 1).iRowInSection= 14;
segments(4, 1).parent_joint = 'ankle_r';

segments(5, 1).name = 'femur_l';
segments(5, 1).iSection     = 1;
segments(5, 1).iRowInSection= 12;
segments(5, 1).parent_joint = 'hip_l';

segments(6, 1).name = 'tibia_l';
segments(6, 1).iSection     = 1;
segments(6, 1).iRowInSection= 13;
segments(6, 1).parent_joint = 'knee_l';

segments(7, 1).name = 'foot_l';
segments(7, 1).iSection     = 1;
segments(7, 1).iRowInSection= 14;
segments(7, 1).parent_joint = 'ankle_l';

% Get mass, mass_center, inertia and length from excel file (loaded to
% obj.parameter)
for iSeg = 1 : length(segments)
    idx = idxParamSections(segments(iSeg, 1).iSection) + segments(iSeg, 1).iRowInSection;
    segments(iSeg, 1).mass         = obj.parameter(idx , 2);
    segments(iSeg, 1).mass_center  = obj.parameter(idx , 4:5);
    segments(iSeg, 1).inertia      = obj.parameter(idx , 3);
    segments(iSeg, 1).length       = obj.parameter(idx , 6);
end

% Make table
segments = struct2table(segments);
% Set name as row name
segments.Properties.RowNames = segments.name;
segments.name = [];
% Set object
obj.segments = segments;

% Set muscles:
% Hardcode name
muscles(1, 1).name = 'Iliopsoas_r';
muscles(2, 1).name = 'Glutei_r';
muscles(3, 1).name = 'Hamstrings_r';
muscles(4, 1).name = 'Rectus_r';
muscles(5, 1).name = 'Vasti_r';
muscles(6, 1).name = 'Gastroc_r';
muscles(7, 1).name = 'Soleus_r';
muscles(8, 1).name = 'TibialisAnt_r';
muscles(9, 1).name = 'Iliopsoas_l';
muscles(10, 1).name = 'Glutei_l';
muscles(11, 1).name = 'Hamstrings_l';
muscles(12, 1).name = 'Rectus_l';
muscles(13, 1).name = 'Vasti_l';
muscles(14, 1).name = 'Gastroc_l';
muscles(15, 1).name = 'Soleus_l';
muscles(16, 1).name = 'TibialisAnt_l';
[muscles.iSection] = deal(2);

% Get muscle parameters from excel file (loaded to obj.parameter) and save
% index
for iMus = 1 : length(muscles)
    muscles(iMus, 1).iRowInSection= 1 + iMus;
    idx = idxParamSections(muscles(iMus, 1).iSection) + muscles(iMus, 1).iRowInSection;
    muscles(iMus, 1).fmax         =  obj.parameter(idx , 2);
    muscles(iMus, 1).lceopt       =  obj.parameter(idx , 3);
    muscles(iMus, 1).width        =  obj.parameter(idx , 4);
    muscles(iMus, 1).peeSlack     =  obj.parameter(idx , 5);
    muscles(iMus, 1).seeSlack     =  obj.parameter(idx , 6);
    muscles(iMus, 1).L0           =  obj.parameter(idx , 7);
    muscles(iMus, 1).dRhip        =  obj.parameter(idx , 8);
    muscles(iMus, 1).dRknee       =  obj.parameter(idx , 9);
    muscles(iMus, 1).dRankle      =  obj.parameter(idx , 10);
    muscles(iMus, 1).dLhip        =  obj.parameter(idx , 11);
    muscles(iMus, 1).dLknee       =  obj.parameter(idx , 12);
    muscles(iMus, 1).dLank        =  obj.parameter(idx , 13);
    muscles(iMus, 1).kPEE         =  obj.parameter(idx , 14);
    muscles(iMus, 1).umax         =  obj.parameter(idx , 15);
    muscles(iMus, 1).vmax         =  obj.parameter(idx , 16);
    muscles(iMus, 1).tact         =  obj.parameter(idx , 17);
    muscles(iMus, 1).tdeact       =  obj.parameter(idx , 18);
    muscles(iMus, 1).gmax         =  obj.parameter(idx , 19);
    muscles(iMus, 1).arel         =  obj.parameter(idx , 20);
    muscles(iMus, 1).FT           =  obj.parameter(idx , 21);
end

% Make table
muscles = struct2table(muscles);
% Set name as row name
muscles.Properties.RowNames = muscles.name;
muscles.name = [];
% Set object
obj.muscles = muscles;

% Set CPs:
% Search in whole section 3 for CPs 
% (An empty line in between is not allowed)
ncontacts = 0;
for iRow = (idxParamSections(3) + 2) : (idxParamSections(4) -1) 
    if ~isnan(obj.parameter(iRow, 2)) % there is an entry
        ncontacts = ncontacts +1;
    end 
end   

% Get CP names
if ~isempty(obj.excelfile)
    % Names can be obtained from excel
    [~, names, ~] = xlsread(obj.excelfile);
    names = names((1 : ncontacts)+headerOffset+idxParamSections(3)+1, 1);
    names = strrep(names, 'R', ''); % Remove leading R and L. This is not nice, but needed for update_isymmetry
    names = strrep(names, 'L', '');
else
    % Names are not given => Number them consecutively
    names = cell(ncontacts, 1);
    ncontactsR = 0;
    ncontactsL = 0;
    for iCPs = 1 : ncontacts
        idx = idxParamSections(3) + 1 + iCPs;
        if obj.parameter(idx, 2) == 0 % right foot
            ncontactsR = ncontactsR +1;
            names{iCPs} = sprintf('CP_%02i', ncontactsR);
        elseif obj.parameter(idx, 2) == 1 % left foot
            ncontactsL = ncontactsL +1;
            names{iCPs} = sprintf('CP_%02i', ncontactsL);
        end
    end
end

% Get CP parameters
for iCPs = 1 : ncontacts
    CPs(iCPs, 1).iSection   = 3; 
    CPs(iCPs, 1).iRowInSection  = 1 + iCPs;
    idx = idxParamSections(CPs(iCPs, 1).iSection) + CPs(iCPs, 1).iRowInSection;
    if obj.parameter(idx, 2) == 0 % right foot
        CPs(iCPs, 1).name = [names{iCPs} '_r']; % This ending is needed
        CPs(iCPs, 1).segment = 'foot_r'; 
    elseif obj.parameter(idx, 2) == 1 % left foot
        CPs(iCPs, 1).name = [names{iCPs} '_l']; % This ending is needed
        CPs(iCPs, 1).segment = 'foot_l'; 
    else
        error(['Error in excel parameter file ' obj.excelfile ...
            ': Only segment 0 or segment 1 are valid entries for the ground contact points']);
    end
    CPs(iCPs, 1).segmentID  = obj.parameter(idx, 2);
    CPs(iCPs, 1).x  = obj.parameter(idx, 3);
    CPs(iCPs, 1).y  = obj.parameter(idx, 4);
    CPs(iCPs, 1).k1 = obj.parameter(idx, 5);
    CPs(iCPs, 1).k2 = obj.parameter(idx, 6);
    CPs(iCPs, 1).a  = obj.parameter(idx, 7);
    CPs(iCPs, 1).c  = obj.parameter(idx, 8);
    CPs(iCPs, 1).b  = obj.parameter(idx, 9);
    % Set default range for all CPs
    CPs(iCPs, 1).range = [-5 -0.5 -10 -0.5; 7 5 10 10];

end

% Make table
CPs = struct2table(CPs);
% Set name as row name
CPs.Properties.RowNames = CPs.name;
CPs.name = [];
% Set object
obj.CPs = CPs;


% Set joints:
% Hardcode name, index, parent_segment and segment
joints(1, 1).name = 'hip_r';
joints(1, 1).iSection = 4;
joints(1, 1).iRowInSection = 2;
joints(1, 1).parent_segment = 'pelvis';
joints(1, 1).segment = 'femur_r';

joints(2, 1).name = 'knee_r';
joints(2, 1).iSection = 4;
joints(2, 1).iRowInSection = 3;
joints(2, 1).parent_segment = 'femur_r';
joints(2, 1).segment = 'tibia_r';

joints(3, 1).name = 'ankle_r';
joints(3, 1).iSection = 4;
joints(3, 1).iRowInSection = 4;
joints(3, 1).parent_segment = 'tibia_r';
joints(3, 1).segment = 'foot_r';

joints(4, 1).name = 'hip_l';
joints(4, 1).iSection = 4;
joints(4, 1).iRowInSection = 5;
joints(4, 1).parent_segment = 'pelvis';
joints(4, 1).segment = 'femur_l';

joints(5, 1).name = 'knee_l';
joints(5, 1).iSection = 4;
joints(5, 1).iRowInSection = 6;
joints(5, 1).parent_segment = 'femur_l';
joints(5, 1).segment = 'tibia_l';

joints(6, 1).name = 'ankle_l';
joints(6, 1).iSection = 4;
joints(6, 1).iRowInSection = 7;
joints(6, 1).parent_segment = 'tibia_l';
joints(6, 1).segment = 'foot_l';

% Get joint parameters from excel file (loaded to obj.parameter) and save
% index
for iJoin = 1 : length(joints)
    idx = idxParamSections(joints(iJoin, 1).iSection) + joints(iJoin, 1).iRowInSection;
    joints(iJoin, 1).jointK2      =  obj.parameter(idx, 4);
    joints(iJoin, 1).jointD       =  obj.parameter(idx, 5);
    joints(iJoin, 1).jointB       =  obj.parameter(idx, 6);
    joints(iJoin, 1).jointK1      =  obj.parameter(idx, 7);
    joints(iJoin, 1).jointPhi0    =  obj.parameter(idx, 8);
end

% Make table
joints = struct2table(joints);
% Set name as row name
joints.Properties.RowNames = joints.name;
joints.name = [];
% Set object
obj.joints = joints;

% Set dofs:
% Hardcode name and get range from excel file (loaded to obj.parameter)
dofs(1, 1).name = 'pelvis_tx';
dofs(1, 1).joint = 'ground_pelvis';
dofs(1, 1).range = [-5 7]; % no range available in excel file

dofs(2, 1).name = 'pelvis_ty';
dofs(2, 1).joint = 'ground_pelvis';
dofs(2, 1).range = [0.5 2]; % no range available in excel file

dofs(3, 1).name = 'pelvis_tilt';
dofs(3, 1).joint = 'ground_pelvis';
dofs(3, 1).range = [-pi pi]; % no range available in excel file

dofs(4, 1).name = 'hip_flexion_r';
dofs(4, 1).joint = 'hip_r';
idx = idxParamSections(obj.joints{dofs(4, 1).joint, 'iSection'}) + obj.joints{dofs(4, 1).joint, 'iRowInSection'};
dofs(4, 1).range = obj.parameter(idx, 2:3)/180*pi;

dofs(5, 1).name = 'knee_angle_r';
dofs(5, 1).joint = 'knee_r';
idx = idxParamSections(obj.joints{dofs(5, 1).joint, 'iSection'}) + obj.joints{dofs(5, 1).joint, 'iRowInSection'};
dofs(5, 1).range = obj.parameter(idx, 2:3)/180*pi;

dofs(6, 1).name = 'ankle_angle_r';
dofs(6, 1).joint = 'ankle_r';
idx = idxParamSections(obj.joints{dofs(6, 1).joint, 'iSection'}) + obj.joints{dofs(6, 1).joint, 'iRowInSection'};
dofs(6, 1).range = obj.parameter(idx, 2:3)/180*pi;

dofs(7, 1).name = 'hip_flexion_l';
dofs(7, 1).joint = 'hip_l';
idx = idxParamSections(obj.joints{dofs(7, 1).joint, 'iSection'}) + obj.joints{dofs(7, 1).joint, 'iRowInSection'};
dofs(7, 1).range = obj.parameter(idx, 2:3)/180*pi;

dofs(8, 1).name = 'knee_angle_l';
dofs(8, 1).joint = 'knee_l';
idx = idxParamSections(obj.joints{dofs(8, 1).joint, 'iSection'}) + obj.joints{dofs(8, 1).joint, 'iRowInSection'};
dofs(8, 1).range = obj.parameter(idx, 2:3)/180*pi;

dofs(9, 1).name = 'ankle_angle_l';
dofs(9, 1).joint = 'ankle_l';
idx = idxParamSections(obj.joints{dofs(9, 1).joint, 'iSection'}) + obj.joints{dofs(9, 1).joint, 'iRowInSection'};
dofs(9, 1).range = obj.parameter(idx, 2:3)/180*pi;

% Make table
dofs = struct2table(dofs);
% Set name as row name
dofs.Properties.RowNames = dofs.name;
dofs.name = [];
% Set object
obj.dofs = dofs;



% Set strainEnergyTerms:
obj.strainEnergyTerms = obj.parameter((1:20)+idxParamSections(5)+1, 2:8);



end

