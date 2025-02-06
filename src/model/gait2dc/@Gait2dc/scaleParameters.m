%======================================================================
%> @file @Gait2dc/scaleParameters.m
%> @brief Gait2dc function to scale the parameters defined in the excelfile
%> @details
%> Details: Gait2dc::scaleParameters()
%> 
%> @author Ton, Eva, Iris, Ann-Kristin, Marlies
%> @date July, 2017
%======================================================================

%======================================================================
%> @brief Function to scale the parameters defined in the excelfile
%> @protected
%>
%> @details 
%> Uses bodymass and height according to Winter 2005. Bodymass and
%> heigt have to be set during object initialization. Changes later changes
%> are not possible. Please create a new object in this case.
%>
%> @param   obj           Gait2dc class object
%======================================================================
function scaleParameters(obj)

% 1) Scale segment parameters
segments = obj.segments;

% Pelvis:
% Scale parameters of the hat (=pelvis)
segments('pelvis', 'mass')    = table( 0.678 * obj.bodymass );
segments('pelvis', 'length')  = table( (0.81 - 0.53) * obj.bodyheight );
segments('pelvis', 'inertia') = table( segments{'pelvis', 'mass'} * (0.496 * segments{'pelvis', 'length'})^2 );
segments('pelvis', 'mass_center') = table( [0, 0.626 * segments{'pelvis', 'length'}] );


% Thigh:
% Scale parameters of the right thigh (=femur_r)
segments('femur_r', 'mass')    = table( 0.1 * obj.bodymass );
segments('femur_r', 'length')  = table( (0.53 - 0.285) * obj.bodyheight );
segments('femur_r', 'inertia') = table( segments{'femur_r', 'mass'} * (0.323 * segments{'femur_r', 'length'})^2 );
segments('femur_r', 'mass_center') = table( [0, -0.433 * segments{'femur_r', 'length'}] );

% Copy scaled parameters to left thigh (=femur_l)
segments('femur_l', 'mass')    = segments('femur_r', 'mass');
segments('femur_l', 'length')  = segments('femur_r', 'length');
segments('femur_l', 'inertia') = segments('femur_r', 'inertia');
segments('femur_l', 'mass_center') =segments('femur_r', 'mass_center');


% Shank:
% Scale parameters of the right shank (=tibia_r)
segments('tibia_r', 'mass')    = table( 0.0465 * obj.bodymass );
segments('tibia_r', 'length')  = table( (0.285 - 0.039) * obj.bodyheight );
segments('tibia_r', 'inertia') = table( segments{'tibia_r', 'mass'} * (0.302 * segments{'tibia_r', 'length'})^2 );
segments('tibia_r', 'mass_center') = table( [0, -0.433 * segments{'tibia_r', 'length'}] );

% Copy scaled parameters to left shank (=tibia_l)
segments('tibia_l', 'mass')    = segments('tibia_r', 'mass');
segments('tibia_l', 'length')  = segments('tibia_r', 'length');
segments('tibia_l', 'inertia') = segments('tibia_r', 'inertia');
segments('tibia_l', 'mass_center') =segments('tibia_r', 'mass_center');


% Foot:
% Scale parameters of the right foot (=foot_r)
segments('foot_r', 'mass')    = table( 0.0145 * obj.bodymass );
segments('foot_r', 'length')  = table( 0.152 * obj.bodyheight );
segments('foot_r', 'inertia') = table( segments{'foot_r', 'mass'} * (0.475 * segments{'foot_r', 'length'})^2 );
segments('foot_r', 'mass_center') = table( [0.5 * segments{'foot_r', 'length'} - Gait2dc.HEELDEFAULTOFFSET,-0.5 * Gait2dc.FOOTSOLEOFFSET] );

% Copy scaled parameters to left foot (=foot_l)
segments('foot_l', 'mass')    = segments('foot_r', 'mass');
segments('foot_l', 'length')  = segments('foot_r', 'length');
segments('foot_l', 'inertia') = segments('foot_r', 'inertia');
segments('foot_l', 'mass_center') =segments('foot_r', 'mass_center');

obj.segments = segments;

% 2) Scale CPs
CPs = obj.CPs;

% Get default body height
mc = ?Gait2dc;
defaultBodyheight = mc.PropertyList(ismember({mc.PropertyList.Name}, 'bodyheight')).DefaultValue; % in m

if ~isa(obj, 'Gait2dc_CPmov')
    % Scale heel CPs
    CPs('heel_r', 'x') = table( -Gait2dc.HEELDEFAULTOFFSET / defaultBodyheight * obj.bodyheight );
    CPs('heel_l', 'x') = table( -Gait2dc.HEELDEFAULTOFFSET / defaultBodyheight * obj.bodyheight );

    % Scale toe CPs
    CPs('toe_r', 'x')  = table( segments{'foot_r', 'length'} + CPs{'heel_r', 'x'} - Gait2dc.TOEOFFSET );
    CPs('toe_l', 'x')  = table( segments{'foot_l', 'length'} + CPs{'heel_l', 'x'} - Gait2dc.TOEOFFSET );
end

obj.CPs = CPs;


% 3) Scale foot length
foot = obj.foot;

% Scale heel
foot{'heel_r', 'x'} = -Gait2dc.HEELDEFAULTOFFSET / defaultBodyheight * obj.bodyheight;
foot{'heel_l', 'x'} = -Gait2dc.HEELDEFAULTOFFSET / defaultBodyheight * obj.bodyheight;

% Scale toe
foot{'toe_r', 'x'} = segments{'foot_r', 'length'} + foot{'heel_r', 'x'};
foot{'toe_l', 'x'} = segments{'foot_l', 'length'} + foot{'heel_l', 'x'};

obj.foot = foot;


end