%======================================================================
%> @file @Gait2dc/getCoM.m
%> @brief Gait2dc function to calculate the center of mass
%> @details
%> Details: Model::getCoM()
%>
%> @author Anne Koelewijn
%> @date December, 2019
%======================================================================

%======================================================================
%> @brief Function to calculate the center of mass
%>
%> @details
%> It solves the center of mass (CoM) coordinates from the locations of
%> the individual centers of mass
%> 
%> @param   obj      Model class object
%> @param   x        Double vector: joint angles (1xnStates or nStatesx1)
%> @retval  CoM      Double vector: x, y location of center of mass
%======================================================================
function [CoM, dCoMdx] = getCoM(obj, x)

if size(x,1) == 1
    x = x';
end

part_mass = obj.parameter(12:15,2); %mass of each body part
CoM_loc = obj.parameter(12:15, [4 5]);

%Calculate global angles --> could be indiced automatically
trunk_angle = x(3);
RThigh_angle = x(3)+x(4);
LThigh_angle = x(3)+x(7);
RShank_angle = x(3)+x(4)+x(5);
LShank_angle = x(3)+x(7)+x(8);
RFoot_angle = x(3)+x(4)+x(5)+x(6);
LFoot_angle = x(3)+x(7)+x(8)+x(9);

angle_L = [LThigh_angle;LShank_angle;LFoot_angle];
angle_R = [trunk_angle;RThigh_angle;RShank_angle;RFoot_angle];

% CoMx = (sum(sin(angle_R).*CoM_loc.*part_mass) + sum(sin(angle_L).*CoM_loc(2:end).*part_mass(2:end)))./sum(bodymass);
CoM = (sum(cos(angle_R).*CoM_loc(:,2).*part_mass - sin(angle_R).*CoM_loc(:,1).*part_mass) + ...
    sum(cos(angle_L).*CoM_loc(2:end,2).*part_mass(2:end)- sin(angle_L).*CoM_loc(2:end,1).*part_mass(2:end)))./obj.bodymass;

if nargout == 2
    dCoMdx = spalloc(1,length(x),7);
    dCoMdangleR = -(CoM_loc(:,2).*part_mass.*sin(angle_R)+CoM_loc(:,1).*part_mass.*cos(angle_R))/obj.bodymass;
    dCoMdangleL = -(CoM_loc(2:end,2).*part_mass(2:end).*sin(angle_L) + CoM_loc(2:end,1).*part_mass(2:end).*cos(angle_L))./obj.bodymass;
    dCoMdx(3) = sum(dCoMdangleR)+sum(dCoMdangleL);
    dCoMdx(4) = sum(dCoMdangleR(2:end));
    dCoMdx(5) = sum(dCoMdangleR(3:end));
    dCoMdx(6) = sum(dCoMdangleR(4:end));
    dCoMdx(7) = sum(dCoMdangleL);
    dCoMdx(8) = sum(dCoMdangleL(2:end));
    dCoMdx(9) = sum(dCoMdangleL(3:end));
end