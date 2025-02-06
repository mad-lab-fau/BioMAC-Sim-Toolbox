%======================================================================
%> @file @TrackingData/inclineGRF.m
%> @brief TrackingData function to incline force plates
%> @details
%> Details: TrackingData::inclineGRF()
%>
%> @author Markus Gambietz
%> @date October, 2022
%======================================================================

% ======================================================================
%> @brief Function to change a treadmill's incline
%>
%> @param   obj     TrackingData class object which should be corrected
%> @param   incline Angle of treadmill incline, rad
% ======================================================================
function inclineGRF(obj, incline)

% Rotation matrix for x any incline
R = [cosd(incline) -sind(incline); sind(incline) cosd(incline)];

% right side
if ismember('GRF_x_r',obj.variables.name) && ismember('GRF_y_r',obj.variables.name)
    idxy = find(string(obj.variables.name) == "GRF_y_r");
    idxx = find(string(obj.variables.name) == "GRF_x_r");   
    for i = 1:obj.nSamples
        tmp = (R * [obj.variables.mean{idxx}(i);obj.variables.mean{idxy}(i)])';
        obj.variables.mean{idxx}(i)=tmp(1);
        obj.variables.mean{idxy}(i)=tmp(2);
    end
end

%left side
if ismember('GRF_x_l',obj.variables.name) && ismember('GRF_y_l',obj.variables.name)
    idxy = find(string(obj.variables.name) == "GRF_y_l");
    idxx = find(string(obj.variables.name) == "GRF_x_l");   
    for i = 1:obj.nSamples
        tmp = (R * [obj.variables.mean{idxx}(i);obj.variables.mean{idxy}(i)])';
        obj.variables.mean{idxx}(i)=tmp(1);
        obj.variables.mean{idxy}(i)=tmp(2);
    end
end
end