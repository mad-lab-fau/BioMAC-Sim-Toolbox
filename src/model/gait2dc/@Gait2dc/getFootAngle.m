%======================================================================
%> @file @Gait2dc/getFootAngle.m
%> @brief Gait2dc function to calculate the angle between foot and ground
%> @details
%> Details: Gait2dc::getFootAngle()
%>
%> @author Marlies Nitschke
%> @date January, 2019
%======================================================================

%======================================================================
%> @brief Function to calculate the angle between foot and ground
%>
%> @details
%> Function uses gait2dc('Stick', x(:, iTime)) to get the position of the
%> foot segment.
%> 
%> This implementation is only suited for our 2D model containing two
%> contact points per foot (heel and toe)!
%>
%> @param   obj      Gait2dc class object
%> @param   x        Double matrice: State vector of model for n time points (Gait2dc.nStates x n)
%> @retval  angle_r  Double vector: Angle between right foot and ground in degree. Angles ranges from -180 to 180 deg. (n x 1)
%> @retval  angle_l  Double vector: Angle between left foot and ground in degree. Angles ranges from -180 to 180 deg. (n x 1)
%======================================================================
function [angle_r, angle_l] = getFootAngle(obj, x)

    % Check if the function can be used
    assert(obj.nCPs == 4, 'The function is implemented using 2 contact points for each foot. You have %i contact points in total.', obj.nCPs);

    % Preallocate output
    angle_r = nan(size(x, 2), 1);
    angle_l = nan(size(x, 2), 1);

    % Do it for all nodes
    for iNode = 1 : size(x,2)

        % Get stick figure
        [R, L] = gait2dc('Stick', x(:, iNode));

        % Get points defining the sole of the feet
        foot_heel_r = R(end-2, :);
        foot_toe_r  = R(end-1, :);
        foot_heel_l = L(end-2, :);
        foot_toe_l  = L(end-1, :);

        % Get angles
        angle_r(iNode) = calculateFootAngle(foot_heel_r, foot_toe_r);
        angle_l(iNode) = calculateFootAngle(foot_heel_l, foot_toe_l);

    end

end


%> @cond DO_NOT_DOCUMENT
%======================================================================
%> @brief Function to calculate the angle between foot and ground
%>
%> @param  cord_heel  Double vector: X and Y coordinate of the heel (1x2) or (2x1)
%> @param  cord_toe   Double vector: X and Y coordinate of the toe (1x2) or (2x1)
%> @retval angle      Double: Foot angle 
%======================================================================
function angle = calculateFootAngle(cord_heel, cord_toe)

    % Get differences of heel and toe
    delta_y = (cord_toe(2)-cord_heel(2));
    delta_x = (cord_toe(1)-cord_heel(1));
    
    % Get angle
    if     delta_x > 0  && delta_y >= 0  % 1. Quadrant
        angle = atand(delta_y / delta_x);
    elseif delta_x < 0  && delta_y >= 0  % 2. Quadrant
        angle = atand(delta_y / delta_x) + 180;
    elseif delta_x < 0  && delta_y <= 0  % 4. Quadrant
        angle = atand(delta_y / delta_x) - 180;
    elseif delta_x > 0  && delta_y <= 0  % 3. Quadrant
        angle = atand(delta_y / delta_x);
    elseif delta_x == 0 && delta_y > 0   % Vertical between 1. Quadrant and 2. Quadrant
        angle = 90;
    elseif delta_x == 0 && delta_y < 0   % Vertical between 3. Quadrant and 4. Quadrant
        angle = -90;
    elseif delta_x == 0 && delta_y == 0  % Horizontal
        angle = 0;
    else
        error('Gait2dc:calculateFootAngle', ...
            'There is an unknown case during the computation for heel (%d, %d) and toe (%d, %d)', ...
            coord_heel(1), coord_heel(2), coord_toe(1), coord_toe(2));
    end

end
%> @endcond