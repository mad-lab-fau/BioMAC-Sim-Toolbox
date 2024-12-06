% ======================================================================
%> @file @Gait3d/simuMarker.m
%> @brief Gait3d function to simulates markers placed on 2D model
%> @details
%> Details: Gait3d::simuMarker()
%>
%> @author Marlies
%> @date June, 2021
% ======================================================================

%======================================================================
%> @brief Matlab function to simulates markers placed on 3D model
%> @public
%>
%> @details
%> Do not change "obj" or "variables" in this function!  This will allow
%> Matlab to "pass by reference" and avoid function call overhead.
%>
%> Position has to given as X, Y, and Z position of marker in segment coordinate 
%> system in m.
%> Direction has to be given in axis definition, e.g [1, 0, 0].
%>
%> @param obj         Gait3d class object
%> @param variables   Table: Variable table containing the marker data
%>                    with at least the columns: type, name, segment, position, direction
%> @param q           Double array: Generalized coordinates (Gait3d.nDofs x 1)
%> @param idxSegment  (optional) Double array: Indices of segments fitting to respective markers (variables.nVars.all x 1)
%>                    => It is faster to pass them than to recalculate each time.
%> @param plocalAll   (optional) Double matrix: Local position of the markers saved in variables (variables.nVars.all x 3)
%>                    => It is faster to pass them than to extract them each time.
%> @param dlocalAll   (optional) Double matrix: Local direction of the markers saved in variables (variables.nVars.all x 3)
%>                    => It is faster to pass them than to extract them each time.
%>
%> @retval m          Simulated marker signals in meter (variables.nVars.all x 1)
%> @retval dm_dq      Sparse Jacobian dm/dq (variables.nVars.all x Gait3d.nDofs)
%======================================================================
function [m, dm_dq] = simuMarker(obj, variables, q, idxSegment, plocalAll, dlocalAll)

% Get segment indices and check if all segments sepecified in the variables table exist
if nargin < 4
    [~, idxSegment] = ismember(variables.segment, obj.segments.Properties.RowNames);
end
if any(idxSegment==0)
    nonExistingSegments = unique(variables.segment(idxSegment==0));
    error('The following segments are not defined in the model: %s',  sprintf('%s ', nonExistingSegments{:}));
end

% Total number of variables
nVars = length(idxSegment); % much faster then size(variables, 1) or height(variables)

% Initialize the outputs
m = zeros(nVars, 1);
if nargout == 2
    nDofs = obj.nDofs;
    dm_dq = zeros(nVars, nDofs);
end
if nVars == 0; return; end

% Call the MEX function to compute the forward kinematics
if nargout == 1
    FK = obj.getFkin(q);
elseif nargout == 2
    [FK, dFKdq] = obj.getFkin(q);
end

% Get all local positions and directions of the markers 
% (According to the Matlab Profiler it is faster to do it this way instead
% of doing it in the for-loop.)
if nargin < 5
    plocalAll = variables.position;
end
if nargin < 6
    dlocalAll = variables.direction;
end

% Simulate marker signal
for iVar = 1 : nVars
    
    % Local position and direction of the marker on this segment   
    plocal = plocalAll(iVar, :);
    dlocal = dlocalAll(iVar, :);
    
    % Position and orientation of the segment
    iSegment = idxSegment(iVar);
    idxp = (iSegment-2)*12 + (1:3);	 % indices of position in FK
    idxR = (iSegment-2)*12 + (4:12); % indices of orientation in FK

    % Marker position: m = dlocal * (p + R*plocal)
    % (This implementation was faster compared to reshaping R and iterating over the DoFs or to compute it as vector and multiplying later with dlocal. -> tested for 2D only)
    % (The brackets around (R*plocal) showed in an example simulation to be important. The number of iterations and the result changed otherwise. -> tested for 2D only)
    if dlocal(1)
        m(iVar) = FK(idxp(1)) + (FK(idxR(1), :) * plocal(1) + FK(idxR(2), :) * plocal(2) + FK(idxR(3), :) * plocal(3));
    elseif dlocal(2)
        m(iVar) = FK(idxp(2)) + (FK(idxR(4), :) * plocal(1) + FK(idxR(5), :) * plocal(2) + FK(idxR(6), :) * plocal(3));
    elseif dlocal(3)
        m(iVar) = FK(idxp(3)) + (FK(idxR(7), :) * plocal(1) + FK(idxR(8), :) * plocal(2) + FK(idxR(9), :) * plocal(3));
    end
    
    % Jacobian
    if nargout == 2
        % dm_dq = dlocal * (dp_dq + dR_dq * plocal)
        % (This implementation was faster compared to reshaping dRdq and iterating over the DoFsor to compute it as vector and multiplying later with dlocal. -> tested for 2D only)
        % (The brackets around (R*plocal) showed in an example simulation to be important. The number of iterations and the result changed otherwise. -> tested for 2D only)
        if dlocal(1)
            dm_dq(iVar, :) = dFKdq(idxp(1), :) + (dFKdq(idxR(1), :) * plocal(1) + dFKdq(idxR(2), :) * plocal(2) + dFKdq(idxR(3), :) * plocal(3));
        elseif dlocal(2)
            dm_dq(iVar, :) = dFKdq(idxp(2), :) + (dFKdq(idxR(4), :) * plocal(1) + dFKdq(idxR(5), :) * plocal(2) + dFKdq(idxR(6), :) * plocal(3));
        elseif dlocal(3)
            dm_dq(iVar, :) = dFKdq(idxp(3), :) + (dFKdq(idxR(7), :) * plocal(1) + dFKdq(idxR(8), :) * plocal(2) + dFKdq(idxR(9), :) * plocal(3));
        end

    end
    
end

end
