% ======================================================================
%> @file @Gait2dc/simuMarker.m
%> @brief Gait2dc function to simulates markers placed on 2D model
%> @details
%> Details: Gait2dc::simuMarker()
%>
%> @author Marlies
%> @date June, 2021
% ======================================================================

%======================================================================
%> @brief Matlab function to simulates markers placed on 2D model
%> @public
%>
%> @details
%> Do not change "obj" or "variables" in this function!  This will allow
%> Matlab to "pass by reference" and avoid function call overhead.
%>
%> Position has to given as X, Y (and Z) position of marker in segment coordinate 
%> system in m.
%> Direction has to be given in axis definition, e.g [1, 0, 0].
%> For position and direction, we will use only the first two entries for
%> the 2D model.
%>
%> @param obj         Gait2dc class object
%> @param variables   Table: Variable table containing the marker data
%>                    with at least the columns: type, name, segment, position, direction
%> @param q           Double array: Generalized coordinates (Gait2dc.nDofs x 1)
%> @param idxSegment  (optional) Double array: Indices of segments fitting to respective markers (variables.nVars.all x 1)
%>                    => It is faster to pass them than to recalculate each time.
%> @param plocalAll   (optional) Double matrix: Local position of the markers saved in variables (variables.nVars.all x 2(or3))
%>                    => It is faster to pass them than to extract them each time.
%> @param dlocalAll   (optional) Double matrix: Local direction of the markers saved in variables (variables.nVars.all x 2(or3))
%>                    => It is faster to pass them than to extract them each time.
%>
%> @retval m          Simulated marker signals in meter (variables.nVars.all x 1)
%> @retval dm_dq      Sparse Jacobian dm/dq (variables.nVars.all x Gait2dc.nDofs)
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
x = zeros(obj.nStates, 1); 
x(obj.hidxStates.q) = q; % gait2dc.c wants x and not q
if nargout == 1
    FK = obj.getFkin(x);
elseif nargout == 2
    [FK, dFKdq] = obj.getFkin(x);
end

% Get all local positions and directions of the markers 
if nargin < 5
    plocalAll = variables.position(:,1:2);
end
if nargin < 6
    dlocalAll = variables.direction(:,1:2);
end

% Simulate marker signal
for iVar = 1 : nVars
    
    % Local position and direction of the marker on this segment   
    plocal = plocalAll(iVar, :);
    dlocal = dlocalAll(iVar, :);
    
    % Position and orientation of the segment
    iSegment = idxSegment(iVar);
    idxp = (iSegment-1)*6 + (1:2);	% indices of position in FK
    idxR = (iSegment-1)*6 + (3:6);	% indices of orientation in FK
    
    % Marker position: m = dlocal * (p + R*plocal)
    % (This implementation was faster compared compared to reshaping of FK or to compute it as vector and multiplying later with dlocal.)
    % (The brackets around (R*plocal) showed in an example simulation to be important. The number of iterations and the result changed otherwise.)
    if dlocal(1)
        m(iVar) = FK(idxp(1)) + (FK(idxR(1), :) * plocal(1) + FK(idxR(2), :) * plocal(2));
    elseif dlocal(2)
        m(iVar) = FK(idxp(2)) + (FK(idxR(3), :) * plocal(1) + FK(idxR(4), :) * plocal(2));
    end
    
    % Jacobian
    if nargout == 2
        % dm_dq = dlocal * (dp_dq + dR_dq * plocal)
        % (This implementation was faster compared to reshaping dRdq and iterating over the DoFs or to compute it as vector and multiplying later with dlocal.)
        % (The brackets around (R*plocal) showed in an example simulation to be important. The number of iterations and the result changed otherwise.)
        if dlocal(1)
            dm_dq(iVar, :) = dFKdq(idxp(1), :) + (dFKdq(idxR(1), :) * plocal(1) + dFKdq(idxR(2), :) * plocal(2));
        elseif dlocal(2)
            dm_dq(iVar, :) = dFKdq(idxp(2), :) + (dFKdq(idxR(3), :) * plocal(1) + dFKdq(idxR(4), :) * plocal(2));
        end
    end
    
end

end
