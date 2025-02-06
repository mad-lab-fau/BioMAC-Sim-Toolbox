%======================================================================
%> @file @Gait3d/getDynamicsCPOffset.m
%> @brief Gait3d computing implicit differential equation for 3D musculoskeletal model
%> which allows an adaption of the CP position
%> @details
%> Details: Gait3d::getDynamicsCPOffset()
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================


%======================================================================
%> @brief Function computing implicit differential equation for 3D musculoskeletal model
%> which allows an adaption of the CP position
%>
%> @details
%> I calls Gait3d.getDynamics() and then overrights the entries to
%> allow apply a vertical shift to the CP position.
%>
%> @param   obj       Gait3d_CPOffset class object
%> @param   x         Double array: State of the model (Gait3d.nStates x 1)
%> @param   xdot      Double array: State derivatives (Gait3d.nStates x 1)
%> @param   u         Double array: Controls of the model (Gait3d.nControls x 1)
%> @param   CPYOffset Double: Vertical offset of original CP positions applied to all CPs
%>
%> @retval  f            Double array: Dynamic residuals (Gait3d.nConstraints x 1)
%> @retval	dfdx	     (optional) Double matrix: Transpose of Jacobian matrix df/dx 		(Gait3d.nStates x Gait3d.nConstraints)
%> @retval	dfdxdot	     (optional) Double matrix: Transpose of Jacobian matrix df/dxdot 	(Gait3d.nStates x Gait3d.nConstraints)
%> @retval	dfdu	     (optional) Double matrix: Transpose of Jacobian matrix df/du 		(Gait3d.nControls x Gait3d.nConstraints)
%> @retval  dfdCPYOffset (optional) Double matrix: Transpose of Jacobian matrix df/dCPYOffset (1 x Gait3d.nConstraints)
%======================================================================
function [f, dfdx, dfdxdot, dfdu, dfdCPYOffset] = getDynamicsCPOffset(obj,x,xdot,u, CPYOffset)

% Call standard getDynamics function
if nargout > 3
    [f, dfdx, dfdxdot, dfdu] = obj.getDynamics(x,xdot,u);
elseif nargout > 1
    [f, dfdx, dfdxdot] = obj.getDynamics(x,xdot,u);
else
    f = obj.getDynamics(x,xdot,u);
end

% Adapt constraints of the ground contact which are influenced
% by the position offset
idxCPEqu2 = find(strcmp(obj.constraints.type, 'CP') & strcmp(obj.constraints.equation, 'equality 2'));
f(idxCPEqu2) = f(idxCPEqu2) - CPYOffset;
dfdCPYOffset = zeros(1, obj.nConstraints);
dfdCPYOffset(idxCPEqu2) = -1;

end