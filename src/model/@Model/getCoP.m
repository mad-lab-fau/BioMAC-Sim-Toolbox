%======================================================================
%> @file @Model/getCoP.m
%> @brief Model function to calculate the center of pressure
%> @details
%> Details: Model::getCoP()
%>
%> @author Marlies Nitschke
%> @date February, 2019
%======================================================================

%======================================================================
%> @brief Function to calculate the center of pressure
%>
%> @details
%> It solves the center of pressure (COP) coordinates from: COP x F + Ty = M, where COP
%> is a point in the XZ plane and Ty is the "free moment" on the Y axis.
%> expand the equation:
%>   1. COPy * Fz - COPz * Fy + 0  = Mx
%>   2. COPz * Fx - COPx * Fz + Ty = My
%>   3. COPx * Fy - COPy * Fx + 0  = Mz
%>
%> Use COPy = 0, to get COPx = Mz / Fy  and COPz = -Mx / Fy
%>
%> For 2D, COPz will be 0.
%> 
%> @param   obj      Model class object
%> @param   grf      Double vector: Ground contact vector containing Fx, Fy, Fz, Mx, My, Mz (right), and 
%>                   Fx, Fy, Fz, Mx, My, Mz (left) (12)
%> @retval  CoP_r    Double vector: Center of pressure in x, y and z direction for right foot (3)
%> @retval  CoP_l    Double vector: Center of pressure in x, y and z direction for left foot (3)
%======================================================================
function [CoP_r, CoP_l] = getCoP(obj, grf)

% Get forces and moments of right and left foot
F_r = grf(1:3);
M_r = grf(4:6);
F_l = grf(7:9);
M_l = grf(10:12);

% Compute CoP of right foot
CoP_r(1) = M_r(3)/F_r(2);
CoP_r(2) = 0;
CoP_r(3) = -M_r(1)/F_r(2);

% Compute CoP of left foot
CoP_l(1) = M_l(3)/F_l(2);
CoP_l(2) = 0;
CoP_l(3) = -M_l(1)/F_l(2);
        
end



