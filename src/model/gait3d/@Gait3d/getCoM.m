%======================================================================
%> @file @Gait3d/getCoM.m
%> @brief Gait3d function to calculate the center of pressure
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
%> It solves the center of pressure (CoM) coordinates from the locations of
%> the individual centers of mass
%>
%> 
%> @param   obj      Model class object
%> @param   x        Double vector: state of system
%> @retval  CoM      Double vector: x, y, z location of center of mass
%======================================================================
function CoM = getCoM(obj, x)

partmass = obj.segments.mass(2:obj.nSegments);
com_loc = obj.segments.mass_center;

fk = obj.getFkin(x(1:obj.nDofs));
fk = reshape(fk, 12, obj.nSegments-1);  
CoM_glo = zeros(size(com_loc,2),size(com_loc,1)-1);
for i=1:obj.nSegments-1
    O = fk(1:3,i);							% position of origin
    R = reshape(fk(4:12,i)',3,3)';			% rotation matrix
    
    CoM_glo(:,i) = O+R*com_loc(i+1,:)';
end

CoM = sum(CoM_glo.*partmass',2)/sum(partmass);





