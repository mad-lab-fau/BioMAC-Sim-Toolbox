%======================================================================
%> @file getComMotion.m
%> @brief Collocation function to compuate center of mass motion
%> @details
%> Details: Collocation::getComMotion()
%>
%> @author Anne Koelewijn
%> @date April, 2020
%======================================================================

%======================================================================
%> @brief Computes center of mass motion
%> @param obj           Collocation class object
%> @param X             Double array: State vector containing at least speed of the movement
%======================================================================
function com_mot = getComMotion(obj,X)
com_mot = 0;
for iNode=1:obj.nNodes
    x1 = X(obj.idx.states(:,iNode));
    x2 = X(obj.idx.states(:,iNode+1));
    CoM1 = obj.model.getCoM(x1);
    CoM2 = obj.model.getCoM(x2);
    com_mot = com_mot+(CoM1-CoM2)^2;%(CoM1/obj.model.bodyheight-0.6)^2;%
end

