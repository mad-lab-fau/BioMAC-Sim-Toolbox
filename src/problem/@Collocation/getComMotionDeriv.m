%======================================================================
%> @file getComMotionDeriv.m
%> @brief Collocation function to compuate center of mass motion
%> @details
%> Details: Collocation::getComMotionDeriv()
%>
%> @author Anne Koelewijn
%> @date April, 2020
%======================================================================

%======================================================================
%> @brief Function to compute center of mass motion
%>
%> @param obj           Collocation class object
%> @param X             Double array: State vector containing at least speed of the movement
%======================================================================
function dCom_motdx = getComMotionDeriv(obj,X)
% com_mot = 0;
nNodesDur = obj.nNodesDur;
dCom_motdx = spalloc(1,length(X),7*(nNodesDur-1));
for iNode=1:nNodesDur-1
    ix1 = obj.idx.states(:,iNode);
    ix2 = obj.idx.states(:,iNode+1);
    x1 = X(ix1);
    x2 = X(ix2);
    [CoM1,dcomdx1] = obj.model.getCoM(x1);
    [CoM2,dcomdx2] = obj.model.getCoM(x2);
%     com_mot = com_mot+(CoM1-CoM2)^2;%(CoM1/obj.model.bodyheight-0.6)^2;%

    dCom_motdx(1,ix1) = dCom_motdx(1,ix1) + 2*(CoM1-CoM2)*dcomdx1; %-2*(CoM1/obj.model.bodyheight-0.6)*dcomdx1/obj.model.bodyheight;%
    dCom_motdx(1,ix2) = dCom_motdx(1,ix2) - 2*(CoM1-CoM2)*dcomdx2;
end   
