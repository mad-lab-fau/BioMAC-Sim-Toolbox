%======================================================================
%> @file getStepLength.m
%> @brief Collocation function to compute step length
%> @details
%> Details: Collocation::getStepLength()
%>
%> @author Anne Koelewijn
%> @date May, 2020
%======================================================================

%======================================================================
%> @brief Computes step length
%> @param obj           Collocation class object
%> @param X             Double array: State vector containing at least speed of the movement
%> @param name          Optional: left or right
%======================================================================

function SL = getStepLength(obj,X,name)

if nargin < 3
    if ~obj.isSymmetric
        error('Please provide left or right')
    end
else
    if ~(strcmpi(name, 'left') || strcmpi(name, 'right'))
        error('Incorrect input for name, please specify left or right')
    end
end
x = X(obj.idx.states);

len_GRF = length(obj.model.GRFNAMES);
GRF = zeros(len_GRF,obj.nNodes);
for iNode = 1:obj.nNodes %accumulate tracking variables for nNodes
    [GRF(:,iNode), ~] = obj.model.getGRF(X(obj.idx.states(:,iNode)));
end

if strcmp(name, 'left')
    %start at right leg, end at left leg
    ind_grf_start = find(strcmp(obj.model.GRFNAMES, 'rightFy'));
    ind_grf_end = find(strcmp(obj.model.GRFNAMES, 'leftFy'));
    %Now find state index of heel marker
    if isa(obj.model, 'Gait3d')
        ind_state_start = find(strcmp(obj.model.states.name, 'CPM_r') & strcmp(obj.model.states.type,'xc'));
        ind_state_end = find(strcmp(obj.model.states.name, 'CPM_l') & strcmp(obj.model.states.type,'xc'));
    elseif isa(obj.model, 'Gait2dc')
        ind_state_start = find(strcmp(obj.model.states.name, 'heel_r') & strcmp(obj.model.states.type,'xc'));
        ind_state_end = find(strcmp(obj.model.states.name, 'heel_l') & strcmp(obj.model.states.type,'xc'));
    else 
        error('wrong model type')
    end
else
    %start at left leg, end at right leg
    ind_grf_start = find(strcmp(obj.model.GRFNAMES, 'leftFy'));
    ind_grf_end = find(strcmp(obj.model.GRFNAMES, 'rightFy'));
    if isa(obj.model, 'Gait3d')
        ind_state_start = find(strcmp(obj.model.states.name, 'CPM_l') & strcmp(obj.model.states.type,'xc'));
        ind_state_end = find(strcmp(obj.model.states.name, 'CPM_r') & strcmp(obj.model.states.type,'xc'));
    elseif isa(obj.model, 'Gait2dc')
        ind_state_start = find(strcmp(obj.model.states.name, 'heel_l') & strcmp(obj.model.states.type,'xc'));
        ind_state_end = find(strcmp(obj.model.states.name, 'heel_r') & strcmp(obj.model.states.type,'xc'));
    else 
        error('wrong model type')
    end
end

if obj.isSymmetric
    %Always exactly half the stride length
    ind_state = find(strcmp(obj.model.states.name, 'pelvis_tx') & strcmp(obj.model.states.type,'q'));
    SL = x(ind_state,end)-x(ind_state,1);
else
    %create two cycles
    GRF(:,obj.nNodes+(1:obj.nNodes)) = GRF;

    %First find when the first heel strike of the opposite leg is
    ind_swing = find(GRF(ind_grf_start,:) < 0.01,1);
    HS_opp = find(GRF(ind_grf_start,ind_swing:end) > 0.05,1);
    HS_opp = HS_opp + ind_swing;
    if HS_opp > obj.nNodes
        HS_opp = HS_opp - obj.nNodes;
    end

    %Find where actual leg heel strikes after this heel strike
    ind_swing = find(GRF(ind_grf_end,HS_opp:end) < 0.01,1);
    ind_swing = ind_swing + HS_opp;
    HS_leg = find(GRF(ind_grf_end,ind_swing:end) > 0.05,1);
    HS_leg = HS_leg+ind_swing;

    if HS_leg > obj.nNodes
        HS_leg = HS_leg - obj.nNodes;
        SL = x(ind_state_end, obj.nNodes)-x(ind_state_start, HS_opp) ;
        SL = SL + (x(ind_state_start, HS_leg)-x(ind_state_start, 1));
    else
        SL = x(ind_state_end, HS_leg)-x(ind_state_start, HS_opp);
    end
end
% if obj.isSymmetric
%     if HS_leg > obj.nNodes
%         HS_leg = 
%     SL = 1;  
% else
%     SL = x(ind_state_end,ind_grf_end)-x(ind_state_start,ind_grf_end);
% end
    


% add_length = 0; % for handling symmetry
% if obj.isSymmetric
%     
% else
%     if isa(obj.model, 'Gait3d')
%         ind_state_start = find(strcmp(obj.model.states.name, 'CPM_l') & strcmp(obj.model.states.type,'xc'));
%         if end_ind > obj.nNodes
%             end_ind = end_ind - obj.nNodes;
%             ind_state_end = find(strcmp(obj.model.states.name, 'CPM_r') & strcmp(obj.model.states.type,'xc'));
%             add_length = x(ind_state_start,obj.nNodes)-x(ind_state_start,1);
%         else
%             ind_state_end = find(strcmp(obj.model.states.name, 'CPM_l') & strcmp(obj.model.states.type,'xc'));
%         end
%     elseif isa(obj.model, 'Gait2dc')
%         ind_state_start = find(strcmp(obj.model.states.name, 'heel_l') & strcmp(obj.model.states.type,'xc'));
%         if end_ind > obj.nNodes
%             end_ind = end_ind - obj.nNodes;
%             ind_state_end = find(strcmp(obj.model.states.name, 'heel_r') & strcmp(obj.model.states.type,'xc'));
%             add_length = x(ind_state_start,obj.nNodes)-x(ind_state_start,1);
%         else
%             ind_state_end = find(strcmp(obj.model.states.name, 'heel_l') & strcmp(obj.model.states.type,'xc'));
%         end
%     else 
%         error('wrong model type')
%     end
% end
% SL = x(ind_state_end,end_ind)-x(ind_state_start,start_ind)+add_length;