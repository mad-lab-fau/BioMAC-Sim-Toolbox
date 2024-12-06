%======================================================================
%> @file resampleInitialContact.m
%> @brief Function to shift the result.X vector to start with an initial contact
%> @details
%> Details: Collocation::resampleInitialContact()
%>
%> @author Markus Gambietz
%> @date April, 2023
%======================================================================

%======================================================================
%> @brief
%> Resamples X so that the Heelstrike is at the first node. Naive
%>   implementation that only considers one Heelstrike per foot
%> @todo Implement symmetric case, 3d-Curved running (?)
%>
%> @param  obj           Collocation class object
%> @param  X             Double array: State vector containing at least 'speed'
%> @param  foot          left ('l') or right (default, 'r')
%> @param  contact_names [2 x n] cell containing all the contact names in the
%model, right_side and left side separated (see default value)
%======================================================================
function X = resampleInitialContact(obj, X, foot, contact_names)
% Check number of input arguments, set default values
if nargin <= 2
    foot = 'r';
end
if nargin <= 3
    contact_names = {'heel_r','front_r';'heel_l','front_l'};
end
% Check symmetry: Assymetric Gait not yet implemented
if ~obj.isSymmetric
    error('Collocation.resampleHeelstrike: Only symmetric gait is supported')
end

contact_size = size(contact_names);
idx_r = [];
idx_l = [];

% Extract GRF
for i = 1:contact_size(2)
    idx_r = [idx_r obj.model.extractState('Fy', contact_names{1,i})];
    idx_l = [idx_l obj.model.extractState('Fy', contact_names{2,i})];
end

% Concat GRF to a single vector
if strcmp(foot,'r')
    Fy = [sum(X(obj.idx.states(idx_r,1:obj.nNodes))) sum(X(obj.idx.states(idx_l,1:obj.nNodes)))];
elseif strcmp(foot, 'l')
    Fy = [sum(X(obj.idx.states(idx_l,1:obj.nNodes))) sum(X(obj.idx.states(idx_r,1:obj.nNodes)))];
end

% Threshold Fy to 0.03, find the 
F_bool = (Fy > 0.03);
idx_start = strfind(F_bool,[0 1]);
X = shift_index(X, idx_start);

    function X_tmp = shift_index(X_tmp, amount)
        idx_states = obj.idx.states;
        idx_controls = obj.idx.controls;
        idx_tx = obj.model.extractState('q','pelvis_tx');
        for k = 1:amount   
            %% States
            % Remember how much the pelvis moves on the first node
            d_pelvis_dx = X_tmp(idx_states(idx_tx,2)) - X_tmp(idx_states(idx_tx,1));
            % Move everything one to the left
            X_tmp(idx_states(:,1:end-1)) = X_tmp(idx_states(:,2:end));
            % Take the new first and make it the last on the other side
            X_tmp(idx_states(:,end)) = X_tmp(idx_states(obj.model.idxSymmetry.xindex,1));
            % Adjust the last pelvis_x: Last node is the old first node, copy
            % the distance
            X_tmp(idx_states(idx_tx,end)) = X_tmp(idx_states(idx_tx,end-1)) + d_pelvis_dx;
            % Shift back, so that the X_coordinate
            X_tmp(idx_states(idx_tx,:)) = X_tmp(idx_states(idx_tx,:)) - d_pelvis_dx;
            %% Controls
            % Move everything one to the left
            X_tmp(idx_controls(:,1:end-1)) = X_tmp(idx_controls(:,2:end));
            % Take the new first and make it the last on the other side
            X_tmp(idx_controls(:,end)) = X_tmp(idx_controls(obj.model.idxSymmetry.uindex',1));            
        end
    end
end

