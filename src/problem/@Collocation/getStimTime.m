%======================================================================
%> @file @Collocation/getStimTime.m
%> @brief Collocation function to calculate the stimulation time for one
%> gait cycle
%> @details
%> Details: Collocation::getStimTime()
%>
%> @author Anne Koelewijn, Markus Gambietz
%> @date December, 2019
%======================================================================

%======================================================================
%> @brief Collocation function to calculate the stimulation time for one
%> gait cycle
%> @details
%> The stimulation time is a variable that is necessary for some of the
%> metabolic models, because the energy expenditure is dependent on the
%> stimulation time. This function calculates the stimulation time at the
%> beginning of the gait cycle by looping to the gait cycle once.
%> 
%> To get the cost of a result, you can call:
%> @code
%> t_stim = problem.getStimTime(X);
%> @endcode
%>
%>
%> @param  obj            Collocation object
%> @param  X              Double matrix: State vector (i.e. result) of the problem
%> @retval t_stim         Double matrix: Activation Time (nMus x nNodesDur)
%======================================================================
function t_stim = getStimTime(obj, X)
% Get muscle activations:
act = X(obj.idx.states(obj.model.extractState('a'),1:obj.nNodes));
% Threshold for activation:
t = 0.1;
% Time step
h = X(obj.idx.dur)/obj.nNodes;

% In the periodic case, we can expand the 'act' array by a full gait cycle
if obj.isPeriodic
    if obj.isSymmetric % In the symmetric case, first expand it to a full gait cycle
        % Get the symmetry indices
        sym_idx = obj.model.idxSymmetry.xindex(obj.model.extractState('a'));
        % Expand act by the symmetric condition
        act = [X(obj.idx.states(sym_idx,1:obj.nNodes)) act];
    end
    % Expand by a full gait cycle
    act = [act act];
end
% Set mask where
act_bool = act > t;
act_bool_interp = double(act_bool); % For interpolating zeros

% Interpolate when an activation just started, see old 'calcStimTime'
iThresh = find(max(0,diff(act_bool,[],2))); % Find the indices before the thresholds
valThresh = (0.1-act(iThresh))./ ...
    (act(iThresh+height(act))-act(iThresh)); % Linear Interpolate the values at the beginning of the activations
act_bool_interp(iThresh+height(act)) = valThresh; % Set the interpolated values

% Calculate t_stim using logical operators
t_stim = cumsum(act_bool_interp,2);
t_mask = t_stim .* ~act_bool; % Mask 0 activations to remove them later
t_stim = t_stim - cummax(t_mask,2);

% Finalize output: t_stim is h * the last nNodes
t_stim = t_stim(:,end-obj.nNodes+1:end) * h;

% Warn when t_stim is constantly > threshold in any muscle in periodic
% movements
if any(min(act_bool,[],2)) && obj.isPeriodic
    warning(['Collocation.getStimTime: A muscle counts as constantly activated.' ...
        'The resulting metabolic cost is wrong.'])
end
end

