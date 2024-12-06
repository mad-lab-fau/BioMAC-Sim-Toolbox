%======================================================================
%> @file @Collocation/checkXScaling.m
%> @brief Collocation function to check the scaling withing X
%> @details
%> Details:  Collocation::checkXScaling()
%>
%> @author Marlies Nitschke
%> @date March, 2019
%======================================================================

%======================================================================
%> @brief Matlab function to check the scaling withing X
%>
%> @details
%> It plots the ranges (min and max) of the absolute values in X for each
%> type of states and controls. This function can be used to check whether
%> the scaling of the states and controls
%>
%> The function can be called using: 
%> @code
%> result.problem.checkXScaling(result.X);
%> @endcode
%>
%> @param  obj          Collocation class object
%> @param  X            Double matrix: State vector (i.e. result) of the problem
%======================================================================
function checkXScaling(obj, X)

% Get minimum and maximum absolute states and controls
states = X(obj.idx.states);
controls = X(obj.idx.controls);
maxStates = max(abs(states), [], 2);
maxControls = max(abs(controls), [], 2);
minStates = min(abs(states), [], 2);
minControls = min(abs(controls), [], 2);

minStatesControls = [minStates; minControls];
maxStatesControls = [maxStates; maxControls];

% Get types for each entry
statesControlsTypes = [obj.model.states.type; obj.model.controls.type];
types = unique(statesControlsTypes, 'stable');
nTypes = length(types);

% Get minimum and maximum among all variables of each type
minStatesControlsEachType = nan(nTypes, 1);
maxStatesControlsEachType = nan(nTypes, 1);
for iType = 1 : nTypes
    idxCurType = find(strcmp(statesControlsTypes, types{iType}));
    minStatesControlsEachType(iType) = min(minStatesControls(idxCurType));
    maxStatesControlsEachType(iType) = max(maxStatesControls(idxCurType));
end

% Plot a graph to show the range of absolute values for each tyoe
figure();
h = bar([minStatesControlsEachType, maxStatesControlsEachType-minStatesControlsEachType],'stacked','BaseValue',min(minStatesControlsEachType));
h(1).Visible = 'off';
set(gca,'YLim',[min(minStatesControlsEachType),max(maxStatesControlsEachType)+1]);
xticks(1:nTypes);
xticklabels(types); xtickangle(90);
title('Range of absolute values of states and controls');

end

