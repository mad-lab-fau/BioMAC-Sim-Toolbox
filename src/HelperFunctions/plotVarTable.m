%======================================================================
%> @file plotVarTable.m
%> @brief Function to plot all variables within a variables table
%> @details
%> Details: plotVarTable()
%>
%> @author Marlies Nitschke
%> @date January, 2019
%======================================================================

% ======================================================================
%> @brief Function to plot all variables within a variables table
%>
%> @param   varTable       Table: Table containing the data of different variables
%>                         with the following columns:
%>                         - type
%>                         - name
%>                         - unit
%>                         - direction (only required for acc, gyro and marker data)
%>                         - sim (optional)
%>                         - simVar (optional)
%>                         - mean (optional)
%>                         - var (optional)
%>                         - mean_extra (optional)
%>                         - var_extra (optional)
%>                         .
%>                         The data of the available columns will be plotted.
%> @param  style           (optional) Struct: Settings defining the style for plotting.
%>                         See plotVarType() for details.
%> @param  plotStance      (optional) Boolean: If true and if the tables contain rows of type "standing"
%>                         where the stance is defined in the column "sim", the stance phase(s)
%>                         will be plotted. (default: false)
% ======================================================================
function plotVarTable(varTable, style, plotStance)

if nargin < 2 || isempty(style)
   style = struct();
end

if nargin < 3
   plotStance = 0;
end

% Get standing entries if it should be ploted
standing_r = [];
standing_l = [];
if plotStance
    iStandingR = find(strcmp(varTable.type, 'standing') & strcmp(varTable.name, 'standing_r'));
    if ~isempty(iStandingR)
        standing_r = varTable.sim{iStandingR};
    end
    iStandingL = find(strcmp(varTable.type, 'standing') & strcmp(varTable.name, 'standing_l'));
    if ~isempty(iStandingL)
        standing_l = varTable.sim{iStandingL};
    end
end

% Remove rows with types which are not plotted from tables
excludedTypes = {'standing', 'standingEvent', 'footAngleEvent', 'dur', 'duration', 'speed', 'traveledSpeed'};
varTable = varTable(~ismember(varTable.type, excludedTypes), :);

% Plot all
types = unique(varTable.type);
for iType = 1 : numel(types)

    varTableType = varTable(ismember(varTable.type, types{iType}), :);
    [hFig, style] = plotVarType(varTableType, style, standing_r, standing_l);

end

end
