%======================================================================
%> @file Collocation/plotMultSimVarTables.m
%> @brief Collocation function to plot the content of multiple simulation variable tables
%> @details
%> Details: Collocation::plotMultSimVarTables()
%>
%> @author Marlies Nitschke
%> @date November, 2018
%======================================================================

%======================================================================
%> @brief Function to plot the content of multiple simulation variable tables
%> @static
%>
%> @details
%> Usage:
%> @code
%> Collocation.plotMultSimVarTables({simVarTable1, simVarTable2; simVarTable3, simVarTable4}, style);
%> @endcode
%>
%> It is also possible to plot simVarTables with different content.
%>
%> When plotting IMU and marker data separate subplots are created for different
%> directions and positions.
%>
%> @param  tables          Cell matrix: Each entry is a simVarTable containing a full gait cycle returned
%>                         from Collocation.report or Collocation.extractData.
%>                         The results in dimension 1 can for example be from another subject and in dimension 2 from another condition, e.g., walking with an exoskeleton.
%> @param  style           (optional) Struct: Style of the figure containing the fields:
%>                           - figureSize
%>                           - xLabelFontSize
%>                           - yLabelFontSize
%>                           - xTickFontSize
%>                           - yTickFontSize
%>                           - legendFontSize
%>                           - lineWidth
%>                           - trackColor
%>                           - extraColor
%>                           - trackFaceAlpha
%>                           - extraFaceAlpha
%>                           - standingRightColor
%>                           - standingLeftColor
%>                           - standingFaceAlpha
%>                           - subFigSettings
%>                           .
%>                         If a field is not defined, default values specified in the description of plotVarType() will be used.
%>                         Use empty to skip the input.
%> @param  plotStance      (optional) Boolean: If true and if the tables contain rows of type "standing",
%>                         the stance phase(s) will be plotted. (default: false)
%======================================================================
function plotMultSimVarTables(tables, style, plotStance)

if nargin < 2 || isempty(style)
   style = struct();
end

if nargin < 3
   plotStance = 0;
end

nTables = numel(tables);

%% Remove rows with types which are not plotted from tables
% Get standing entries if it should be ploted
if plotStance
    standings = cell(size(tables));
    for iTab = 1 : nTables
        standings{iTab} = tables{iTab}(ismember(tables{iTab}.type, 'standing'), {'type', 'name', 'sim'});
    end
end
% Remove rows
excludedTypes = {'standing', 'standingEvent', 'footAngleEvent', 'dur', 'speed'};
for iTab = 1 : nTables
   tables{iTab} = tables{iTab}(~ismember(tables{iTab}.type, excludedTypes), :);
end


%% Find unique mean and var to plot mean and var only once if tracking data was equal
% match table rows by identifier and do not assume the same order
% (Attention: Not working for IMU data if mean of multiple axis has the
% same mean which will basically not be the case.)
for iTab = 1 : nTables
    curTab = tables{iTab}(:, {'type', 'name', 'mean', 'var'});
    if isempty([curTab.mean{:}]); continue; end % All empty => nothing to compare
    
    for iTabCompare = iTab+1 : nTables % do not compare twice
        curTabComp = tables{iTabCompare}(:, {'type', 'name', 'mean', 'var'});
        
        % compare each row with each row
        for iRowCur = 1:height(curTab)
            for iRowComp = 1:height(curTabComp)
                
                if strcmp(curTab.type{iRowCur}, curTabComp.type{iRowComp}) ...
                        && strcmp(curTab.name{iRowCur}, curTabComp.name{iRowComp}) ...
                        && isequaln(curTab.mean{iRowCur}, curTabComp.mean{iRowComp}) ...
                        && isequaln(curTab.var{iRowCur}, curTabComp.var{iRowComp})
                    
                    tables{iTabCompare}.mean(iRowComp, :) = cell(1, 1);  % clear mean column of respective table row
                    tables{iTabCompare}.var(iRowComp, :) = cell(1, 1);  % clear var column of respective table row
                    break; % assume that a single row is unique within a table
                end

            end
        end        
    end
end

%% Find unique extra mean and var to plot mean and var only once if extra data was equal
% match table rows by identifier and do not assume the same order
% (Attention: Not working for IMU data if mean of multiple axis has the
% same mean which will basically not be the case.)
for iTab = 1 : nTables
    curTab = tables{iTab}(:, {'type', 'name', 'mean_extra', 'var_extra'});
    if isempty([curTab.mean_extra{:}]); continue; end % All empty => nothing to compare
    
    for iTabCompare = iTab+1 : nTables % do not compare twice
        curTabComp = tables{iTabCompare}(:, {'type', 'name', 'mean_extra', 'var_extra'});
        
        % compare each row with each row
        for iRowCur = 1:size(curTab,1)
            for iRowComp = 1:size(curTabComp,1)
                
                if strcmp(curTab.type{iRowCur}, curTabComp.type{iRowComp}) ...
                        && strcmp(curTab.name{iRowCur}, curTabComp.name{iRowComp}) ...
                        && isequaln(curTab.mean_extra{iRowCur}, curTabComp.mean_extra{iRowComp}) ...
                        && isequaln(curTab.var_extra{iRowCur}, curTabComp.var_extra{iRowComp})

                    tables{iTabCompare}.mean_extra(iRowComp, :) = cell(1, 1);  % clear mean column of respective table row
                    tables{iTabCompare}.var_extra(iRowComp, :) = cell(1, 1);  % clear var column of respective table row
                    break; % assume that a single row is unique within a table
                end

            end
        end        
    end
end

%% find all possible entries for type and name 
% -> determine required amount of subplots for each figure

% Determine the columns which are used as identifier based on the fact that
% IMU and/or marker data is contained or not
containsIMUMarker = false;
for iTab = 1 : nTables
    if any(ismember(tables{iTab}.type, {'acc', 'gyro', 'marker'}))
        containsIMUMarker = true;
        break
    end
end

if ~containsIMUMarker
    identifiers = {'type', 'name'};
else
    identifiers = {'type', 'name', 'direction', 'position'};
end

% Get all variables which are in all tables
allEntries = table();
for iTab = 1 : nTables
    allEntries = [allEntries; tables{iTab}(:, identifiers)];
end

% Find unique entries of all tables
isIMUMarker = ismember(allEntries.type, {'acc', 'gyro', 'marker'});
allEntriesIMUMarker = unique(allEntries(isIMUMarker, :));
allEntriesOthers = unique(allEntries(~isIMUMarker, {'type', 'name'}));

allTypesIMUMarker = unique(allEntriesIMUMarker.type);
allTypesOthers = unique(allEntriesOthers.type);
allTypes = [allTypesIMUMarker; allTypesOthers];

nTypesIMUMarker = numel(allTypesIMUMarker);
nTypesOthers = numel(allTypesOthers);
nTypes = numel(allTypes);

% Define a new table containing the unique types of all tables and the
% required number of subplots -> pass that table to the plot function!
nPlotTable = table(allTypes, cell(nTypes,1), 'VariableNames', {'type', 'nPlots'});
for iType = 1 : nTypesIMUMarker
    nPlotTable.nPlots{iType} = sum(ismember(allEntriesIMUMarker.type, allTypesIMUMarker{iType}));
end
for iType = 1 : nTypesOthers
    nPlotTable.nPlots{nTypesIMUMarker+iType} = sum(ismember(allEntriesOthers.type, allTypesOthers{iType}));
end
    
%% Get colors
nRows = size(tables, 1);
nCols = size(tables, 2);
colors = lines(nCols);  
intensity = linspace(0.6, 1, nRows);


%% Plot all
for iTab = 1 : nTables
    curTab = tables{iTab};
    curTypes = unique(curTab.type);
    [iRow, iCol] = ind2sub(size(tables), iTab);
    style.simColor = intensity(iRow) * colors(iCol, :);

    % Get standing
    standing_r = [];
    standing_l = [];
    if plotStance
        iStandingR = find(strcmp(standings{iTab}.name, 'standing_r'));
        if ~isempty(iStandingR)
            standing_r = standings{iTab}.sim{iStandingR};
        end
        iStandingL = find(strcmp(standings{iTab}.name, 'standing_l'));
        if ~isempty(iStandingL)
            standing_l = standings{iTab}.sim{iStandingL};
        end
    end

    % Plot
    for iType = 1 : numel(curTypes)
    
        simVarTable = curTab(ismember(curTab.type, curTypes{iType}), :);
        style.nSubFig = nPlotTable.nPlots{ismember(nPlotTable.type, curTypes{iType})};

        % Plot
        [hFig, style] = plotVarType(simVarTable, style, standing_r, standing_l);
        
    end
 
end

%% Plot legend for tables
type = 'legend';
hFig = figure(sum(double(type)));
set(hFig, 'Name', type);
if ~strcmp(get(gcf,'WindowStyle'), 'docked')   % prevents warning from outerposition in next line
    set(gcf,'units','centimeters','outerposition',style.figureSize);
end
set(gcf,'defaultTextInterpreter','latex');
hold on;
legendEntries = cell(nTables, 1);
for iTab = 1 : nTables
    x = 1:100;
    [iRow, iCol] = ind2sub(size(tables), iTab);
    color = intensity(iRow) * colors(iCol, :);
    plot(x, zeros(size(x))-iTab, 'Color', color, 'LineWidth', 3);    
    legendEntries{iTab} = sprintf('Row: %i, Col: %i', iRow, iCol);
end
box off;
axis off;
legend(legendEntries, 'Location', 'eastoutside');

end