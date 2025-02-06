%======================================================================
%> @file Collocation/plotMeanSimVarTable.m
%> @brief Collocation function to plot the mean and SD of multiple simulation variable tables
%> @details
%> Details: Collocation::plotMeanSimVarTable()
%>
%> @author Marlies Nitschke
%> @date November, 2018
%======================================================================

%======================================================================
%> @brief Computes function to plot the mean and SD of multiple simulation variable tables
%> @static
%>
%> @details
%> The function computes and plots the mean and variance for each column of tables.
%>
%> When plotting IMU and marker data separate subplots are created for different
%> directions and positions.
%>
%> The function is not yet able to plot stance phases.
%>
%> @code
%> Collocation.plotMeanSimVarTable({simVarTable1, simVarTable2; simVarTable3, simVarTable4}, style);
%> @endcode
%>
%> @note
%> The function is currently assuming that tracking and reference data is identical over all tables
%> within each column and is therefore only plotting the tracking and reference data of the first row.
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
%>                           - subFigSettings
%>                           .
%>                         If a field is not defined, default values specified in the description of plotVarType() will be used.
%>                         Use empty to skip the input.
%> @retval meanTablesOut   Cell vector: SimVarTables reflecting the mean of each column in tables. The simulated means and variances
%>                         are saved in the table columns 'sim' and 'simVar', respectively. The columns 'mean', 'var',
%>                         'mean_extra', and 'var_extra' reflect the tracking and reference data which is currently assumed to be
%>                         equal for all tables within one column.
%======================================================================
function meanTablesOut = plotMeanSimVarTable(tables, style)

if nargin < 2 || isempty(style)
   style = struct();
end

%% Remove rows with types which are not plotted from meanTables
% Remove rows
excludedTypes = {'standing', 'standingEvent', 'footAngleEvent', 'dur', 'speed'};
for iTab = 1 : numel(tables)
   tables{iTab} = tables{iTab}(~ismember(tables{iTab}.type, excludedTypes), :);
end

%% Get mean table (assuming same order! and same variables!)
nMeanTables = size(tables, 2);
meanTables = cell(1, nMeanTables);
warning('off', 'MATLAB:table:RowsAddedNewVars');
% Get means for each column
for iTab = 1 : nMeanTables
    % Initialize using the first row
    if any(ismember(tables{1, iTab}.type, {'acc', 'gyro', 'marker'})) % IMU and marker data
        curMeanTable = tables{1, iTab}(:, {'type', 'name', 'unit', 'direction', 'position'});
    else
        curMeanTable = tables{1, iTab}(:, {'type', 'name', 'unit'});
    end
    % Add columns for mean and variance of tracking data (We assume here
    % that tracking data is equal for all tables)
    curMeanTable(:, 'mean') = tables{1, iTab}(:, 'mean');
    curMeanTable(:, 'var')  = tables{1, iTab}(:, 'var');
    % extra data
    curMeanTable(:, 'mean_extra') = tables{1, iTab}(:, 'mean_extra');
    curMeanTable(:, 'var_extra')  = tables{1, iTab}(:, 'var_extra');
    
    for iTabRow = 1 : height(curMeanTable)
        
        for iRow = 1 : size(tables, 1)
            curTab = tables{iRow, iTab};
            simMatrix(:, iRow) = curTab.sim{iTabRow};
        end
        curMeanTable.sim{iTabRow} = mean(simMatrix, 2); % along columns
        curMeanTable.simVar{iTabRow} = var(simMatrix, 0, 2); % along columns
               
    end
    
    meanTables{iTab} = curMeanTable;
    clear simMatrix % Has to be cleared to use data with different number of nodes
end
warning('on', 'MATLAB:table:RowsAddedNewVars');
meanTablesOut = meanTables;

%% Find unique mean and var to plot mean and var only once if tracking data was equal
% match table rows by identifier and do not assume the same order
% (Attention: Not working for IMU data if mean of multiple axis has the
% same mean which will basically not be the case.)
for iTab = 1 : nMeanTables
    curTab = meanTables{iTab}(:, {'type', 'name', 'mean', 'var'});
    if isempty([curTab.mean{:}]); continue; end % All empty => nothing to compare
    
    for iTabCompare = iTab+1 : nMeanTables % do not compare twice
        curTabComp = meanTables{iTabCompare}(:, {'type', 'name', 'mean', 'var'});
        
        % compare each row with each row
        for iRowCur = 1:height(curTab)
            for iRowComp = 1:height(curTabComp)
                
                if strcmp(curTab.type{iRowCur}, curTabComp.type{iRowComp}) ...
                        && strcmp(curTab.name{iRowCur}, curTabComp.name{iRowComp}) ...
                        && isequaln(curTab.mean{iRowCur}, curTabComp.mean{iRowComp}) ...
                        && isequaln(curTab.var{iRowCur}, curTabComp.var{iRowComp})
                    
                    meanTables{iTabCompare}.mean(iRowComp, :) = cell(1, 1);  % clear mean column of respective table row
                    meanTables{iTabCompare}.var(iRowComp, :) = cell(1, 1);  % clear var column of respective table row
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
for iTab = 1 : nMeanTables
    curTab = meanTables{iTab}(:, {'type', 'name', 'mean_extra', 'var_extra'});
    if isempty([curTab.mean_extra{:}]); continue; end % All empty => nothing to compare
    
    for iTabCompare = iTab+1 : nMeanTables % do not compare twice
        curTabComp = meanTables{iTabCompare}(:, {'type', 'name', 'mean_extra', 'var_extra'});

        % compare each row with each row
        for iRowCur = 1:size(curTab,1)
            for iRowComp = 1:size(curTabComp,1)

                if strcmp(curTab.type{iRowCur}, curTabComp.type{iRowComp}) ...
                        && strcmp(curTab.name{iRowCur}, curTabComp.name{iRowComp}) ...
                        && isequaln(curTab.mean_extra{iRowCur}, curTabComp.mean_extra{iRowComp}) ...
                        && isequaln(curTab.var_extra{iRowCur}, curTabComp.var_extra{iRowComp})

                    meanTables{iTabCompare}.mean_extra(iRowComp, :) = cell(1, 1);  % clear mean column of respective table row
                    meanTables{iTabCompare}.var_extra(iRowComp, :) = cell(1, 1);  % clear var column of respective table row
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
for iTab = 1 : nMeanTables
    if any(ismember(meanTables{iTab}.type, {'acc', 'gyro', 'marker'}))
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
for iTab = 1 : nMeanTables
    allEntries = [allEntries; meanTables{iTab}(:, identifiers)];
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
colors = lines(nMeanTables);

%% Plot all
for iTab = 1 : nMeanTables
   curTab = meanTables{iTab};
   curTypes = unique(curTab.type);
   style.simColor = colors(iTab, :);
    
   % Plot
   for iType = 1 : numel(curTypes)
       meanSimVarTable = curTab(ismember(curTab.type, curTypes{iType}), :);
       style.nSubFig = nPlotTable.nPlots{ismember(nPlotTable.type, curTypes{iType})};
       
       % Make the figure
       [hFig, style] = plotVarType(meanSimVarTable, style);
       
   end
end

%% Plot legend for tables
type = 'legend';
hFig = figure(sum(double(type)));
set(hFig, 'Name', type);
set(gcf,'units','centimeters','outerposition', style.figureSize);
set(gcf,'defaultTextInterpreter','latex');
hold on;
legendEntries = cell(nMeanTables, 1);
for iTab = 1 : nMeanTables
    x = 1:100;
    color = colors(iTab, :);
    plot(x, zeros(size(x))-iTab, 'Color', color, 'LineWidth', 3);    
    legendEntries{iTab} = sprintf('Col: %i', iTab);
end
box off;
axis off;
legend(legendEntries, 'Location', 'eastoutside');

end