%======================================================================
%> @file HelperFunctions/plotVarType.m
%> @brief Function to plot one type of variable of a variables table
%> @details
%> Details: plotVarType()
%>
%> @author Marlies Nitschke
%> @date July, 2021
%======================================================================

%======================================================================
%> @brief Function to plot one type of variable of a variables table
%>
%> @details
%> Default values of the fields of style:
%> @code
%> style.figureSize         = [0 0 16 20]; % in cm
%> style.xLabelText         = 'Gait Cycle in \%';
%> style.xLabelFontSize     = 10;
%> style.yLabelFontSize     = 10;
%> style.xTickFontSize      = 11;
%> style.yTickFontSize      = 11;
%> style.legendFontSize     = 9;
%> style.lineWidth          = 0.5;
%> style.simColor           = 'k';   % or a rgb vector
%> style.trackColor         = 'r--'; % or a rgb vector
%> style.extraColor         = 'k--'; % or a rgb vector
%> style.trackFaceAlpha     = 0.1;
%> style.extraFaceAlpha     = 0.1;
%> style.standingRightColor = 'k';
%> style.standingLeftColor  = 'k';
%> style.standingFaceAlpha  = 0.1;
%> style.subFigSettings.nCol        = 2;
%> style.subFigSettings.width       = 5;   % in cm
%> style.subFigSettings.height      = 3.5; % in cm
%> style.subFigSettings.originUp    = 4;   % in cm
%> style.subFigSettings.originRight = 2;   % in cm
%> style.subFigSettings.oneUp       = 1.5; % in cm
%> style.subFigSettings.oneRight    = 2.5; % in cm
%> @endcode
%> The default values showed to work well when using the figure for DIN A4. It might however not be
%> well suited for your screen or your purposes. Especially the option nCol might be useful.
%>
%> @param   varTable       Table: Table containing the data of ONE type of variable
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
%> @param   style          (optional) Struct: Style of the figure defined by variouse optional fields (see description above).
%>                         If a field is not defined, default values specified in the description above will be used.
%>                         If the content of multiple function calls should be plotted in the same figures,
%>                         additionally the field "nSubFig" has to be specified with the number of subfigures.
%>                         Use empty to skip the input.
%> @param   standing_r     (optional) Double vector: Stance phase of right foot with 0: no standing; 1: standing (size of other data)
%>                         Use empty vector to skip it.
%> @param   standing_l     (optional) Double vector: Stance phase of left foot with 0: no standing; 1: standing (size of other data)
%> @retval  hFig           Handle of figure
%> @retvl   style          Struct: Style of figure containing all fields
%======================================================================
function [hFig, style] = plotVarType(varTable, style, standing_r, standing_l)

    % Check if table contains only one type
    type = unique(varTable.type);
    if numel(type) > 1
        msg = ['The varTable can only have one type. ', ...
               'Currently there are the types: ', ...
               sprintf('%s  ', type{:})];
        error(msg);
    end
    
    % Get default values for style if not given
    if nargin < 2 || isempty(style); style = struct(); end
    if ~isfield(style, 'figureSize'); style.figureSize = [0 0 16 20]; end
    if ~isfield(style, 'xLabelText'); style.xLabelText = 'Gait Cycle in \%'; end
    if ~isfield(style, 'xLabelFontSize'); style.xLabelFontSize = 10; end
    if ~isfield(style, 'yLabelFontSize'); style.yLabelFontSize = 10; end
    if ~isfield(style, 'xTickFontSize'); style.xTickFontSize = 11; end
    if ~isfield(style, 'yTickFontSize'); style.yTickFontSize = 11; end
    if ~isfield(style, 'legendFontSize'); style.legendFontSize = 9; end
    if ~isfield(style, 'lineWidth'); style.lineWidth = 0.5; end
    if ~isfield(style, 'simColor'); style.simColor = 'k'; end
    if ~isfield(style, 'trackColor'); style.trackColor = 'r--'; end
    if ~isfield(style, 'extraColor'); style.extraColor = 'k--'; end
    if ~isfield(style, 'trackFaceAlpha'); style.trackFaceAlpha = 0.1; end
    if ~isfield(style, 'extraFaceAlpha'); style.extraFaceAlpha = 0.1; end
    if ~isfield(style, 'standingRightColor'); style.standingRightColor = 'k'; end
    if ~isfield(style, 'standingLeftColor'); style.standingLeftColor = 'k'; end
    if ~isfield(style, 'standingFaceAlpha'); style.standingFaceAlpha = 0.1; end

    if ~isfield(style, 'subFigSettings'); style.subFigSettings = struct(); end
    if ~isfield(style.subFigSettings, 'width'); style.subFigSettings.width = 5; end
    if ~isfield(style.subFigSettings, 'height'); style.subFigSettings.height = 3.5; end
    if ~isfield(style.subFigSettings, 'originUp'); style.subFigSettings.originUp = 4; end
    if ~isfield(style.subFigSettings, 'originRight'); style.subFigSettings.originRight = 2; end
    if ~isfield(style.subFigSettings, 'oneUp'); style.subFigSettings.oneUp = 1.5; end
    if ~isfield(style.subFigSettings, 'oneRight'); style.subFigSettings.oneRight = 2.5; end
    if ~isfield(style.subFigSettings, 'nCol'); style.subFigSettings.nCol = 5; end

    % Get stuff from style struct
    figureSize         = style.figureSize;
    xLabelText         = style.xLabelText;
    xLabelFontSize     = style.xLabelFontSize;
    yLabelFontSize     = style.yLabelFontSize;
    xTickFontSize      = style.xTickFontSize;
    yTickFontSize      = style.yTickFontSize;
    legendFontSize     = style.legendFontSize;  
    lineWidth          = style.lineWidth;
    simColor           = style.simColor;
    trackColor         = style.trackColor;
    extraColor         = style.extraColor;
    trackFaceAlpha     = style.trackFaceAlpha;
    extraFaceAlpha     = style.extraFaceAlpha;
    standingRightColor = style.standingRightColor;
    standingLeftColor  = style.standingLeftColor;
    standingFaceAlpha  = style.standingFaceAlpha;
    subFigSettings     = style.subFigSettings;
    nCol               = style.subFigSettings.nCol;

    % Set stuff which is need
    % Name of type
    type = type{:};
    % Number of subfigures
    nSubFigNow = height(varTable); % Plot all entries in separate figure
    if isfield(style, 'nSubFig')
        nSubFig = style.nSubFig; % Additional input if function is called multiple times with different number of rows
        enableMultipleCalls = 1;
    else
        nSubFig = nSubFigNow;
        enableMultipleCalls = 0;
    end
    allEntries = {'Simulated signal', 'Simulated mean $\pm$ SD', 'Tracked signal', 'Tracked mean $\pm$ SD', 'Extra signal', 'Extra mean $\pm$ SD', 'Stance Phase of Right Foot', 'Stance Phase of Left Foot'}; % legend entries
    hSubFigs = cell(nSubFig, 1); % We need all handles to assign the legend the subplot containing the most data (also tracking). Otherwise this is not working matlab2tikz()
    allHandels = cell(nSubFig, numel(allEntries)); % Handles which should be in legend for each plot
    allFaceAlphas = [nan, trackFaceAlpha, nan, trackFaceAlpha, nan, extraFaceAlpha, standingFaceAlpha, standingFaceAlpha];

    supportedColumns = {'sim', 'simVar', 'mean', 'var', 'mean_extra', 'var_extra'};
    dataColumns = intersect(varTable.Properties.VariableNames, supportedColumns, 'stable');
    if isempty(dataColumns)
        msgSupColumns   = sprintf('''%s'', ', supportedColumns{:});
        error('None of the supported data columns is defined in varTable. The following columns are supported: %s', msgSupColumns(1:end-2));
    end

    % Time vector
    nNodes = 0;
    iCol = 1;
    while nNodes == 0 && iCol <= numel(dataColumns) % Check all columns which might contain data (but assume equal size)
        idxNonEmpty = find(~cellfun(@ isempty, varTable.(dataColumns{iCol}))); % Check all rows
        if ~isempty(idxNonEmpty)
            nNodes = numel(varTable.(dataColumns{iCol}){idxNonEmpty(1)}); % take first non empty
        end
        iCol = iCol +1;
    end
    if nNodes == 0 % no data found so far
        if nargin > 2 && ~isempty(standing_r)
            nNodes = numel(standing_r);
        elseif nargin > 3 && ~isempty(standing_l)
            nNodes = numel(standing_l);
        end
    end
    if nNodes > 1 
        time = 100*(1:nNodes)/nNodes;
    else
        time = 1:100;
    end
    
    % Initialize figure
    if enableMultipleCalls
        hFig = figure(sum(double(type)));
        set(hFig, 'Name', type);
        % If figure already exists, get axes tags to find later respective
        % subplot to plot into!
        nOldSubPlots = size(get(gcf, 'children'),1) - 1; % get existing number of subplots (-1 for legend)
        if nOldSubPlots > 0
            subAxes = findobj( get(gcf,'Children'), '-depth', 1, 'type', 'axes'); % all axes without legend
            subAxes = flip(subAxes); % flip to get the matching order of subplots
            subTags = {subAxes.Tag};
        end
    else
        hFig = figure('Name', type);
        clf;
    end
    if ~strcmp(get(gcf,'WindowStyle'), 'docked')   % prevents warning from outerposition in next line
        set(gcf,'units','centimeters','outerposition',figureSize);
    end
    set(gcf,'defaultTextInterpreter','latex');
    
    % Initialize signals to plot with defaul value for the case that
    % columns are not given
    % (variance is by default 0 to not skrew up compuation of the range
    % when mean is given)
    sim      = nan(nNodes, 1);
    simVar   = zeros(nNodes, 1);
    trackAV  = nan(nNodes, 1);
    trackVar = zeros(nNodes, 1);
    extraAV  = nan(nNodes, 1);
    extraVar = zeros(nNodes, 1);

    % Do subplots
    for iSubFig = 1 : nSubFigNow

        % Define signals to plot
        if any(ismember(dataColumns, 'sim'));        sim      = varTable.sim{iSubFig}; end
        if any(ismember(dataColumns, 'simVar'));     simVar   = varTable.simVar{iSubFig}; end
        if any(ismember(dataColumns, 'mean'));       trackAV  = varTable.mean{iSubFig}; end
        if any(ismember(dataColumns, 'var'));        trackVar = varTable.var{iSubFig}; end
        if any(ismember(dataColumns, 'mean_extra')); extraAV  = varTable.mean_extra{iSubFig}; end
        if any(ismember(dataColumns, 'var_extra'));  extraVar = varTable.var_extra{iSubFig}; end
        name = varTable.name{iSubFig};
        unit = varTable.unit{iSubFig};

        if nNodes == 1
            % Make it also fitting to 100 %
            sim = repmat(sim, 100, 1);
            simVar = repmat(simVar, 100, 1);
            assert((isempty(trackAV) || all(isnan(trackAV))) && (isempty(trackVar) || all(isnan(trackVar))), 'Tracking data can not be plotted for nNodes = 1.');
            assert((isempty(extraAV) || all(isnan(extraAV)) || numel(extraAV) == 100) && (isempty(extraVar) || all(isnan(extraVar)) || numel(extraVar) == 100), 'Extra data has to be resampled to 100 nodes before calling report().');
        end

        % For IMU and marker data, add direction to the name
        if ismember(type, {'acc', 'gyro', 'marker'})
            % Get direction
            direction = varTable.direction(iSubFig, :);

            % Add direction to the name
            if isequal(direction, [1, 0, 0]) || isequal(direction, [1, 0])
                name = [name ' x'];
            elseif isequal(direction, [0, 1, 0]) || isequal(direction, [0, 1])
                name = [name ' y'];
            elseif isequal(direction, [0, 0, 1])
                name = [name ' z'];
            else
                error('Currently only unit vectors for the direction are supported by this function.');
            end
        end

        % Get unique tag for the subplot
        tag = [type ' ' name];

        % Find the subfigure that matches the signal to plot.
        if enableMultipleCalls && nOldSubPlots > 0
            iSubFigTag = find(strcmp(subTags, tag));
            if isempty(iSubFigTag) % If the subfigure does not yet exist, add one!
                iSubFigTag = numel(findobj( get(gcf,'Children'), '-depth', 1, 'type', 'axes')) +1; % Number of axes + 1
            end
        else
            iSubFigTag = iSubFig;
        end

        % Create subfigure
        nRow = ceil(nSubFigNow/nCol);
        if enableMultipleCalls && nOldSubPlots > 0 && iSubFigTag <= nOldSubPlots
            % Use old axes
            axes(subAxes(iSubFigTag));
            hSubFigs{iSubFigTag} = subAxes(iSubFigTag);
        else
            % Create new axes
            hSubFigs{iSubFigTag} = axes('Units', 'centimeters', 'Position', getSubfigurePosition(iSubFigTag, nRow, nCol, subFigSettings), 'Tag', tag);
        end
        hold on;

        % Get min and max of axis
        minY = min([trackAV-sqrt(trackVar);extraAV-sqrt(extraVar); sim-sqrt(simVar)]);
        maxY = max([trackAV+sqrt(trackVar);extraAV+sqrt(extraVar); sim+sqrt(simVar)]);
        margin1 = 0.05*(maxY-minY); % 5 % of range
        margin2 = 1e-3; % for the case that minY and maxY are 0
        margin = max([margin1, margin2]);
        minY = minY - margin;
        maxY = maxY + margin;
        if isempty(minY); minY = 0; end % no data is plotted
        if isempty(maxY); maxY = 1; end
        % if subplot already exists, compare the axis with existing and set
        % max accordingly
        if enableMultipleCalls && nOldSubPlots > 0 && iSubFigTag <= nOldSubPlots
            minYOld = subAxes(iSubFigTag).YLim(1);
            if minY > minYOld
                minY = minYOld;
            end
            maxYOld = subAxes(iSubFigTag).YLim(2);
            if maxY < maxYOld
                maxY = maxYOld;
            end
        end
        axis([0 100 minY maxY]);

        % Plot standing
        if nargin > 2 && ~isempty(standing_r)
            hStandingR = plotStancePhase(time, standing_r, minY, maxY, 0, 100, standingRightColor, standingFaceAlpha);
            allHandels{iSubFig, 7} = hStandingR; % 7th legend entry
        end
        if nargin > 3 && ~isempty(standing_l)
            hStandingL = plotStancePhase(time, standing_l, minY, maxY, 0, 100, standingLeftColor, standingFaceAlpha);
            allHandels{iSubFig, 8} = hStandingL; % 8th legend entry
        end
        
        % Plot tracking
        [hTrackAV, hTrackVar] = plotVariable(time, trackAV, trackVar, trackColor, trackColor, trackFaceAlpha, lineWidth);
        allHandels{iSubFig, 3} = hTrackAV; % 3nd legend entry
        allHandels{iSubFig, 4} = hTrackVar; % 4rd legend entry

        % Plot extra variables
        [hExtraAV, hExtraVar] = plotVariable(time, extraAV, extraVar, extraColor, extraColor, extraFaceAlpha, lineWidth);
        allHandels{iSubFig, 5} = hExtraAV; % 5nd legend entry
        allHandels{iSubFig, 6} = hExtraVar; % 6rd legend entry

        % Plot simulation
        [hSim, hSimVar] = plotVariable(time, sim, simVar, simColor, simColor, trackFaceAlpha, lineWidth);
        allHandels{iSubFig, 1} = hSim; % 1st legend entry
        allHandels{iSubFig, 2} = hSimVar; % 2st legend entry

        % Set labels
        set(gca,'TickLabelInterpreter','latex')
        yTickHandle = get(gca,'YAxis');
        set(yTickHandle,'FontSize', yTickFontSize);
        if ~isempty(unit)
            ylabel([strrep(name, '_', '\textunderscore ') ' in ' strrep(strrep(unit, '%', '\%'),'^2','\textsuperscript{2}')], 'FontSize', yLabelFontSize);
        else
            ylabel(strrep(name, '_', '\textunderscore '), 'FontSize', yLabelFontSize);
        end
        
        xTickHandle= get (gca,'XAxis');
        set (xTickHandle,'FontSize', xTickFontSize);
        xlabel(xLabelText, 'FontSize', xLabelFontSize);

        box on;
    end

    % Get the handles for legend from different subplots (also tracking)
    idxLegNotEmptyAll = ~cellfun('isempty', allHandels);
    idxLegNotEmpty = sum(idxLegNotEmptyAll, 1)>0;
    handels = cell(1, numel(allEntries));
    for iEntry = find(idxLegNotEmpty)
        iSubFig = find(idxLegNotEmptyAll(:, iEntry), 1); % First subplot
        handels{iEntry} = allHandels{iSubFig, iEntry};
    end

    if enableMultipleCalls && nOldSubPlots > 0
        % Get existing legend
        oldLegend = findobj( get(gcf,'Children'), '-depth', 1, 'type', 'legend');

        % Make legend if not all are contained already
        % (It is not combining the legends from different varTables but makes a
        % new one for every table.)
        if ~all(ismember(allEntries(idxLegNotEmpty), oldLegend.String))
            lgnd = addLegend([handels{idxLegNotEmpty}], allEntries(idxLegNotEmpty), allFaceAlphas(idxLegNotEmpty), legendFontSize);
            % Set legend centered below all subfigures
            setLegendBelowSubfigures(lgnd, subFigSettings, nCol)
        end
    else
        % Make legend
        lgnd = addLegend([handels{idxLegNotEmpty}], allEntries(idxLegNotEmpty), allFaceAlphas(idxLegNotEmpty), legendFontSize);
        % Set legend centered below all subfigures
        setLegendBelowSubfigures(lgnd, subFigSettings, nCol)
    end
    
end