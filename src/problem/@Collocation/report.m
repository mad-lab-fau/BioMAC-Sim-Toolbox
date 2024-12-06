%======================================================================
%> @file Collocation/report.m
%> @brief Collocation function to report a current state vector of the problem
%> @details
%> Details: Collocation::report()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Computes function to report a current state vector of the problem
%>
%> @details
%> It calls extractData() (see for details) to obtain the simulated data.
%>
%> If you want to plot a static solution with one node (e.g. standing), you
%> can still compare this to extra data using settings.variableTable. However,
%> the variables table has to be resampled to 100 nodes before calling report().
%>
%> @param  obj             Collocation class object
%> @param  X               Double matrix: State vector (i.e. result) of the problem
%> @param  settings        (optional) Struct: Settings defining what to report. Additionally to
%>                         the fields specified in extractData(), it supports the
%>                         following fields:
%>                          - variableTable     (see extractData())
%>                          - plotObjectives    (boolean, default: 1)
%>                          - plotConstraints   (boolean, default: 1)
%>                          - plotStick         (boolean, default: 1)
%>                          - plotInitialGuess  (boolean, default: 0)
%>                          - plotStance        (boolean, default: 0)
%>                          .
%>                         (use empty to skip)
%> @param  style           (optional) Struct: Settings defining the style for reporting.
%>                         See plotVarType() for details. (use empty to skip)
%> @param  resultFilename  (optional) String: Filename (without extension) to save 
%>                         Matlab figures and to save figures as tikz standalone. 
%> @retval conString       String: Formated text to be printed in the console 
%>                         containing problem informations
%> @retval texString       String: Formated LaTeX text containing problem
%>                         informations and tikz code for the figures.
%> @retval simVarTable     Table: Summarizing all data which was reported
%>                         (see Collocation.extractData() for details)
%======================================================================
function [conString, texString, simVarTable] = report(obj,X,settings,style,resultFilename)

if nargin < 3
    settings = struct();
end
if ~isfield(settings, 'plotObjectives'); settings.plotObjectives = 1; end
if ~isfield(settings, 'plotConstraints'); settings.plotConstraints = 1; end
if ~isfield(settings, 'plotStick'); settings.plotStick = 1; end

if nargin < 4
    style = struct();
end

if nargin > 4
    saveIt = 1;
else
    saveIt = 0;
end

%% Read general information
nNodes = obj.nNodes;
model = obj.model;
isSymmetric = obj.isSymmetric;

%% Add general information to the conString and texString
logicalString = {'false', 'true'};
if ~isfield(obj.idx,'speed') 
    if obj.nNodes == 1
        speed = 0;
    else
        speed = nan; 
    end
    speedStr = sprintf('%i', speed);
else
    speed = X(obj.idx.speed);
    speedStr = sprintf('%.3f, ', speed);
    speedStr = speedStr(1:end-2); % Remove last comma and space
end
if ~isfield(obj.idx,'dur') && obj.nNodes == 1
    dur = 0;
    metCost = 0;
elseif ~isfield(obj.idx,'dur')
    dur = nan;
    metCost = nan;
elseif isa(model, 'quad_11DOF') || isa(model, 'Quadruped')
    dur = X(obj.idx.dur);
    metCost = 1;
else
    dur = X(obj.idx.dur);
    metCost = obj.getMetabolicCost(X);
end 
objTerms = struct2table(obj.objectiveTerms, 'AsArray', 1); % need 'AsArray' if there is only one row
objTerms = objTerms(:, {'name', 'weightedValue', 'weight', 'unweightedValue'});
conTerms = struct2table(obj.constraintTerms, 'AsArray', 1); % need 'AsArray' if there is only one row
conTerms = conTerms(:,  {'name', 'normc'});
if isfield(obj.initialguess, 'info') && isfield(obj.initialguess.info, 'type') % true if makeinitialguess() was called
    initialguessType = obj.initialguess.info.type;
else
    initialguessType = 'Was set in addOptimVar() and not in makeinitialguess()';
end

% get the tracked variables and add the information to objective terms
[trackedVarCon, trackedVarTex] = getTrackedVariables(obj.objectiveTerms, '     -> ');

% Fill the information in a formated string to be printed in
% the console
conString = [newline, ...
             sprintf('**** Problem **** \n'), ...
             sprintf('** General Information ** \n'), ...
             sprintf('   - Model: %s \n', class(model)), ...
             sprintf('   - Number of nodes: %u \n', nNodes), ...
             sprintf('   - Symmetry: %s \n', logicalString{isSymmetric+1}), ...
             sprintf('   - Euler Method: %s \n', obj.Euler), ...
             sprintf('   - Translation speed: %s (m/s) \n', speedStr), ...
             sprintf('   - Movement duration: %.3f (s) \n', dur), ...
             sprintf('   - Metabolic cost: %.3f (J/m/kg) \n', metCost), ...
             sprintf('   - Type of initial guess: %s \n', initialguessType), ...
             sprintf('   - Objective Terms: \n'), ...
             evalc('disp(objTerms)'), ...
             sprintf('     tracked Variables: \n'), ...
             sprintf('%s \n', trackedVarCon), ...
             sprintf('   - Constraint Terms: \n'), ...
             evalc('disp(conTerms)'), ...
             newline];
         
% Modifying the (possibly long) path for Linux and Windows       
initialguessType = strrep(initialguessType, '\', '/'); % Windows would have backslash which cannot be used in LaTeX
initialguessType = strrep(initialguessType, '/', '/\linebreak[3]'); % Breaking a long word if necessary
initialguessType = strrep(initialguessType, '_','\textunderscore '); % Inserting \textunderscore since LaTeX doesn't recognize it 
         
         
% Fill the information in a formated latex string
texString = [newline, ...
             sprintf('\\section{Problem} \n'), ...
             sprintf('\\subsection{General Information} \n'), ...
             sprintf('\\begin{itemize} \n'), ...
             sprintf('   \\item Model: %s \n', strrep(class(model), '_', '\textunderscore ')), ...
             sprintf('   \\item Number of nodes: %u \n', nNodes), ...
             sprintf('   \\item Symmetry: %s \n', logicalString{isSymmetric+1}), ...
             sprintf('   \\item Euler Method: %s \n', obj.Euler), ...
             sprintf('   \\item Translation speed: %s (m/s) \n', speedStr), ...
             sprintf('   \\item Movement duration: %.3f (s) \n', dur), ...
             sprintf('   \\item Metabolic cost: %.3f (J/m/kg) \n', metCost), ...
             sprintf('   \\item Type of initial guess: %s \n', initialguessType), ...
             sprintf('   \\item Objective Terms: \n') ...
             tableToLaTeX(objTerms, {'%s','%e','%e','%e'}), ...
             sprintf('tracked Variables: \n') ...
             sprintf('%s \n', trackedVarTex), ...
             sprintf('   \\item Constraint Terms: \n'), ...
             tableToLaTeX(conTerms, {'%s', '%e'}), ...
             sprintf('\\end{itemize} \n')];

%% Get simulated data
if isfield(settings, 'plotStance') && settings.plotStance
    settings.standing = {'standing_r', 'standing_l'}; % Do always both
end
if ~isfield(settings, 'variableTable')
   settings.variableTable = [];
end
simVarTable = obj.extractData(X, settings, settings.variableTable, 1);

%% Do plotting, saving, strings for all fields
% Prepare some stuff
if isfield(settings, 'plotStance') && settings.plotStance
    doPlotStancePhase  = 1;
    standing_r = simVarTable.sim{strcmp(simVarTable.type, 'standing') & strcmp(simVarTable.name, 'standing_r')};
    if ~isSymmetric
        standing_l = simVarTable.sim{strcmp(simVarTable.type, 'standing') & strcmp(simVarTable.name, 'standing_l')};
    else
        standing_l = [];
    end
else
    doPlotStancePhase  = 0;
end

% Do it for all types
excludedTypes = {'standing', 'standingEvent', 'footAngleEvent', 'dur', 'speed'};
types = setdiff(unique(simVarTable.type), excludedTypes);
for iType = 1 : length(types)
    
    type = types{iType};
    
    % Get only entries of type
    simVarTableType = simVarTable(strcmp(simVarTable.type, type), :);
    
    % Make the figure
    if doPlotStancePhase
        [hFig, style] = plotVarType(simVarTableType, style, standing_r, standing_l);
    else
        [hFig, style] = plotVarType(simVarTableType, style);
    end
    
    % Save figure
    if saveIt % standalone
        texStringFigure = saveStandaloneFig(hFig, [resultFilename '_' type]);
    else % no standalone
        texStringFigure = getTikzCode(hFig);
    end
    
    % Create tex string for the subsection of the current type
    texString = [texString, newline, ...
        sprintf('\\subsection{%s} \n', type), ...
        sprintf('\\begin{figure}[H] \n'), ...
        sprintf('\\centering \n'),...
        texStringFigure, newline, ...
        sprintf('\\end{figure} \n')];
    
end

%% Get histories of objective terms and plot them
fieldName = 'plotObjectives';
if isfield(settings, fieldName) && ~isempty(settings.(fieldName)) && settings.(fieldName)
    
    typeStr = 'Objectives';
    
    % Get all objectives which are there
    objTab = struct2table(obj.objectiveTerms, 'AsArray', 1); % need 'AsArray' if there is only one row
    
    % Make the figure
    hFig = plotObjTable(objTab, style);
    
    % Save figure
    if saveIt % standalone
        filename = [resultFilename '_' typeStr];
        % Save .fig
        savefig(hFig, [filename, '.fig']);
        % Save tikz (standalone file)
        matlab2tikz('figurehandle',hFig,'filename',[filename '.tex'] ,'standalone', true, 'showInfo', false, 'strictFontSize',true);
        % Save the pdf
        [onlyFolderPath,onlyFileName]=fileparts(filename);
        % check whether the tex file has no errors so the pdf can be safely created
        try
            callPdflatex(onlyFolderPath,onlyFileName);
        catch
            warning('Could not make the pdf since the associated figure file has errors. Check the file:\n %s.tex\n', filename);
        end
        % Get line to include the standalone file
        texStringFigure =  sprintf('\\includestandalone {%s} ',onlyFileName);
    else % no standalone
        % Write tikz tex string without saving it (no standalone)
        texStringFigure = matlab2tikz('figurehandle',hFig,'standalone', false, 'showInfo', false, 'strictFontSize',true); % no filename given = not saving
    end
    % Create tex string for "Objective Terms" subsection
    texString = [texString, newline, ...
        sprintf('\\subsection{Objective Terms} \n'), ...
        sprintf('\\begin{figure}[H] \n'), ...
        sprintf('\\centering \n'),...
        texStringFigure, newline, ...
        sprintf('\\end{figure} \n')];

end

%% Get histories of objective terms and plot them
fieldName = 'plotConstraints';
if isfield(settings, fieldName) && ~isempty(settings.(fieldName)) && settings.(fieldName)
    
    typeStr = 'Constraints';
    
    % Get all constraint which are there
    conTab = struct2table(obj.constraintTerms, 'AsArray', 1); % need 'AsArray' if there is only one row
    
    % Make the figure
    hFig = plotConTable(conTab, style);
    
    % Save figure
    if saveIt % standalone
        filename = [resultFilename '_' typeStr];
        % Save .fig
        savefig(hFig, [filename, '.fig']);
        % Save tikz (standalone file)
        matlab2tikz('figurehandle',hFig,'filename',[filename '.tex'] ,'standalone', true, 'showInfo', false, 'strictFontSize',true);
        % Save the pdf
        [onlyFolderPath,onlyFileName]=fileparts(filename);
        % check whether the tex file has no errors so the pdf can be safely created
        try
            callPdflatex(onlyFolderPath,onlyFileName);
        catch
            warning('Could not make the pdf since the associated figure file has errors. Check the file:\n %s.tex\n', filename);
        end
        % Get line to include the standalone file
        texStringFigure =  sprintf('\\includestandalone {%s} ',onlyFileName);
    else % no standalone
        % Write tikz tex string without saving it (no standalone)
        texStringFigure = matlab2tikz('figurehandle',hFig,'standalone', false, 'showInfo', false, 'strictFontSize',true); % no filename given = not saving
    end
    % Create tex string for "Constraint Terms" subsection
    texString = [texString, newline, ...
        sprintf('\\subsection{Constraint Terms} \n'), ...
        sprintf('\\begin{figure}[H] \n'), ...
        sprintf('\\centering \n'),...
        texStringFigure, newline, ...
        sprintf('\\end{figure} \n')];

end

%% Plot the stick figure of the result
fieldName = 'plotStick';
if isfield(settings, fieldName) && ~isempty(settings.(fieldName)) && settings.(fieldName)
    
    typeStr = 'StickFigure';
    
    % Make the figure
    hFig = plotStickFigures(obj.model, X(obj.idx.states), style, typeStr);
       
    % Save figure
    if saveIt % standalone
        filename = [resultFilename '_' typeStr];
        % Save .fig
        savefig(hFig, [filename, '.fig']);
        % Save tikz (standalone file)
        matlab2tikz('figurehandle',hFig,'filename',[filename '.tex'] ,'standalone', true, 'showInfo', false, 'strictFontSize',true);
        % Save the pdf
        [onlyFolderPath,onlyFileName]=fileparts(filename);
        % check whether the tex file has no errors so the pdf can be safely created
        try
            callPdflatex(onlyFolderPath,onlyFileName);
        catch
            warning('Could not make the pdf since the associated figure file has errors. Check the file:\n %s.tex\n', filename);
        end
        % Get line to include the standalone file
        texStringFigure =  sprintf('\\includestandalone {%s} ',onlyFileName);
    else % no standalone
        % Write tikz tex string without saving it (no standalone)
        texStringFigure = matlab2tikz('figurehandle',hFig,'standalone', false, 'showInfo', false, 'strictFontSize',true); % no filename given = not saving
    end
    % Create tex string for "Stick Figure" subsection
    texString = [texString, newline, ...
        sprintf('\\subsection{Stick Figure} \n'), ...
        sprintf('\\begin{figure}[H] \n'), ...
        sprintf('\\centering \n'),...
        texStringFigure, newline, ...
        sprintf('\\end{figure} \n')];
    
end


%% Plot the stick figure of the initial guess
fieldName = 'plotInitialGuess';
if isfield(settings, fieldName) && ~isempty(settings.(fieldName)) && settings.(fieldName)
    
    typeStr = 'StickFigureInitialGuess';
    
    % Make the figure
    hFig = plotStickFigures(obj.model, obj.initialguess.X(obj.idx.states), style, typeStr);
        
    % Save figure
    if saveIt % standalone
        filename = [resultFilename '_' typeStr];
        % Save .fig
        savefig(hFig, [filename, '.fig']);
        % Save tikz (standalone file)
        matlab2tikz('figurehandle',hFig,'filename',[filename '.tex'] ,'standalone', true, 'showInfo', false, 'strictFontSize',true);
        % Save the pdf
        [onlyFolderPath,onlyFileName]=fileparts(filename);
        % check whether the tex file has no errors so the pdf can be safely created
        try
            callPdflatex(onlyFolderPath,onlyFileName);
        catch
            warning('Could not make the pdf since the associated figure file has errors. Check the file:\n %s.tex\n', filename);
        end
        % Get line to include the standalone file
        texStringFigure =  sprintf('\\includestandalone {%s} ',onlyFileName);
    else % no standalone
        % Write tikz tex string without saving it (no standalone)
        texStringFigure = matlab2tikz('figurehandle',hFig,'standalone', false, 'showInfo', false, 'strictFontSize',true); % no filename given = not saving
    end
    % Create tex string for "Stick Figure of Initial Guess" subsection
    texString = [texString, newline, ...
        sprintf('\\subsection{Stick Figure of Initial Guess} \n'), ...
        sprintf('\\begin{figure}[H] \n'), ...
        sprintf('\\centering \n'),...
        texStringFigure, newline, ...
        sprintf('\\end{figure} \n')];
    
end




end


%> @cond DO_NOT_DOCUMENT  
%======================================================================
%> @brief Function to make the whole figure for objective terms
%>
%> @param   objTable       Table: Result of struct2table(problem.objectiveTerms).
%>                         The table has to contain the fields: 
%>                         name, weightedValueHist, weight.
%> @param   style          Struct: Style of the figure containing the fields:
%>                           - figureSize
%>                           - xLabelFontSize
%>                           - yLabelFontSize
%>                           - xTickFontSize
%>                           - yTickFontSize
%>                           - legendFontSize
%>                           - lineWidth
%>                           - subFigSettings
%>                         See plotVarType() for details.
%> @retval  hFig           Handle of figure
%======================================================================
function hFig = plotObjTable(objTable, style)

    % Set values
    type = 'Objectives';
    nRow = 1; % We plot only two subplots
    nCol = 2;
    nTerms = height(objTable);
    allEntries = cell(nTerms, 1); % legend entries
    allHandels = cell(nTerms, 1); % handles for legend
    
    % Get stuff from style struct
    % Figure size
    figureSize         = style.figureSize;
    % X Label Font size
    xLabelFontSize         = style.xLabelFontSize;
    % Y Label Font size
    yLabelFontSize         = style.yLabelFontSize;
    % X Tick Font Size
    xTickFontSize           = style.xTickFontSize ;
    % Y Tick Font Size
    yTickFontSize           = style.yTickFontSize ;
    % Legend Font Size
    legendFontSize         = style.legendFontSize;  
    % Line Width
    lineWidth          = style.lineWidth;
    % Subfigure size and positioning
    subFigSettings     = style.subFigSettings;

    
    % Initialize figure
    hFig = figure('name', type);
    clf;
    set(gcf,'units','centimeters','outerposition',figureSize);
    set(gcf,'defaultTextInterpreter','latex');
    
    % Plot the weighted terms
    hSubFig1 = axes('Units', 'centimeters', 'Position', getSubfigurePosition(1, nRow, nCol, subFigSettings));
    hold on;
    for iTerm = 1 : nTerms   
        allHandels{iTerm} = plot(objTable.weightedValueHist{iTerm}, 'LineWidth', lineWidth);
        allEntries{iTerm} = objTable.name{iTerm};
    end
    
    % Set ticks
    xTickHandle= get (gca,'XAxis');
    set (xTickHandle,'FontSize', xTickFontSize );
    yTickHandle= get (gca,'YAxis');
    set (yTickHandle,'FontSize', yTickFontSize);
    
    % Set labels
    set(gca,'TickLabelInterpreter','latex')
    ylabel('Weighted Terms', 'FontSize', yLabelFontSize);
    xlabel('Function Evaluations', 'FontSize', xLabelFontSize);
    
    
    box on;
    
    % Plot unweighted terms 
    % => These are not saved. However, as long as the weight of a term is
    % not 0, we can compute it.
    hSubFig2 = axes('Units', 'centimeters', 'Position', getSubfigurePosition(2, nRow, nCol, subFigSettings));
    hold on;
    for iTerm = 1 : nTerms  
        if objTable.weight(iTerm) ~= 0
            plot(objTable.weightedValueHist{iTerm}/objTable.weight(iTerm), 'LineWidth', lineWidth);
        end
    end
    % Set ticks
    xTickHandle= get (gca,'XAxis');
    set (xTickHandle,'FontSize', xTickFontSize);
    yTickHandle= get (gca,'YAxis');
    set (yTickHandle,'FontSize', yTickFontSize);
    
    % Set labels
    set(gca,'TickLabelInterpreter','latex')
    ylabel('Unweighted Terms', 'FontSize', yLabelFontSize);
    xlabel('Function Evaluations', 'FontSize', xLabelFontSize);
    
    
    box on;
    
    % Get subplot containing the most data (also tracking)
    lgnd = addLegend([allHandels{:}], allEntries,legendFontSize);
    % Set legend centered below all subfigures
    setLegendBelowSubfigures(lgnd, subFigSettings, nCol)
    
end


%======================================================================
%> @brief Function to make the whole figure for constraint terms
%>
%> @param   objTable       Table: Result of struct2table(problem.constraintTerms).
%>                         The table has to contain the fields: name, normcHist
%> @param   style          Struct: Style of the figure containing the fields:
%>                           - figureSize
%>                           - xLabelFontSize
%>                           - yLabelFontSize
%>                           - xTickFontSize
%>                           - yTickFontSize
%>                           - legendFontSize
%>                           - lineWidth
%>                           - subFigSettings
%>                           - simColor
%>                         See plotVarType() for details.
%> @retval  hFig           Handle of figure
%======================================================================
function hFig = plotConTable(conTable, style)

    % Set values
    type = 'Constraints';
    nRow = 1; % We plot only one subplots
    nCol = 1;
    nTerms = height(conTable);
    allEntries = cell(nTerms, 1); % legend entries
    allHandels = cell(nTerms, 1); % handles for legend
    
    % Get stuff from style struct
    % Figure size
    figureSize         = style.figureSize;
    % X Label Font size
    xLabelFontSize         = style.xLabelFontSize;
    % Y Label Font size
    yLabelFontSize         = style.yLabelFontSize;
    % X Tick Font Size
    xTickFontSize           = style.xTickFontSize ;
    % Y Tick Font Size
    yTickFontSize           = style.yTickFontSize ;
    % Legend Font Size
    legendFontSize         = style.legendFontSize;  
    % Line Width
    lineWidth          = style.lineWidth;
    % Subfigure size and positioning
    subFigSettings     = style.subFigSettings;

    
    % Initialize figure
    hFig = figure('name', type);
    clf;
    set(gcf,'units','centimeters','outerposition',figureSize);
    set(gcf,'defaultTextInterpreter','latex');
    
    % Plot the weighted terms
    hSubFig1 = axes('Units', 'centimeters', 'Position', getSubfigurePosition(1, nRow, nCol, subFigSettings));
    hold on;
    for iTerm = 1 : nTerms   
        allHandels{iTerm} = plot(conTable.normcHist{iTerm}, 'LineWidth', lineWidth);
        allEntries{iTerm} = conTable.name{iTerm};
    end
    % Set ticks
    xTickHandle= get (gca,'XAxis');
    set (xTickHandle,'FontSize', xTickFontSize);
    yTickHandle= get (gca,'YAxis');
    set (yTickHandle,'FontSize', yTickFontSize);
    % Set labels
    set(gca,'TickLabelInterpreter','latex')
    ylabel('Constraint Violations', 'FontSize', yLabelFontSize);
    xlabel('Function Evaluations', 'FontSize', xLabelFontSize);
    
    
    box on;
    
    % Get subplot containing the most data (also tracking)
    lgnd = addLegend([allHandels{:}], allEntries,legendFontSize);
    % Set legend centered below all subfigures
    setLegendBelowSubfigures(lgnd, subFigSettings, nCol)
    
end


%======================================================================
%> @brief Function to plot the stick figures from different views
%>
%> @param   model      
%> @param   x          Double matrice: State vector of model for n time points (Gait3d.nStates x n)
%> @param   style      Struct: Style of the figure containing the fields:
%>                       - figureSize
%>                       - xLabelFontSize
%>                       - yLabelFontSize
%>                       - xTickFontSize
%>                       - yTickFontSize
%>                       - legendFontSize
%>                       - lineWidth
%>                       - subFigSettings
%>                       - simColor
%>                      See plotVarType() for details.
%> @param   typeStr     String: Name of the figure
%> @retval  hFig        Handle of figure
%======================================================================
function hFig = plotStickFigures(model, x, style, typeStr)

    % Initialize figure
    hFig = figure('name', typeStr);
    clf;
    set(gcf,'units','centimeters','outerposition',style.figureSize);
    set(gcf,'defaultTextInterpreter','latex');
    
    % Get specific subfigure settings
    subFigSettings = style.subFigSettings;
    subFigSettings.width = style.figureSize(3);
    subFigSettings.height = style.figureSize(4)/3 - 2*subFigSettings.oneUp; % 3 is here maximum number of subplots below each other
        
    % Make the plot
    if isa(model,'Gait2dc')
        
        % Plot it
        h1 = subplot(3, 1, 1); % Use also subplot to have same scaling as for 3D
        set(h1, 'Units', 'centimeters');
        model.showStick(x); axis normal;
        
        % Edit tick and label font sizes
        set(get(gca,'XAxis'),  'FontSize', style.xTickFontSize);
        set(get(gca,'YAxis'),  'FontSize', style.yTickFontSize);
        set(get(gca,'XLabel'), 'FontSize', style.xLabelFontSize);
        set(get(gca,'YLabel'), 'FontSize', style.yLabelFontSize);        
        
        % Get ranges of plots in meter 
        xlimits = xlim;
        ylimits = ylim;
        sizesMeter(1) = xlimits(2) - xlimits(1); % used as width
        sizesMeter(2) = ylimits(2) - ylimits(1); % used as height
        
        % Get sizes for plots
        idxWidth = 1;
        [~, idxMaxWidth] = max(sizesMeter(idxWidth));
        scaleWidth = subFigSettings.width/sizesMeter(idxWidth(idxMaxWidth));
        
        idxHeight = 2;
        [~, idxMaxHeight] = max(sizesMeter(idxHeight));
        scaleHeight = subFigSettings.height/sizesMeter(idxHeight(idxMaxHeight));
        
        scale = min(scaleWidth, scaleHeight); % That's the one fitting in the figure
        
        h1SubFigPos = h1.Position;
        pos = getSubfigurePosition(1, 3, 1, subFigSettings);
        h1SubFigPos(1) = pos(1);
        h1SubFigPos(2) = pos(2); 
        h1SubFigPos(3) = sizesMeter(1)*scale; % width
        h1SubFigPos(4) = sizesMeter(2)*scale; % height
        set(h1, 'Position', h1SubFigPos);
        
    elseif isa(model, 'Gait3d')
        
        % Plot front view
        h1 = subplot(3, 1, 1);
        set(h1, 'Units', 'centimeters');
        model.showStick(x,[],0,0,0,0,0,0); axis normal;
        
        % Edit tick and label font sizes
        set(get(gca,'XAxis'),  'FontSize', style.xTickFontSize);
        set(get(gca,'YAxis'),  'FontSize', style.yTickFontSize);
        set(get(gca,'ZAxis'),  'FontSize', style.yTickFontSize);  % Use the style of y also for z
        set(get(gca,'XLabel'), 'FontSize', style.xLabelFontSize);
        set(get(gca,'YLabel'), 'FontSize', style.yLabelFontSize);
        set(get(gca,'ZLabel'), 'FontSize', style.yLabelFontSize); % Use the style of y also for z
        
        % Plot side view
        h2 = subplot(3, 1, 2);
        set(h2, 'Units', 'centimeters');
        model.showStick(x); axis normal;
        
        % Edit tick and label font sizes
        set(get(gca,'XAxis'),  'FontSize', style.xTickFontSize);
        set(get(gca,'YAxis'),  'FontSize', style.yTickFontSize);
        set(get(gca,'ZAxis'),  'FontSize', style.yTickFontSize);  % Use the style of y also for z
        set(get(gca,'XLabel'), 'FontSize', style.xLabelFontSize);
        set(get(gca,'YLabel'), 'FontSize', style.yLabelFontSize);
        set(get(gca,'ZLabel'), 'FontSize', style.yLabelFontSize); % Use the style of y also for z
        
        % Plot top view
        h3 = subplot(3, 1, 3);
        set(h3, 'Units', 'centimeters');
        model.showStick(x,[],0,0,0,0,90,90); axis normal;
        
        % Edit tick and label font sizes
        set(get(gca,'XAxis'),  'FontSize', style.xTickFontSize);
        set(get(gca,'YAxis'),  'FontSize', style.yTickFontSize);
        set(get(gca,'ZAxis'),  'FontSize', style.yTickFontSize);  % Use the style of y also for z
        set(get(gca,'XLabel'), 'FontSize', style.xLabelFontSize);
        set(get(gca,'YLabel'), 'FontSize', style.yLabelFontSize);
        set(get(gca,'ZLabel'), 'FontSize', style.yLabelFontSize); % Use the style of y also for z 
        
        % Get ranges of plots in meter (attention: labels to not fit to axes)
        xlimits = ylim;
        ylimits = zlim;
        zlimits = xlim;
        sizesMeter(1) = xlimits(2) - xlimits(1); % used as width
        sizesMeter(2) = ylimits(2) - ylimits(1); % used as height
        sizesMeter(3) = zlimits(2) - zlimits(1); % used as width and height
        
        % Get sizes for plots
        idxWidth = [1, 3];
        [~, idxMaxWidth] = max(sizesMeter(idxWidth));
        scaleWidth = subFigSettings.width/sizesMeter(idxWidth(idxMaxWidth));
        
        idxHeight = [2, 3];
        [~, idxMaxHeight] = max(sizesMeter(idxHeight));
        scaleHeight = subFigSettings.height/sizesMeter(idxHeight(idxMaxHeight));
        
        scale = min(scaleWidth, scaleHeight); % That's the one fitting in the figure
        
        h1SubFigPos = h1.Position;
        pos = getSubfigurePosition(1, 3, 1, subFigSettings);
        h1SubFigPos(1) = pos(1)+(sizesMeter(1)-sizesMeter(3))/2*scale; % center the figure
        h1SubFigPos(2) = pos(2); 
        h1SubFigPos(3) = sizesMeter(3)*scale; % width
        h1SubFigPos(4) = sizesMeter(2)*scale; % height
        set(h1, 'Position', h1SubFigPos);
        
        h2SubFigPos = h2.Position;
        pos = getSubfigurePosition(2, 3, 1, subFigSettings);
        h2SubFigPos(1) = pos(1); 
        h2SubFigPos(2) = pos(2); 
        h2SubFigPos(3) = sizesMeter(1)*scale; % width
        h2SubFigPos(4) = sizesMeter(2)*scale; % height
        set(h2, 'Position', h2SubFigPos);
        
        h3SubFigPos = h3.Position;
        pos = getSubfigurePosition(3, 3, 1, subFigSettings);
        h3SubFigPos(1) = pos(1); 
        h3SubFigPos(2) = pos(2);
        h3SubFigPos(3) = sizesMeter(1)*scale; % width
        h3SubFigPos(4) = sizesMeter(3)*scale; % height
        set(h3, 'Position', h3SubFigPos);
        
    elseif isa(model,'Quadruped')
        % Plot it
        h1 = subplot(3, 1, 1); % Use also subplot to have same scaling as for 3D
        set(h1, 'Units', 'centimeters');
        model.showStick(x); axis normal;
        
        % Edit tick and label font sizes
        set(get(gca,'XAxis'),  'FontSize', style.xTickFontSize);
        set(get(gca,'YAxis'),  'FontSize', style.yTickFontSize);
        set(get(gca,'XLabel'), 'FontSize', style.xLabelFontSize);
        set(get(gca,'YLabel'), 'FontSize', style.yLabelFontSize);        
        
        % Get ranges of plots in meter 
        xlimits = xlim;
        ylimits = ylim;
        sizesMeter(1) = xlimits(2) - xlimits(1); % used as width
        sizesMeter(2) = ylimits(2) - ylimits(1); % used as height
        
        % Get sizes for plots
        idxWidth = 1;
        [~, idxMaxWidth] = max(sizesMeter(idxWidth));
        scaleWidth = subFigSettings.width/sizesMeter(idxWidth(idxMaxWidth));
        
        idxHeight = 2;
        [~, idxMaxHeight] = max(sizesMeter(idxHeight));
        scaleHeight = subFigSettings.height/sizesMeter(idxHeight(idxMaxHeight));
        
        scale = min(scaleWidth, scaleHeight); % That's the one fitting in the figure
        
        h1SubFigPos = h1.Position;
        pos = getSubfigurePosition(1, 3, 1, subFigSettings);
        h1SubFigPos(1) = pos(1);
        h1SubFigPos(2) = pos(2); 
        h1SubFigPos(3) = sizesMeter(1)*scale; % width
        h1SubFigPos(4) = sizesMeter(2)*scale; % height
        set(h1, 'Position', h1SubFigPos);
    else
        warning('The class ''%s'' of the model is unknown. Only models of type ''Gait2dc'' or ''Gait3d'' or classes inheriting from them can be used.', class(model));
    end
    

    
end


%======================================================================
%> @brief Function to extract information about tracked variables
%>
%> @details
%> It extracts all tracked variables for each objective and returns them as
%> String for visualization in the console and for passing them to latex
%>
%> @param   objectiveTerms    Struct: containing the measured objectives and all further details
%> @param   prefixCon         String: Prefix used at beginning of every line for console format
%> @retval  trackedVarCon     String: containing all tracked variables in console format
%> @retval  trackedVarTex     String: containing all tracked variabled in tex format
%======================================================================
function [trackedVarCon, trackedVarTex] = getTrackedVariables(objectiveTerms, prefixCon)

% Get indices of objectives which contain "track"
idx = find(contains({objectiveTerms.name}, 'track'));

% Initialize
nTrackingTerms = numel(idx);
trackedVar = cell(nTrackingTerms, 1);
iTrack = 0;

for iTerm = idx
    iTrack = iTrack +1;
    
    % delete 'track' from objectives name for nicer visualization
    if strfind(objectiveTerms(iTerm).name, 'track')
        trackedVar{iTrack} = cat(2, trackedVar{iTrack}, objectiveTerms(iTerm).name(6:end), ': ');
    else
        trackedVar{iTrack} = cat(2, trackedVar{iTrack}, objectiveTerms(iTerm).name, ': ');
    end
    
    % get the variables as a cell array
    c = sort(table2cell(objectiveTerms(iTerm).varargin{1,1}.variables(:,{'name'})));
    
    if ~isempty(c{1})
        % check for variables that were tracked on both sides
        sym = arrayfun(@(x) (c{x}(end) == 'l' || c{x}(end) == 'r') ... % ends with l or r
            && sum(contains(c, c{x}(1:end-1)))==2, ... % and there are two elements starting with this name
            1:numel(c)); % for all elements of c

        % if the variable is tracked on both sides replace the side by
        % (left and right)
        for k = 1:numel(sym)
            if sym(k)
                c{k} = cat(2,c{k}(1:end-2),' (left and right)');
            end
        end
        c = unique(c); % delete duplicate entries

        % add the variable to the string
        for k = 1:size(c,1)
            trackedVar{iTrack} = cat(2,trackedVar{iTrack},c{k});
            if k < size(c,1)
                trackedVar{iTrack} = cat(2,trackedVar{iTrack},', ');
            end
        end
    end
    
end

% Get output for console
trackedVarCon = [];
if nTrackingTerms == 0 % no tracking
    trackedVarCon = '     No tracking';
else
    for iTrack = 1 : nTrackingTerms
        trackedVarCon = [trackedVarCon, prefixCon, trackedVar{iTrack}];
        if iTrack ~= nTrackingTerms
            trackedVarCon = [trackedVarCon, sprintf(' \n')];
        end
    end
end

% Get output for latex
if nTrackingTerms == 0 % no tracking
    trackedVarTex = 'No tracking';
else
    trackedVarTex = ['\begin{itemize}' sprintf(' \n')];
    for iTrack = 1 : nTrackingTerms
        trackedVarTex = [trackedVarTex, '\item ', strrep(trackedVar{iTrack}, '_','\textunderscore '), sprintf(' \n')];
    end
    trackedVarTex = [trackedVarTex, '\end{itemize}'];
end

end

%======================================================================
%> @brief Function to get figure as tikz code
%>
%> @details
%> It writes tex code as string without saving the figure as standalone tikz.
%>
%> @param   hFig             Handle: Figure handle from which the tikz code
%>                           should be generated.
%> @retval  texStringFigure  String: Tex code of the figure.
%======================================================================
function texStringFigure = getTikzCode(hFig)

texStringFigure = matlab2tikz('figurehandle',hFig,'standalone', false, 'showInfo', false, 'strictFontSize',true); % no filename given = not saving

end

%======================================================================
%> @brief Function to save figure as standalone tikz file
%>
%> @details
%> It saves the figure as .fig and standalone tikz file. It directly
%> creates a pdf from the tikz code.
%>
%> @param   hFig             Handle: Figure handle from which the tikz code
%>                           should be generated.
%> @param   filename         String: Filename including path to save the figure
%> @retval  texStringFigure  String: Tex code to include the standalone tex
%>                           file into the overall tex document.
%======================================================================
function texStringFigure = saveStandaloneFig(hFig, filename)

% Save .fig
savefig(hFig, [filename, '.fig']);

% Save tikz (standalone file)
matlab2tikz('figurehandle',hFig,'filename',[filename '.tex'] ,'standalone', true, 'showInfo', false, 'strictFontSize',true);
% Save the pdf
[onlyFolderPath,onlyFileName]=fileparts(filename);
% check whether the tex file has no errors so the pdf can be safely created
try
    callPdflatex(onlyFolderPath,onlyFileName);
catch
    warning('Could not make the pdf since the associated figure file has errors. Check the file:\n %s.tex\n', filename);
end

% Get line to include the standalone file
texStringFigure =  sprintf('\\includestandalone {%s} ',onlyFileName);

end

%> @endcond
