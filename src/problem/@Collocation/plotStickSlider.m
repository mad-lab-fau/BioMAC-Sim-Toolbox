%======================================================================
%> @file plotStickSlider.m
%> @brief Collocation function to create plot stick figure with slider
%> @details
%> Details: Collocation::plotStickSlider()
%>
%> @author Marlies Nitschke
%> @date January, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to create plot stick figure with slider
%>
%> @details
%> Uses Model::showStick() to plot the result.
%>
%> If you are adding more content to the plot, it will take noticable more
%> time to initialize the figure. Be patient!
%>
%> Example usage:
%> @code
%> markerTable = result.problem.objectiveTerms(1).varargin{1}.variables; % assuming trackMarker is your first objective
%> result.problem.plotStickSlider(result.X, 0, {1, 0, 1, 0}, markerTable);
%> @endcode
%>
%> @param  obj          Collocation class object
%> @param  X            Double matrix: State vector (i.e. result or initial guess) of the problem
%> @param  subTrans     (optional) Bool: Subtract the horizontal translation
%>                      of the pelvis from the movement. (default: false)
%> @param  showStickOpt (optional) Cell: Cell with optional input for model.showStick. Cell can
%>                      contain all possible input options after ''range''. (default: empty cell)
%> @param  markerTable  (optional) Table: Variables table specifying markers with at least the columns:
%>                      type, name, segment, position, direction to call Model.showMarker. If the
%>                      column "mean" is also give, it will be plotted as reference. (default: empty table)
%======================================================================
function plotStickSlider(obj, X, subTrans, showStickOpt, markerTable)

if nargin < 4
    showStickOpt = {};
end

% get some information
model = obj.model;
isSymmetric = obj.isSymmetric;
N = obj.nNodes;

% get states
x = X(obj.idx.states(:, 1:N))';
[nrows,ncolumns]=size(x);
if ncolumns < model.nDofs
    error('Collocation:plotStickSlider', 'Not enough columns in x.');
end
if N ~= nrows
    error('Collocation:plotStickSlider', 'The motion data does not fit to the nNodes.');
end

% mirrow if symmetric
if isSymmetric
    % get duration and speed of movement
    duration = X(obj.idx.dur);
    speed = X(obj.idx.speed);

    % mirrow
    flip =  x(1:N,:);
    unitdisplacement = zeros(model.nStates,1);
    unitdisplacement(model.idxForward) = 1;
    
    flip = repmat(model.idxSymmetry.xsign,1,N)'.*flip(:,model.idxSymmetry.xindex) + repmat(unitdisplacement,1,N)'*duration*speed;
    x = [x(1:N,:);flip];
    
    N = 2* N;
end

% subtract the horizontal translation of the pelvis from the
% movement
idxPelvisX = model.extractState('q', 'pelvis_tx'); % forwards
idxPelvisZ = model.extractState('q', 'pelvis_tz'); % sidewards
if nargin > 2 && subTrans
    if ~ismember(idxPelvisX, model.idxForward) || (~isempty(idxPelvisZ) && ~ismember(idxPelvisZ, model.idxSideward))
        error('Collocation:plotStickSlider', 'Something wrong with indices. Check states of the model.');
    end
    
    motPelForward = x(:,idxPelvisX); % save in case we need it again
    motPelSideward = x(:,idxPelvisZ);
    x(:, model.idxForward)  = x(:, model.idxForward)  - motPelForward;
    x(:, model.idxSideward) = x(:, model.idxSideward) - motPelSideward;
end

xrange = [min(x(:,idxPelvisX))-1  , max(x(:,idxPelvisX))+1];
yrange = [-0.2, 2];
zrange = [min(x(:,idxPelvisZ))-1  , max(x(:,idxPelvisZ))+1]; % only for 3D otherwise empty
if isempty(zrange); zrange = []; end % change size of emtpy matrix to avoid warning
range = [xrange; yrange; zrange];

nframes = size(x,1);

% prepare marker data
if nargin > 4 && ~isempty(markerTable)
    % Get only the rows with marker data
    markerTable = markerTable(strcmp(markerTable.type, 'marker'), :);

    % get measured mean data and adapt it if we adapted the motion
    if any(strcmp(markerTable.Properties.VariableNames, 'mean'))

        % get data and convert it to meter
        markerMean = cell2mat(markerTable.mean');
        if any(strcmp(markerTable.Properties.VariableNames, 'unit')) && strcmp(markerTable.unit{1}, 'mm')
            markerMean = markerMean / 1000;
        end

        % add NaNs if symmetric to plot only the first half (this is the easiest for now)
        if isSymmetric
            markerMean = [markerMean; nan(size(markerMean))];
        end

        % subtract the horizontal translation of the pelvis from the movement
        if nargin > 2 && subTrans
            % get rows which have forward and sideward motion
            idxForward = find(markerTable.direction(:, 1)); % x direction = [1, 0, 0]
            if size(markerTable.direction, 2) > 2 % for 2D, direction could only be two entries
                idxSideward = find(markerTable.direction(:, 3));
            else
                idxSideward = [];
            end
            % subtract forward motion of pelvis
            markerMean(:, idxForward)  = markerMean(:, idxForward)  - motPelForward;
            markerMean(:, idxSideward) = markerMean(:, idxSideward) - motPelSideward;
        end
    else
        markerMean = nan(nframes, height(markerTable));
    end
end

% initialize figure window
f = figure();
clf;
% f = gcf; hold on;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
set(gcf, 'color', 'white');

% Plot all stick figures
iFrame = 1;
for frame = iFrame : nframes
    hold on;
    
    % Visualize Treadmill speed if it is defined
    if isprop(model, 'speed_right') && isprop(model, 'speed_left') && model.speed_right ~= 0  && model.speed_left ~= 0
        if frame == 1
            [x_points, y_points] = model.showTreadmill([], [], xrange, yrange, X(obj.idx.dur), nframes);
        else
            [x_points, y_points] = model.showTreadmill(x_points, y_points, xrange, yrange, X(obj.idx.dur), nframes);
        end
    end
    
   % Visualize model and marker 
   model.showStick(x(frame,1:model.nStates)', range, showStickOpt{:});
   if nargin > 4 && ~isempty(markerTable)
      model.showMarker(x(frame,1:model.nStates)', markerTable, markerMean(frame, :));
   end
end
nElementsPerFrame = numel(f.Children(1).Children) / nframes;

% Show only first frame
showSingleStick(iFrame);

% Create slider
slidPosBot = 20;
sliderMin = 1;
sliderMax = nframes;
sliderStep = [1 1];
if sliderMax - sliderMin ~= 0
   sliderStep = sliderStep / (sliderMax - sliderMin);
end
b = uicontrol('Parent',f,'Style','slider','Position',[600,slidPosBot,400,23],...
              'value',iFrame, 'min',sliderMin, 'max',sliderMax, 'SliderStep', sliderStep);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[550,slidPosBot,30,23], 'String','1','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[1020,slidPosBot,30,23], 'String',num2str(nframes),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[750,slidPosBot-25,100,23],'String','Frame','BackgroundColor',bgcolor);
bl4 = uicontrol('Parent',f,'Style','text','Position',[750,slidPosBot+25,100,23],'String',num2str(b.Value),'BackgroundColor',bgcolor);
addlistener(b, 'Value', 'PostSet', @(obj,event)set(b, 'Value', max(b.Min,  min(b.Value,b.Max))) ); % create a listener to check the value are reset appropriately
addlistener(b, 'Value', 'PostSet', @(obj,event)set(b, 'Value', round(b.Value)) ); % create a listener to prevent floating values
addlistener(b, 'Value', 'PostSet', @(obj,event)uicontrol('Parent',f,'Style','text','Position',[750,slidPosBot+25,100,23],'String',num2str(b.Value),'BackgroundColor',bgcolor)); % create a listener to set the current frame 
b.Callback = @showSingleStick; % set callback to slider which plots the stick figure

    function showSingleStick(es, ~)
        
        % error checking
        if isnumeric(es)
            iFrame = es;
        elseif isa(es, 'matlab.ui.control.UIControl')
            iFrame = es.Value;
        else
            error('es is of wrong type');
        end
        
        % make all invisible
        set(f.Children(end).Children, 'Visible','off'); 
        
        % show frame iFrame (order of patches is vise reversed)
        set(f.Children(end).Children([1:nElementsPerFrame]+(nframes-iFrame)*nElementsPerFrame), 'Visible','on') 
        
    end

end


