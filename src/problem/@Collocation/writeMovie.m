%======================================================================
%> @file @Collocation/writeMovie.m
%> @brief Collocation function to create a movie of a simulation result
%> @details
%> Details:  Collocation::writeMovie()
%>
%> @author Marlies Nitschke
%> @date January, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to create a movie of a simulation result
%>
%> @details
%> Uses Model::showStick() to plot the result.
%>
%> Whether the optimal solution is converged is not checked.
%>
%> Example usage:
%> @code
%> markerTable = result.problem.objectiveTerms(1).varargin{1}.variables; % assuming trackMarker is your first objective
%> result.problem.writeMovie(result.X, result.filename, [], 0, {1, 0, 1, 0}, 1, 0, markerTable)
%> @endcode
%>
%> @param  obj          Collocation class object
%> @param  X            Double matrix: State vector (i.e. result) of the problem
%> @param  filename     (optional) String: Filename which is used to save the OpenSim files.
%>                      Use an empty char array '' to not save a file.
%> @param  fps          (optional) Double: Desired frame rate (frames per second)
%>                      To skip the input, you can use []. (default: 5)
%> @param  subTrans     (optional) Bool: Subtract the horizontal translation 
%>                      of the pelvis from the movement. To skip the input, 
%>                      you can use [].  (default: false)
%> @param  showStickOpt (optional) Cell: Cell with optional input for model.showStick. Cell can 
%>                      contain all possible input options after ''range''. (default: empty cell)
%> @param  makeGIF      (optional) Bool: If true, a gif will be saved.
%>                      Otherwise an avi will be saved. (default: 0=avi)
%> @param  playSpeed    (Optional) Double: If specified, the real duration of the simulated motion
%>                      will be taken into account with
%>                      - playSpeed == 1: duration of movie reflects the duration of motion
%>                      - playSpeed  < 1: slow motion (e.g. 0.5 => half the playback speed and double the duration)
%>                      - playSpeed  > 1: fast motion (e.g. 2 => double the playback speed and half the duration)
%>                      (default: 0: fps nodes are displayed per second and duration is not taken into account)
%>                      Attention: The duration of the movie will not be equal to the duration of the motion since
%>                      the function is plotting the motion till node nNodes and not till nNodes+1.
%> @param  markerTable  (optional) Table: Variables table specifying markers with at least the columns:
%>                      type, name, segment, position, direction to call Model.showMarker. If the
%>                      column "mean" is also give, it will be plotted as reference. (default: empty table)
%> @param  rangeNodes   (optional) Array: Specifying first and last node
%>                       which should be used for the movie. This might not work well for symmetric motions
%>                       which are automatically mirrored.
%======================================================================
function writeMovie(obj, X, filename, fps, subTrans, showStickOpt, makeGIF, playSpeed, markerTable, rangeNodes)
	
% Error checking
msgID = 'Collocation:writeMovie';
assert(nargin >= 2, msgID, 'Not enough input arguments.');
assert(isa(obj, 'Collocation'), msgID, '''obj'' should be of type ''Collocation''.');
assert(length(X)==obj.nVars && isa(X,'double'), msgID, '''X'' is not of correct size or/and type.');
assert(nargin < 3 || isempty(fileparts(filename)) || exist(fileparts(filename),'dir')==7, msgID, ...
    'The path specified in ''filename'' does not exist');

% Set default entries
if nargin < 3 || isempty(filename)
    saveMovie = 0;
else
    saveMovie = 1;
end
if nargin < 4 || isempty(fps)
    fps = 5; 
end
if nargin < 5 || isempty(subTrans)
    subTrans = 0;
end
if nargin < 6 
    showStickOpt = {};
end
if nargin < 7
   makeGIF = 0; 
end
if nargin < 8
   playSpeed = 0;
end
if nargin < 10
    N = obj.nNodes;
    rangeNodes = 1:N;
else
    N = rangeNodes(2) - rangeNodes(1)+1;
    rangeNodes = rangeNodes(1):rangeNodes(2);
end


% Get some information
model = obj.model;
isSymmetric = obj.isSymmetric;
h = X(obj.idx.dur)/(obj.nNodesDur-1);
duration = h * N; % adapt for rangeNodes
times = 0:h:(N-1)*h;

x = X(obj.idx.states(:, rangeNodes))';
[nrows,ncolumns]=size(x);
assert(ncolumns >= model.nDofs, msgID, 'Not enough columns in x.');
assert(numel(times) == nrows, msgID, 'The time vector does not have the same length as the motion data.');

% Mirrow if symmetric
if isSymmetric
    flip =  x(1:N,:);
    unitdisplacement = zeros(model.nStates,1);
    unitdisplacement(model.idxForward) = 1;
    speed = X(obj.idx.speed);
    flip = repmat(model.idxSymmetry.xsign,1,N)'.*flip(:,model.idxSymmetry.xindex) + repmat(unitdisplacement,1,N)'*duration*speed;
    x = [x(1:N,:);flip];

    N = 2* N;
    duration = 2*duration;
    times = 0:h:(N-1)*h;

end

% Subtract the horizontal translation of the pelvis from the movement
idxPelvisX = model.extractState('q', 'pelvis_tx'); % forwards
idxPelvisZ = model.extractState('q', 'pelvis_tz'); % sidewards
if subTrans
    if ~ismember(idxPelvisX, model.idxForward) || (~isempty(idxPelvisZ) && ~ismember(idxPelvisZ, model.idxSideward))
        error('Collocation:writeMovie', 'Something wrong with indices. Check states of the model.'); 
    end

    motPelForward  = x(:,idxPelvisX); % save in case we need it again
    motPelSideward = x(:,idxPelvisZ);
    x(:, model.idxForward)  = x(:, model.idxForward)  - motPelForward;
    x(:, model.idxSideward) = x(:, model.idxSideward) - motPelSideward;
end

% Get range of movement
xrange = [min(x(:,idxPelvisX))-1  , max(x(:,idxPelvisX))+1];
yrange = [-0.2, 2];
zrange = [min(x(:,idxPelvisZ))-1  , max(x(:,idxPelvisZ))+1]; % only for 3D otherwise empty
if isempty(zrange); zrange = []; end % change size of emtpy matrix to avoid warning
range = [xrange; yrange; zrange];

% Interpolate at different time points if specific playSpeed is required
if playSpeed ~= 0
    nframes = round((duration * fps -1)/ playSpeed); % round it since we can only have integers
    fps = (nframes * playSpeed +1)/duration; % adapt fps to still match the required time even though we rounded nframes => todo: This is not working properly since only integer fps are working. See issue #69.
    timesOld = times;
    times = (0:nframes-1)'/nframes*duration;
    x = interp1(timesOld, x, times, 'linear', 'extrap');

    % Warning: It is not perfectly working yet since only integer fps can be used
    warning('This functionality is unfortunately not yet working perfectly! The movie will have not have 100% the correct duration.');

    if nframes > N
       warning('With %.f fps and a playSpeed of %.02f, you are upsampling the motion from %d nodes to %d frames.', fps, playSpeed, N, nframes);
    end
    if nframes < 10
       warning('With %.f fps and a playSpeed of %.02f, you have only %d frames in total.', fps, playSpeed, nframes);
    end
end

% Prepare marker data
if nargin > 8 && ~isempty(markerTable)
    % Get only the rows with marker data
    markerTable = markerTable(strcmp(markerTable.type, 'marker'), :);

    % get measured mean data and adapt it if we adapted the motion
    if any(strcmp(markerTable.Properties.VariableNames, 'mean'))

        % get data and convert it to meter
        markerMean = cell2mat(markerTable.mean');
        if any(strcmp(markerTable.Properties.VariableNames, 'unit')) && strcmp(markerTable.unit{1}, 'mm')
            markerMean = markerMean / 1000;
        end
        markerMean = markerMean(rangeNodes, :);

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

        % Interpolate at different time points if specific playSpeed is required
        if playSpeed ~= 0
            markerMean = interp1(timesOld, markerMean, times, 'linear', 'extrap');
        end
    else
        markerMean = nan(size(x,1), height(markerTable));
    end
end

% Initialize figure window
figure();
clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.4]);
set(gcf, 'color', 'white');
    
% Plot all frames
nframes = size(x,1);
for i=1:nframes
    
    % Visualize Treadmill speed if it is defined
    if isprop(model, 'speed_right') && isprop(model, 'speed_left') && model.speed_right ~= 0  && model.speed_left ~= 0
        if i == 1
            [x_points, y_points] = model.showTreadmill([], [], xrange, yrange, X(obj.idx.dur), nframes);
        else
            [x_points, y_points] = model.showTreadmill(x_points, y_points, xrange, yrange, X(obj.idx.dur), nframes);
        end
    end
    
    model.showStick(x(i,1:model.nStates)', range, showStickOpt{:});
    if nargin > 8 && ~isempty(markerTable)
        model.showMarker(x(i,1:model.nStates)', markerTable, markerMean(i, :));
    end
    
    if ~saveMovie
        shg; % Show current figure to look at movie
    elseif saveMovie && ~makeGIF
        if (i==1)
            filename = [filename,'.avi'];
            videoWriter = VideoWriter(filename);
            videoWriter.FrameRate = fps;
            open(videoWriter);
            
            F = getframe(gcf);
            frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
        else
            F = getframe(gcf);%,frame);
        end
        writeVideo(videoWriter,F);
    elseif saveMovie && makeGIF
        if (i==1)
            filename = [filename, '.gif'];
            gif(filename, 'DelayTime', 1/fps, 'frame',gca); % use gcf for entire figure
        else
            gif;
        end
    end
    
end
if saveMovie && ~makeGIF
    close(videoWriter);
end
if saveMovie
    disp(['Animation has been saved on ' filename]);
end

end
