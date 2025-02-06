%======================================================================
%> @file @Collocation/writeMovieCurvedRunning.m
%> @brief Collocation function to create a movie of a simulation result for curved running
%> @details
%> Details:  Collocation::writeMovieCurvedRunning()
%>
%> @author Eva Dorschky, Marlies Nitschke
%> @date January, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to create a movie of a simulation result for curved running
%>
%> @details
%> This function writes movies for simulations which used
%> Collocation::curveConstraint_pelvis213() instead of Collocation::periodicityConstraint().
%>
%> Whether the optimal solution is converged is not checked.
%>
%> @param  obj          Collocation class object
%> @param  X            Double matrix: State vector (i.e. result) of the problem
%> @param  filename     String: Filename which is used to save the OpenSim files
%> @param  radius       Double: Radius in m
%> @param  fps          (optional) Double: Desired frame rate (frames per second)
%>                      To skip the input, you can use []. (default: 5)
%> @param  showStickOpt (optional) Cell: Cell with optional input for model.showStick. Cell can
%>                      contain all possible input options after ''range''. (default: empty cell)
%> @param  makeGIF      (optional) Bool: If true, a gif will be saved.
%>                      Otherwise an avi will be saved. (default: 0=avi)
%======================================================================
function writeMovieCurvedRunning(obj, X, filename, radius, fps, showStickOpt, makeGIF)

% Error checking
msgID = 'Collocation:writeMovieCurvedRunning';
assert(nargin >= 4, msgID, 'Not enough input arguments.');
assert(isa(obj, 'Collocation'), msgID, '''obj'' should be of type ''Collocation''.');
assert(length(X)==obj.nVars && isa(X,'double'), msgID, '''X'' is not of correct size or/and type.');
assert(~isempty(filename), msgID, '''filename'' is empty.');
assert(exist(fileparts(filename),'dir')==7, msgID, 'The path specified in ''filename'' does not exist.');

% Set default entries
if nargin < 5 || isempty(fps)
    fps = 5;
end
if nargin < 6 
   showStickOpt = {};
end
if nargin < 7
   makeGIF = 0; 
end

% Get some information
model = obj.model;
isSymmetric = obj.isSymmetric;
N = obj.nNodes;
duration = X(obj.idx.dur);
times = (0:N-1)'/N*duration;
speed = X(obj.idx.speed);

x = X(obj.idx.states(:, 1:N))';
[nrows,ncolumns]=size(x);
assert(ncolumns >= model.nDofs, msgID, 'Not enough columns in x.');
assert(numel(times) == nrows, msgID, 'The time vector does not have the same length as the motion data.');

% Mirrow if symmetric => Currently not supported!
if isSymmetric
    error('Symmetric results are not supported. Please add the functionality.');
end

% Add gait cycles to get a whole circle
% all states in forward direction (x direction)
idxForwardAll = model.idxForwardAll;
% all states in sideward direction (x direction);
idxSidewardAll = model.idxSidewardAll;
% pelvis rotation
idxPelvisRotation = obj.idx.states(model.extractState('q','pelvis_rotation'));
% angle of circular segment
theta = 2*asin(norm(speed)*duration/(2*radius));
% number of gait need for a whole circle
nTimes = ceil(2*pi/abs(theta));
% initialize
x_firCycle = x; % first cycle
x_curCycle = x; % current cylce
x_trans = zeros(size(x)); % horizontal translation (no pelvis rotation)
% compute horizontal shift during a gait cycle
x_shift = zeros(size(x)); 
x_shift(:, idxForwardAll)  = x_firCycle(:,idxForwardAll)  -repmat(x_firCycle(1,idxForwardAll),N,1);
x_shift(:, idxSidewardAll) = x_firCycle(:,idxSidewardAll) -repmat(x_firCycle(1,idxSidewardAll),N,1);
% add nTimes-1 gait cycles
for iTimes = 1:nTimes-1
   
    % current total angle
    thetaTotal = iTimes*theta;
  
    % compute global translation
    x_tmp_end = zeros(size(x(1, :)));
    x_tmp_end(idxForwardAll)  =  x_firCycle(1, idxForwardAll) * cos(thetaTotal) + x_firCycle(1, idxSidewardAll) * sin(thetaTotal);
    x_tmp_end(idxSidewardAll) = -x_firCycle(1, idxForwardAll) * sin(thetaTotal) + x_firCycle(1, idxSidewardAll) * cos(thetaTotal);
    
    x_tmp_end = repmat(x_tmp_end,N,1);
    x_trans(:,idxForwardAll) = x_tmp_end(:,idxForwardAll);
    x_trans(:,idxSidewardAll) = x_tmp_end(:,idxSidewardAll);              
    
	% compute current cycle
    x_curCycle(:,idxForwardAll)  = x_trans(:,idxForwardAll)  + x_shift(:, idxForwardAll)*cos(thetaTotal) + x_shift(:, idxSidewardAll)*sin(thetaTotal);
    x_curCycle(:,idxSidewardAll) = x_trans(:,idxSidewardAll) - x_shift(:, idxForwardAll)*sin(thetaTotal) + x_shift(:, idxSidewardAll)*cos(thetaTotal);
    x_curCycle(:,idxPelvisRotation) = x_curCycle(:,idxPelvisRotation) + theta; % pelvis rotation
    
    % concatenate all cycles
    x = [x; x_curCycle];
    
end
% figure(); model.showStick(x(1:200, 1:model.nStates)'); view(90, 90);

% Get range of movement
idxPelvisX = model.extractState('q', 'pelvis_tx'); % forwards
idxPelvisZ = model.extractState('q', 'pelvis_tz'); % sidewards
xrange = [min(x(:,idxPelvisX))-1  , max(x(:,idxPelvisX))+1];
yrange = [-0.2, 2];
zrange = [min(x(:,idxPelvisZ))-1  , max(x(:,idxPelvisZ))+1]; % only for 3D otherwise empty
if isempty(zrange); zrange = []; end % change size of emtpy matrix to avoid warning
range = [xrange; yrange; zrange];

% Initialize figure window
figure();
clf;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.4]);
set(gcf, 'color', 'white');
    
% Plot all frames
nframes = size(x,1);
for i=1:nframes
    
    model.showStick(x(i,1:model.nStates)', range, showStickOpt{:});
    
    if ~makeGIF
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
    else  
        if (i==1)
            filename = [filename, '.gif'];
            gif(filename, 'DelayTime', 1/fps, 'frame',gca); % use gcf for entire figure
        else
            gif;
        end
    end
    
end

hold on
px = x(:,model.extractState('q','pelvis_tx'));
py = x(:,model.extractState('q','pelvis_ty'));
pz = x(:,model.extractState('q','pelvis_tz'));
plot3(pz,px,py);
hold off
if ~makeGIF
    F = getframe(gcf);
    writeVideo(videoWriter,F);
else
    gif;
end

if ~makeGIF
    close(videoWriter);
end
disp(['Animation has been saved on ' filename]);

end
