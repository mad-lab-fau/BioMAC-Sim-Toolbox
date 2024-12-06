%======================================================================
%> @file plotStancePhase.m
%> @brief Function to plot a shaded area to indicate the stance phase
%> @details
%> Details: plotStancePhase()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================


%======================================================================
%> @brief Function to plot a shaded area to indicate the stance phase
%>
%> @details
%> Figure has to be opened. 'hold on' has to be active.
%>
%> It is using 0.1 for FaceAlpha to plot the area.
%>
%>
%> @param  t          Double array: time axis corresponding to the boolean standing array
%> @param  standing   Double array: 0: no standing; 1: standing (length of t)
%> @param  minY       Double: min Y value of the area
%> @param  maxY       Double: max Y value of the area
%> @param  minX       Double: min X value to paint (first sample of time signal)
%> @param  maxX       Double: max X value to paint (last sample of time signal)
%> @param  colorSpec  (optional) Double array or String: Color of the area. 
%>                    Can be empty to be skipped. (default: 'k')
%> @param  faceAlpha  (optional) Double: Value for 'FaceAlpha'. (default: 0.1)
%> @retval hStanding  Handle of the last plotted area
%======================================================================
function hStanding = plotStancePhase(t, standing, minY, maxY, minX, maxX, colorSpec, faceAlpha)

if nargin < 7 || isempty(colorSpec)
    colorSpec = 'k';
end

if nargin < 8 || isempty(faceAlpha)
    faceAlpha = 0.1;
end

hStanding = [];

% get the start and stops
starts = find(diff(standing)> 0.5)+1; % = HS
stops = find(diff(standing) < -0.5); % = TO

% sort the indices in ascendig order
indices = [starts; stops];
flags = [zeros(size(starts)); ones(size(stops))]; % 0: starts (HS); 1: stops (TO)

[indices_sort, idx] = sort(indices);
flags_sort = flags(idx);

% plot it
yDimensions = [minY,maxY,maxY,minY];
i = 1;
while  i <= length(indices_sort)
    
    if i == 1 && flags_sort(1) == 1 % starts as first event
        xDimensions = [minX, minX, t(indices_sort(1)), t(indices_sort(1))];
        i = i +1;
    elseif i == length(indices_sort) && flags_sort(end) == 0 % HS as last event
        xDimensions = [t(indices_sort(end)), t(indices_sort(end)), maxX, maxX];
        i = i +1;
    elseif flags_sort(i) == 0 && flags_sort(i+1) == 1 % HS to TO
        xDimensions = [t(indices_sort(i)), t(indices_sort(i)), t(indices_sort(i+1)), t(indices_sort(i+1))];
        i = i +2; % skip the TO
    else % wrong detection
        i = i +1;
        continue; % Continue in loop without plotting
    end
    
    % Plot area
    hStanding = fill(xDimensions, yDimensions, colorSpec,'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
    
end

end