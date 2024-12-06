%======================================================================
%> @file HelperFunctions/plotVariable.m
%> @brief Function to plot the mean and optional also plus minus SD
%> @details
%> Details: plotVariable()
%>
%> @author Marlies Nitschke
%> @date May, 2018
%======================================================================

%======================================================================
%> @brief Function to plot the mean and optional also plus minus SD
%>
%> @details
%> Figure has to be opened. 'hold on' has to be active.
%>
%> Returns an emtpy vector for the handles which were not plotted.
%>
%> @param  t              Double array: Time or sample points of the data (nSamples x 1)
%> @param  dataMean       Double array: Mean of the  tracking data variable. It is only
%>                        plotted if any entry is not nan. (nSamples x 1)
%> @param  dataVar        (optional) Double array: Variance of the  tracking
%>                        data variable. Can be emtpy to be skipped. (nSamples x 1)
%> @param  lineSpec       (optional) String: Line specification for plotting the mean (If you want to skip it use '')
%> @param  colorSpec      (optional) Double array or String: Color of the area (If you want to skip it use '')
%> @param  faceAlpha      (optional) Double: Value for 'FaceAlpha'. (default: 0.1) (If you want to skip it use [])
%> @param  lineWidth      (optional) Double: Value for 'LineWidth'. (default: 0.5)
%> @retval hMean          Handle of line plot
%> @retval hVar           Handle of line plot
%======================================================================
function [hMean, hVar] = plotVariable(t, dataMean,dataVar,lineSpec, colorSpec, faceAlpha, lineWidth)
	
% check input
if nargin < 4
    lineSpec = '';
end
if nargin < 5
    colorSpec = '';
end
if nargin < 6 || isempty(faceAlpha)
    faceAlpha = 0.1;
end
if nargin < 7
    lineWidth = 0.5;
end

% transform input into a column vector
t = t(:);
dataMean = dataMean(:);
if nargin > 2 && ~isempty(dataVar)
    dataVar = dataVar(:);
end

% plot dataMean +- sd 
if nargin > 2 && ~isempty(dataVar) && any(~isnan(dataVar(:))) && sum(dataVar(~isnan(dataVar))) ~=0
    % input is given && not empty && value which is not non && no zero variance for all entries (would be single trial)
    
    % plot area
    dataSD = sqrt(dataVar);
    x = [t; t(end:-1:1)];
    y = [dataMean(:)-dataSD(:); dataMean(end:-1:1)+dataSD(end:-1:1)];	
    x = x(~isnan(y));
    y = y(~isnan(y));
    hVar = fill(x, y, colorSpec,'FaceAlpha',faceAlpha, 'EdgeColor', 'none');

else 
    hVar = [];
end


% plot the mean
if any(~isnan(dataMean(:))) && ~isempty(lineSpec)
    if ischar(lineSpec)
        hMean = plot(t, dataMean, lineSpec, 'LineWidth', lineWidth);
    else
        hMean = plot(t, dataMean, 'Color', lineSpec, 'LineWidth', lineWidth);
    end
else
    hMean = [];
end


end