%======================================================================
%> @file setLegendBelowSubfigures.m
%> @brief Function to position a legend centered below all subfigures
%> @details
%> Details: setLegendBelowSubfigures()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Function to position a legend centered below all subfigures
%>
%> @details
%> The upper corner will be set 1.5 cm below the bottom of the lowest
%> subfigures.
%>
%> @param  lgnd            matlab.graphics.illustration.Legend: Handle of the legend
%> @param  subFigSettings  Struct: Containing the settings which were used to place the
%>                         subfigures. It has to has the following fields: 
%>                         width, originRight, oneRight, originUp which are assumed to 
%>                         be given in centimeters.
%> @param  nCol            Double: Number of columns in the figure
%======================================================================
function setLegendBelowSubfigures(lgnd, subFigSettings, nCol)
    
    if isempty(lgnd)
        return;
    end

    % Compute effective width of figure
    effectiveFigWidth = 2*subFigSettings.originRight + nCol * subFigSettings.width + (nCol-1) * subFigSettings.oneRight;
    
    % Get legendWidth and legendHeight which was already set
    unitPlot = 'centimeters';
    set(lgnd, 'Unit', unitPlot);
    curPosition = get(lgnd, 'Position');
    legendWidth = curPosition(3);
    legendHeight = curPosition(4);
    
    % Compute position from left origin
    legendLeft = effectiveFigWidth /2 - legendWidth/2;
    
    % Set position from bottom origin
    distLegFig = 1.5;
    legendBottom = subFigSettings.originUp - legendHeight - distLegFig;
    
    % Set position of legend
    position = [legendLeft, legendBottom, legendWidth, legendHeight];
    set(lgnd, 'Position', position)
    
end