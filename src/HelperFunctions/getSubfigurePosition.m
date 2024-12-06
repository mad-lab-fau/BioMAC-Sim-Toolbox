%======================================================================
%> @file getSubfigurePosition.m
%> @brief Function to get the position of a subfigure
%> @details
%> Details: getSubfigurePosition()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Function to get the position of a subfigure
%>
%> @details
%> If you want to have a defined spacing between subfigures (e.g. to export
%> it into tikz for a publication), this function can compute the positions
%> of the subfigures for you. 
%>
%> This function can be used like the following to set the position of a
%> subfigure:
%> @code
%> set(hSub, 'Unit', 'centimeters', 'Position', getSubfigurePosition(iAng, nRow, nCol, subFigSettings));
%> @endcode
%> During this call, you have to set the unit which was also used in
%> subFigSettings.
%>
%> @param   iFig            Double: Index of the subfigure
%> @param   nRow            Double: Number of rows in the figure
%> @param   nCol            Double: Number of columns in the figure
%> @param   subFigSettings  Struct: Containing the settings to place the
%>                          subfigures in the fields: width, height, originRight, 
%>                          originUp, oneRight, and oneUp.
%> @retval  pos             Double vector: [left, bottom, width, height] 
%======================================================================
function pos = getSubfigurePosition(iFig, nRow, nCol, subFigSettings)
    % Get how often we have to go up and right
    [iCol, iRow] = ind2sub([nCol, nRow], iFig); % order is different. => row and col has to be switched
    nUp    = nRow - iRow;
    nRight = iCol - 1;
    
    % Compute the position
    left   = subFigSettings.originRight + nRight * (subFigSettings.oneRight + subFigSettings.width );
    bottom = subFigSettings.originUp    + nUp    * (subFigSettings.oneUp    + subFigSettings.height);
    pos = [left, bottom, subFigSettings.width, subFigSettings.height];
end
