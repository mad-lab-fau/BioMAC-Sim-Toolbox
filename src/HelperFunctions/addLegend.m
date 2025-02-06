%======================================================================
%> @file addLegend.m
%> @brief Function to add the legend for a plot
%> @details
%> Details: addLegend()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Function to add the legend for a plot
%>
%> @details
%> By default, Matlab does not set a face alpha for the icon of patches in
%> the legend. In this function the face alpha value is set. However, the 
%> value has to be known and has to be given as input in allFaceAlphas.
%>
%> The legend interpreter is set to Latex.
%>
%> Example how to call this function:
%> @code
%> allHandels = [hSim, hTrackAV, hTrackSD];
%> allEntries = {'Simulated', 'Tracked mean', 'Tracked mean $\pm$ SD'};
%> legendFontSize % by default is legendFontSize = 9;
%> allFaceAlphas = [nan, nan, trackFaceAlpha]; % with e.g. standingFaceAlpha = 0.1;
%> lgnd = addLegend(allHandels, allEntries, allFaceAlphas);
%> @endcode
%>
%> @param  allHandels      Handle array: All handles which should get an entry
%> @param  allEntries      Cell array of Strings: All legend entries (size of allHandels)
%> @param  allFaceAlphas   (optional) Double array: Containing the alpha values for the
%>                         each patch and NaNs for other legend entries. If it is not
%>                         given the alpha values in the legend will not be set. (size of allHandels)
%> @param  legendFontSize  (optional) Double: Legend font size (by default,
%>                         legendFontSize = 9)
%> @retval lgnd            matlab.graphics.illustration.Legend: Handle of the legend
%======================================================================
function lgnd = addLegend(allHandels, allEntries, allFaceAlphas, legendFontSize)
       
    
    if isempty(allHandels)
       lgnd = legend();
       return;
    end

    if nargin < 4
        legendFontSize = 9;
    end
   
    % Add legend
    [lgnd, lgndIcons] = legend(allHandels, allEntries, 'Interpreter', 'latex', 'FontSize', legendFontSize);
        
    % Set transparency for patches
    if  nargin> 2 && ~isempty(allFaceAlphas)
        % Get indices of patches
        useFaceAlpha = zeros(size(allHandels));
        for iHand = 1 : length(allHandels)
            if isa(allHandels(iHand), 'matlab.graphics.primitive.Patch')
                useFaceAlpha(iHand) = 1;
            end
        end
        
        
        useFaceAlpha = find(useFaceAlpha);
        
        % Get face alphas of the patches
        faceAlphas = allFaceAlphas(useFaceAlpha);
        
        % Get handles of legend icons of patches
        patchInLegend = findobj(lgndIcons, 'type', 'patch');
        if numel(patchInLegend) ~= numel(faceAlphas)
           error('Number of patches has to be identical to the number of entries in the faceAlphas array.'); 
        end
        
        % Set face alpha
        for iPatch = 1 : numel(patchInLegend)
            set(patchInLegend(iPatch), 'facea', faceAlphas(iPatch));
        end
    end
   
    
end
