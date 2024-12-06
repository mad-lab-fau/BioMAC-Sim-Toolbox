%======================================================================
%> @file tableToLaTeX.m
%> @brief Function to convert a matlab table to a table formated in LaTeX
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Function to convert a matlab table to a table formated in LaTeX
%>
%> @details
%> Warning: This function was not tested for all possible content of
%> tables. Please extend the functionality of this function if you dicover
%> some problems.
%>
%> This function requires the packages "booktabs" and "float".
%>
%> Currently no special formatting of the table is possible.
%>
%> @param table  Table: Table which can contain various number of columns
%>               and rows.
%> @param colFor Cell array of Strings: Format for each column 
%>               (e.g. {'%s','%f'})
%======================================================================
function texString = tableToLaTeX(table, colFor)

if ~isa(table, 'table')
    error('Input has to be an object of class table.');
end

% General infos
colHeaders = table.Properties.VariableNames;
nCols = length(colHeaders);
nRows = height(table);

% Start
texStringSta = [sprintf('\\begin{table} [H] \n'), ...
                sprintf('\\centering \n'), ...
                sprintf('\\begin{tabular}{%s} \n', repmat('l ', 1, nCols)), ...
                sprintf('\\toprule[1pt] \\midrule[0.3pt] \n')];

% Column Headers
if nCols > 1
    texStringHea = [sprintf('{%s} & ', colHeaders{1:end-1}), ...
                    sprintf('{%s} \\\\ \n', colHeaders{end}), ...
                    sprintf('\\midrule \n')];
else
    texStringHea = [sprintf('{%s} \\\\ \n', colHeaders{1}), ...
                    sprintf('\\midrule \n')];
end


% Create format specification for content
if ~iscell(colFor) 
   error('colFor has to be a cell array.');
end
if length(colFor) ~= nCols
    error('colFor has to have %u elements.', nCols);
end

% Create content string containing all rows
texStringCon = '';
for iRow = 1 : nRows
    texStringCurRow = '';
    for iCol = 1 : nCols
        entry = table{iRow, iCol};
        if iscell(entry) && numel(entry) == 1
            entry = entry{:};
        end
        
        texStringCurRowPart = sprintf(colFor{iCol}, entry);
        texStringCurRowPart = strrep(texStringCurRowPart, '_', '\textunderscore ');
        %> @todo Add other special characters or write a function for that
        
        texStringCurRow = [texStringCurRow, texStringCurRowPart];
        if iCol ~= nCols
            texStringCurRow = [texStringCurRow, ' & '];
        end
    end
texStringCon = [texStringCon, ...
                sprintf('%s \\\\ \n', texStringCurRow)];
end


% End
texStringEnd = [sprintf('\\midrule[0.3pt] \\bottomrule[1pt] \n'), ...
               sprintf('\\end{tabular} \n'), ...
               sprintf('\\label{} \n'), ...
               sprintf('\\end{table} \n')];

% Concatenate parts
texString = [texStringSta, texStringHea, texStringCon, texStringEnd];

end