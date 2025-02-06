%======================================================================
%> @file compareElements.m
%> @brief Function to find differences between two elements
%>
%> @author Marlies Nitschke
%> @date March, 2018
%======================================================================

%======================================================================
%> @brief Function to find differences between data elements
%>
%> @details
%> This function will give you a fast overview what are the differences between
%> two elements. For example, this is usefull to compare two simulation
%> results. 
%> @code
%> diffTabOut = compareElements(result1, result2)
%> @endcode
%>
%> It was implemented recursively. Hence, it is calling itself until it
%> reaches for example a numeric value or a string. The tables containing
%> the differences will be concatenated.
%>
%> Warning:
%> - The function was not fully tested.
%> - The function output may be not suitable for your usecase.
%> - The comparison of some classes could be missing. In this case, an
%>   error will be thrown. => Feel free to extend this function.
%>
%> @param  element1      Element of arbitrary class
%> @param  element2      Element of arbitrary class
%> @param  diffTabIn     (optimal)Table: Differences which were already disconvered
%> @param  elemPath      (optimal)String: Path which was used to reach the element (e.g. test.name)
%> @retval diffTabIn     Table: Differences which were disconvered
%======================================================================
function diffTabOut = compareElements(element1, element2, diffTabIn, elemPath)

if nargin < 3
    % Create an empty table
    variableNames = {'Path', 'element1', 'element2', 'description'};
    diffTabOut = cell2table(cell(0,length(variableNames)), 'VariableNames', variableNames);
else
    % There was already some difference discovered before
    diffTabOut = diffTabIn;
end

if nargin < 4
    % We are at the root of the elements
    elemPath = '';
end

% We only have to do something if they are not equal
if ~isequaln(element1, element2)
    
    % check whether the elements are of the same class
    if ~strcmp(class(element1), class(element2))
        diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, 'Different classes'}];
        return;
    end
    
    % do different test depending on the class
    if istable(element1) 
        % check whether they contain the same columns
        varNames1 = element1.Properties.VariableNames;
        varNames2 = element2.Properties.VariableNames;
        varNamesDiff1 = setdiff(varNames1, varNames2);
        varNamesDiff2 = setdiff(varNames2, varNames1);
        varNamesDiff = {varNamesDiff1{:}, varNamesDiff2{:}};
        if ~isempty(varNamesDiff)
           % there are missing columns
           des = strcat('These columns are not contained in both structs: ', sprintf('  %s',varNamesDiff{:}));
           diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];
        end
        varNamesBoth = intersect(varNames1, varNames2);
        
        % check whether they have the same number of rows
        if height(element1) ~= height(element2)
            des = 'The number of rows is different';
            diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];
        else
            % check whether they contain the same rows
            rows1 = element1.Properties.RowNames;
            rows2 = element2.Properties.RowNames;
            rowsDiff1 = setdiff(rows1, rows2);
            rowsDiff2 = setdiff(rows2, rows1);
            rowsDiff = {rowsDiff1{:}, rowsDiff2{:}};
            if ~isempty(rowsDiff)
                % there are missing rows
                des = strcat('These rows are not contained in both structs: ', sprintf('  %s',rowsDiff{:}));
                diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];
            end
        
            % go through all rows and columns recursivly
            rowsBoth = intersect(rows1, rows2);
            if ~isempty(rowsBoth) % This works only if a row name is given
                for iRow = 1 : length(rowsBoth)
                    for iVar = 1 : length(varNamesBoth)
                        diffTabOut = compareElements(element1{rowsBoth{iRow}, varNamesBoth{iVar}}, element2{rowsBoth{iRow}, varNamesBoth{iVar}} , ...
                            diffTabOut, [elemPath '{' rowsBoth{iRow} ', ' varNamesBoth{iVar} '}']);
                    end
                end
            else
                for iRow = 1 : height(element1)
                    for iVar = 1 : length(varNamesBoth)
                        % We assume here that they have the same order
                        diffTabOut = compareElements(element1{iRow, varNamesBoth{iVar}}, element2{iRow, varNamesBoth{iVar}} , ...
                            diffTabOut, [elemPath '{' num2str(iRow) ', ' varNamesBoth{iVar} '}']);
                    end
                end
            end
        end
        
    elseif isobject(element1) 
        % go through all properties recursivly
        prop = properties(element1);
        for iProp = 1 : length(prop)
            diffTabOut = compareElements(element1.(prop{iProp}), element2.(prop{iProp}) , diffTabOut, [elemPath '.' prop{iProp}]);
        end        
    elseif isstruct(element1)
        % check whether they contain the same fields
        fields1 = fields(element1);
        fields2 = fields(element2);
        fieldsDiff1 = setdiff(fields1, fields2);
        fieldsDiff2 = setdiff(fields2, fields1);
        fieldsDiff = {fieldsDiff1{:}, fieldsDiff2{:}};
        if ~isempty(fieldsDiff)
           % there are missing fields
           des = strcat('These fields are not contained in both structs: ', sprintf('  %s',fieldsDiff{:}));
           diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];
        end
        fieldsBoth = intersect(fields1, fields2);
        
        % check whether they are struct arrays and contain the same number
        % of rows
        if length(element1) ~= length(element2)
           % there are different number of rows
           des = strcat('There are different number of rows in both struct arays: ', sprintf('  %s',fieldsDiff{:}));
           diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];
        end
        
        % go through all fields recursivly
        if length(element1) == 1 && length(element2) == 1
            for iFiel = 1 : length(fieldsBoth)
                diffTabOut = compareElements(element1.(fieldsBoth{iFiel}), element2.(fieldsBoth{iFiel}) , diffTabOut, [elemPath '.' fieldsBoth{iFiel}]);
            end  
        else
            % Process it as table to use the code for tables
            diffTabOut = compareElements(struct2table(element1), struct2table(element2) , diffTabOut, ['struct2table(' elemPath ')']);
        end 

    elseif isnumeric(element1)
        % check whether the size is equal
        if sum(size(element1) ~= size(element2)) ~= 0
           diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, 'Not of same size.'}];
           return;
        end
        
        % find indices where the elements are different
        idxDif= find(element1 ~= element2);
        des = strcat('Different numeric values for these linear indices: ', sprintf('  %s',num2str(idxDif')));
        diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];        
        
    elseif iscell(element1)
        if sum(size(element1) ~= size(element2)) ~= 0
            diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, 'Not of same size.'}];
            return;
        end
        
        % go through all entries recursivly
        for iEnt = 1 : numel(element1)
            diffTabOut = compareElements(element1{iEnt}, element2{iEnt}, diffTabOut, [elemPath '{' num2str(iEnt) '}']);
        end  
        
    elseif ischar(element1)
        des = 'Different string';
        diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];        
    elseif islogical(element1)
        des = 'Different logical value';
        diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}];  
    elseif isa(element1,'function_handle')
        des = 'Different function handles';
        diffTabOut = [diffTabOut; {elemPath, {element1}, {element2}, des}]; 
    else 
        error(['No implementation for class ' class(element1)]);
    end


end