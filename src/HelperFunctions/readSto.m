%======================================================================
%> @file readSto.m
%> @brief Function to read data to a .sto file
%>
%> @author opensim-matlab Toolbox
%======================================================================

%======================================================================
%> @brief Function to read data to a .sto file
%>
%> @details
%> This function assumes that the data is given in the correct units, i.e.
%> it is not performing any conversion.
%>
%> The function is currently in the original status of the opensim-matlab
%> toolbox.
%> 
%> @todo 
%> In future it might be helpful to convert the data struct into a table
%> compatable with simVarTable.
%>
%> @param   filenameAll   String: Filename to save the data including path and file extension ('.sto')
%> @retval  data          Struct: Data where each column in the file is respresented as a struct field
%======================================================================
function [data] = readSto(filenameAll)

tempData = dlmread(filenameAll,'\t',7,0);
[m n] = size(tempData);

% get the header names
fid = fopen(filenameAll,'r');
for i = 1:7
    remain = fgetl(fid);
end
fclose(fid);

for i = 1: n
    [s{i}, remain] = strtok(remain);
end

for i = 1 : n 
    data.(strrep(s{i}, '.', '_')) = tempData(:,i);     
end



