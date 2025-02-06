%======================================================================
%> @file readMot.m
%> @brief Function to read data to a .mot file
%>
%> @author opensim-matlab Toolbox
%======================================================================

%======================================================================
%> @brief Function to read data to a .mot file
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
%> @param   filenameAll   String: Filename to save the data including path and file extension ('.mot')
%> @retval  data          Struct: Data where each column in the file is respresented as a struct field
%======================================================================
function data = readMot(filenameAll)


fid = fopen(filenameAll);
c = textscan(fid, '%s','delimiter', '\t');
fclose(fid);


record = 0;
colHeaders = [];

for i = 1:length(c{1,1})
    
    if strcmp(c{1,1}{i},'time')
        record = 1;
    elseif ~isempty( str2num(c{1,1}{i}) )
    break
    end
    
    
    if record == 1
        colHeaders =  [colHeaders { c{1,1}{i} }];
    end
end


temp.colHeaders = colHeaders;
temp.data    = dlmread(filenameAll, '\t', 11, 0);

for i = 1:length(temp.colHeaders)
    data.(temp.colHeaders{i}) = temp.data(:,i);
end