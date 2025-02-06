%======================================================================
%> @file writeMotSto.m
%> @brief Function to write data to a .mot or .sto file
%>
%> @author Marlies Nitschke
%> @date July, 2019
%======================================================================

%======================================================================
%> @brief Function to write data to a .mot or .sto file
%>
%> @details
%> This function assumes that the data is given in the correct units. 
%> isDegrees is set in the mot file according to the optional input.
%> Independently of isDegrees, units will never be converted.
%>
%> @param  times         Double vector: Time in seconds (nNodes x 1)
%> @param  data          Double matrix: Data, e.g. all generalized coordinates (nNodes x nSig)
%> @param  names         Cell of chars: Names of the data which will be used as headers (1 x nSig)
%> @param  filenameAll   String: Filename to save the data including path and file extension ('.mot' or '.sto')
%> @param  inDegrees     (optional) Boolean: True if the data is given in degrees (default: 0) 
%======================================================================
function writeMotSto(times, data, names, filenameAll, inDegrees)

% Extract and check input
nNodes = numel(times);
assert(nNodes == size(data, 1), 'The number of points in ''times'' have to be consistent with the number of rows in ''data''.');
nSig = numel(names);
assert(nSig == size(data, 2), 'The number of entries in ''names'' have to be consistent with the columns of rows in ''data''.');

if nargin > 4 && inDegrees
    inDegreesStr = 'yes';
else % default
    inDegreesStr = 'no';
end

% Create file
[~, fileName, fileExtension] = fileparts(filenameAll); 
fid = fopen(filenameAll, 'w');
if fid == -1
   error('Connot open %s. Maybe it is open in Opensim', fillnameAll);
end

% Write header
fprintf(fid,'%s\n',[fileName, fileExtension]);
fprintf(fid,'nRows=%i\n',nNodes);
fprintf(fid,'nColumns=%i\n',nSig+1); % number of signals + 1 (for time)
fprintf(fid,'inDegrees=%s\n',inDegreesStr);
fprintf(fid,'endheader\n');
fprintf(fid,'%s','time');
for iSig = 1 : nSig
    fprintf(fid,'\t%s', names{iSig});
end
fprintf(fid,'\n');

% Write data
for iNode = 1 : nNodes
    fprintf(fid,'%10.5f', times(iNode));
    fprintf(fid,' %10.5f',data(iNode, :));
    fprintf(fid,'\n');
end

% Close file
fclose(fid);

end