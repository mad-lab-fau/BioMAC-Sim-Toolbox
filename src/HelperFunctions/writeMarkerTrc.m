%======================================================================
%> @file writeMarkerTrc.m
%> @brief Function to write marker data to a .trc file
%>
%> @author Marlies Nitschke
%> @date November, 2021
%======================================================================

%======================================================================
%> @brief Function to write marker data to a .trc file
%>
%> @details
%> This function assumes that the data is given in the correct units. 
%> The unit is set in the trc file according to the optional input.
%> Independently of inMM, units will never be converted.
%>
%> @param  times         Double vector: Time in seconds (nNodes x 1)
%> @param  data          Double matrix: Data, e.g. all generalized coordinates (nNodes x 3*nMarker)
%> @param  names         Cell of chars: Names of the data which will be used as headers (1 x nMarker)
%> @param  filenameAll   String: Filename to save the data including path and file extension ('.trc')
%> @param  inMM          (optional) Boolean: True if the data is given in millimeter (default: 0) 
%======================================================================
function writeMarkerTrc(times, data, names, filenameAll, inMM)

% Extract and check input
nNodes = numel(times);
assert(nNodes == size(data, 1), 'The number of points in ''times'' have to be consistent with the number of rows in ''data''.');
nMarker = numel(names);
assert(nMarker*3 == size(data, 2), 'The must be 3 entries (X, Y, Z) in number of entries in ''data'' for each entry in ''names''.');
dataRate = round(1/mean(diff(times))); % We assume here that we can round the data rate

if nargin > 4 && inMM
    unitsStr = 'mm';
    formatSpecData = '%.5f\t';
else % default
    unitsStr = 'm';
    formatSpecData = '%.8f\t';
end

% Create file
[~, fileName, fileExtension] = fileparts(filenameAll); 
fid = fopen(filenameAll, 'w');
if fid == -1
   error('Connot open %s. Maybe it is open in Opensim', fillnameAll);
end

% Write header
fprintf(fid,'PathFileType\t%i\t(X/Y/Z)\t%s\n',4, [fileName fileExtension]);
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','DataRate', 'CameraRate', 'NumFrames', 'NumMarkers', 'Units', 'OrigDateRate', 'OrigDataStartFrame', 'OrigNumFrames');
fprintf(fid,'%i\t%i\t%i\t%i\t%s\t%i\t%i\t%i\n',dataRate, dataRate, nNodes, nMarker, unitsStr, dataRate, 1, nNodes);
fprintf(fid,'Frame#\tTime\t');
for iMarker = 1 : nMarker
    fprintf(fid,'%s\t\t\t', names{iMarker});
end
fprintf(fid,'\n\t\t');
for iMarker = 1 : nMarker
    fprintf(fid,'X%i\tY%i\tZ%i\t', iMarker, iMarker, iMarker);
end
fprintf(fid,'\n\n');

% Write data
for iNode = 1 : nNodes
    fprintf(fid,'%i\t', iNode);
    fprintf(fid,'%.5f\t', times(iNode));
    fprintf(fid,formatSpecData, data(iNode, :));
    fprintf(fid,'\n');
end

% Close file
fclose(fid);

end