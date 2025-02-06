%======================================================================
%> @file writeControlXML.m
%> @brief Function to write controls to a .xml file
%>
%> @author Marlies Nitschke
%> @date July, 2019
%======================================================================

%======================================================================
%> @brief Function to write controls to a .xml file
%>
%> @details
%>
%> @param  times         Double vector: Time in seconds (nNodes x 1)
%> @param  controls      Double matrix: Muscle controls (nNodes x nSig)
%> @param  names         Cell of chars: Names of the data which will be used as headers (1 x nSig)
%> @param  filenameAll   String: Filename to save the data including path and file extension ('.xml')
%======================================================================
function writeControlXML(times, controls, names, filenameAll)

% Extract and check input
nNodes = numel(times);
assert(nNodes == size(controls, 1), 'The number of points in ''times'' have to be consistent with the number of rows in ''data''.');
nSig = numel(names);
assert(nSig == size(controls, 2), 'The number of entries in ''names'' have to be consistent with the columns of rows in ''data''.');

% Create file
fid = fopen(filenameAll, 'w');
if fid == -1
    error('Connot open %s. Maybe it is open in Opensim', fillnameAll);
end

% Write data and header to xml file
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" ?>\n');
fprintf(fid,'<OpenSimDocument Version="40000">\n');
fprintf(fid,'	<ControlSet name="Control Set">\n');
fprintf(fid,'		<objects>\n');

for iSig = 1 : nSig % write data for each muscle
    
    fprintf(fid,'			<ControlLinear name="%s">\n', names{iSig});                    % adding muscle's name
    fprintf(fid,'				<is_model_control>true</is_model_control>\n');
    fprintf(fid,'				<extrapolate>true</extrapolate>\n');
    fprintf(fid,'				<default_min>0.02</default_min>\n');
    fprintf(fid,'				<default_max>1</default_max>\n');
    fprintf(fid,'				<filter_on>false</filter_on>\n');
    fprintf(fid,'				<use_steps>false</use_steps>\n');
    fprintf(fid,'				<x_nodes>\n');
    for iNode = 1 : nNodes % write all values (for each time value) for the current muscle
        fprintf(fid,'					<ControlLinearNode>\n');
        fprintf(fid,'						<t>%.10f</t>\n', times(iNode));                 % adding time    value
        fprintf(fid,'						<value>%.10f</value>\n', controls(iNode,iSig)); % adding control value
        fprintf(fid,'					</ControlLinearNode>\n');
    end
    fprintf(fid,'				</x_nodes>\n');
    fprintf(fid,'				<min_nodes />\n');
    fprintf(fid,'				<max_nodes />\n');
    fprintf(fid,'				<kp>100</kp>\n');
    fprintf(fid,'				<kv>20</kv>\n');
    fprintf(fid,'			</ControlLinear>\n');
    
end

fprintf(fid,'		</objects>\n');
fprintf(fid,'		<groups />\n');
fprintf(fid,'	</ControlSet>\n');
fprintf(fid,'</OpenSimDocument>\n');

% Close file
fclose(fid);

end