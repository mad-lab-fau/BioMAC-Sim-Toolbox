%======================================================================
%> @file writeIMUcsv.m
%> @brief Function to write IMU data to a .csv file
%>
%> @author Marlies Nitschke
%> @date July, 2023
%======================================================================

%======================================================================
%> @brief Function to write IMU data to a .csv file
%>
%> @details
%> This functions follows the format definition of APDM.
%>
%> This function assumes that the data is given in the correct units. 
%> The unit for acceleration signals is m/s^2 and for gyroscope signals is
%> rad/s. Units will never be converted.
%>
%> @param  times              Double vector: Time in seconds (nNodes x 1)
%> @param  variableTableIMU   Table: Variables table with at least the columns type, name, direction, and mean.
%> @param  filenameAll        String: Filename to save the data including path and file extension ('.csv')
%======================================================================
function writeIMUcsv(times, variableTableIMU, filenameAll)

% Extract and check input
nNodes = numel(times);
variableTableIMU = variableTableIMU(ismember(variableTableIMU.type, {'acc', 'gyro'}), :);
nSampAll = (cellfun(@length,variableTableIMU.mean));
assert(sum(diff(nSampAll)) == 0, 'The number of samples must be equal for all rows');
assert(nNodes == nSampAll(1), 'The number of points in ''times'' have to be consistent with the samples in ''variableTableIMU''.');
dataRate = round(1/mean(diff(times))); % We assume here that we can round the data rate

% Create cell that is written into csv
dataCell = cell(4+nNodes, 1+height(variableTableIMU)); % 4: header rows; 1: time column
dataCell{1, 1} = 'Test Name:';
[~, fileName, ~] = fileparts(filenameAll); 
dataCell{1, 2} = fileName;
dataCell{2, 1} = 'Sampling Rate:';
dataCell{2, 2} = dataRate;
dataCell{2, 3} = 'Hz';
dataCell{3, 1} = 'Time';
dataCell{4, 1} = 's';
dataCell(5:end, 1) = num2cell(times);

% Fill cell with IMU data
variableTableIMU = sortrows(variableTableIMU, {'name', 'type'}); % Group by sensor
for iVar = 1 : height(variableTableIMU)
    switch variableTableIMU.type{iVar}
        case 'acc'
            name = 'Acceleration';
            unit = 'm/s^2';
        case 'gyro'
            name = 'Angular Velocity';
            unit = 'rad/s';
    end
    if isequal(variableTableIMU.direction(iVar, :), [1, 0, 0])
        axis = 'X';
    elseif isequal(variableTableIMU.direction(iVar, :), [0, 1, 0])
        axis = 'Y';
    elseif isequal(variableTableIMU.direction(iVar, :), [0, 0, 1])
        axis = 'Z';
    else 
        error('Direction used in row %i is not supported', iVar)
    end
    dataCell{3, 1+iVar} = sprintf('%s/%s/%s', variableTableIMU.name{iVar}, name, axis);
    dataCell{4, 1+iVar} = unit;
    dataCell(5:end, 1+iVar) = num2cell(variableTableIMU.mean{iVar});
end

% Save cell as csv
writecell(dataCell,filenameAll);

end