%======================================================================
%> @file @TrackingData/extractData.m
%> @brief TrackingData function to extract only specific data from the TrackingData object
%> @details
%> Details: TrackingData::extractData()
%>
%> @author Marlies Nitschke
%> @date June, 2018
%======================================================================

% ======================================================================
%> @brief Function to extract only specific data from the TrackingData object
%>
%> @details
%> All inputs except type can also be empty. If it is empty it is ignored.
%>
%> Example:
%> @code
%> trackingDataExtracted = trackingData.extractData('gyro', '', [1, 0, 0]);
%> @endcode
%>
%> @param   obj         TrackingData class object which should be used to extract data.
%>                      This object will not be manipulated.
%> @param   type        Char: Type which will be extracted
%> @param   names       (optional) Char or cell of chars: Names which will be extracted
%> @param   directions  (optional) Vector or cell of vectors: Directions which will be extracted (Used for IMUs)
%> @param   positions   (optional) Vector or cell of vectors: Positions which will be extracted (Used for IMUs)
%> @param   units       (optional) Char or cell of chars: Units which will be extracted
%>
%> @retval  obj_out     TrackingData class object with the extracted data
% ======================================================================
function obj_out = extractData(obj, type, names, directions, positions, units)


if nargin < 2
   error('TrackingData:extractData', 'You have to define at least a type which should be extracted');
end

% get logical array for type
if ischar(type)
    condType = strcmp(obj.variables.type, type);
end

% get logical array for names
if nargin > 2 && ~isempty(names) && (ischar(names) || iscell(names))
    condNames = ismember(obj.variables.name, names);
else
   condNames = true(size(condType)); % no restriction by the names
end

% get logical array for directions
if nargin > 3 && ~isempty(directions) && (isvector(directions) || iscell(directions))
    if isvector(directions) && ~iscell(directions)
        condDir = ismember(obj.variables.direction, directions, 'rows');
    elseif iscell(directions)
        condDir = false(size(condType)); % init
        for iDir = 1 : length(directions)
            condDir = condDir | ismember(obj.variables.direction, directions{iDir}, 'rows'); % search for each cell entry
        end
    end
    
else
   condDir = true(size(condType)); % no restriction by the directions
end

% get logical array for positions
if nargin > 4 && ~isempty(positions) && (isvector(positions) || iscell(positions))
    if isvector(positions)
        condPos = ismember(obj.variables.position, positions, 'rows');
    elseif iscell(positions)
        condPos = false(size(condType)); % init
        for iPos = 1 : length(positions)
            condPos = condPos | ismember(obj.variables.position, positions{iPos}, 'rows'); % search for each cell entry
        end
    end
    
else
   condPos = true(size(condType)); % no restriction by the positions
end

% get logical array for units
if nargin > 5 && ~isempty(units) && (ischar(units) || iscell(units))
    condUnits = ismember(obj.variables.unit, units);
else
   condUnits = true(size(condType)); % no restriction by the units
end

% extract the rows which fulfill all conditions
obj_out = copy(obj); % deep copy
obj_out.variables = obj.variables(condType & condNames & condDir & condPos & condUnits, :);

% set processing info
st = dbstack;
fctName = st(1).name;
T = TrackingData.createDefaultProcessingTable();
T.type{1} = fctName;
obj_out.performedProcessing = [obj_out.performedProcessing; T];


end