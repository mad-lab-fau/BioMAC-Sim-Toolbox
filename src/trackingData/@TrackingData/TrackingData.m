% ======================================================================
%> @file TrackingData.m
%> @brief Matlab class to define tracking data
%>
%> @author Marlies Nitschke
%> @date March, 2018
% ======================================================================

% ======================================================================
%> @brief The class to define tracking data
%>
%> @details
%> This class was implemented to ensure a consistent data structure
%> of tracking data
%>
%> This class contains methods which are generally applicable to load or
%> process tracking data. Project specific methods
%> should be saved in our personal scripts folder!
%>
%> This class inherits from matlab.mixin.Copyable. This means that we can
%> make a deep copy of a TrackingData object:
%> @code
%> trackingData2 = copy(trackingData);
%> @endcode
%>
%> The method plotVarTable(trackingData.variables) can be used for visualization:
%> @code
%> plotVarTable(trackingData.variables);
%> @endcode
%> 
%> @todo The property isSymmetric is very confusing. Because it is not
%> telling something about the symmetry of the data. If we would like to
%> track other movements later we should have both isSymmetric and isPeriodic.
%>
%> @todo Test the implementation with IMU data.
% ======================================================================
classdef TrackingData < handle & matlab.mixin.Copyable

    properties(SetAccess = private)
        %> String: Name of the study
        studyName
        %> Datetime: dd-MMM-yyyy of the study created using:
        %> @code
        %> t = datetime(2018,12,31)
        %> @endcode
        studyDate
        %> String: Name of the participant
        participantName
        %> Double: Height of participant in m
        participantHeight
        %> Double: Mass of participant in kg
        participantMass
        %> Double: Age of participant in years
        participantAge
        %> String: Gender of participant ('female' or 'male')
        participantGender
        %> String: Experiment ID of participant
        participantID
        %> Cell array of strings: Description of the single data trials (nTrials x 1)
        trialList
        %> String: Movement type which was recorded (e.g. walking)
        movementType
        %> String: More detailed description of the movement. Like that
        %> they were walking at prefered speed
        movementDescription
        %> Table: Movement events saved with name and index. This is should
        %> contain at least a description of the start and the end of the
        %> data.
        movementEvents
        %> String: Running style of the participant (i.e. 'rearfoot', 'midfoot', 'forefoot')
        runningStyle
        %> Boolean: 1 if symmetry assumption can be used, i.e. only one half of the movement
        %> (e.g. a half gait cycle) is contained in the tracking data and
        %> we assume that the second have is identically to get a periodic
        %> movement. 0 if we have already a periodic movement which can
        %> directy be used in the periodicity constraint.
        isSymmetric
        %> Table: Data with columns type, name, mean, var (=variance), unit, 
        %> and optionally segment, position, and direction for IMU and marker data.
        %> - Each row must have an unique identifier defined by the columns
        %>   specified in TrackingData.VARIABLEIDENTIFIER
        %> - The following types are supported in this class: 'angle',
        %>   'translation', 'moment', 'GRF', 'CoP', 'GRM', 'acc', 'gyro', 'marker', 'duration', 'speed'
        %> - If only one trial was used, the variance of a variable has to be a zero vector.
        %> - Data can be extracted using for example:
        %>   @code
        %>   trackingData.variables.mean{strcmp(trackingData.variables.type, 'speed')}
        %>   @endcode
        variables
        %> Table: Information about performed processing. The default table
        %> is created in TrackingData.createDefaultProcessingTable where all
        %> the fields are explained. If you are writing a new function 
        %> manipulating the data, please also log this in performedProcessing.
        performedProcessing
        
        %> String: User name of the computer 
        %> (Will be automatically set in the constructor)
        userName
        %> String: Name of the computer
        %> (Will be automatically set in the constructor)
        computerName
        %> Datetime: Creation date with 'dd-MMM-yyyy HH:mm:ss'
        %> (Will be automatically set in the constructor)
        creationTime
        %> String: Git hash of the current commit
        %> (Will be automatically set in the constructor)
        gitHashString
        %> String: Url of git repository
        %> (Will be automatically set in the constructor)
        gitURL        
    end
    
    properties (Dependent, SetAccess = protected)
        %> Double: Number of samples (same for each variable)
        nSamples
        %> Double: Number of variables
        nVariables
        %> Double: Number of trials
        nTrials
    end
    
    properties(Constant)
        %> Cell containing Strings: Column names of TrackingData.variables
        %> which are needed to identify an unique variable row. Position,
        %> direction, and segment are only needed for IMU and marker data.
        VARIABLEIDENTIFIER = {'type', 'name', 'unit', 'position', 'direction', 'segment'};
    end
    
    properties
        %> String: Filename with path, without extension
        filename
        %> String: Personal comment
        comment
    end
    
    methods
        %======================================================================
        %> @brief Constructor setting default TrackingData object
        %>
        %> @param   varargin    (optional) Parameter settings: Parameter
        %>                      name as string followed by the parameter value. Multiple
        %>                      parameters can be set.
        %>                      E.g.: TrackingData('sourceInfo', sourceInfo, 'triallist', triallist);
        %> @retval  obj         TrackingData class object
        %======================================================================
        function [obj] = TrackingData(varargin)
            
            % add general information
            [~, obj.userName] = system('whoami');
            try
                obj.computerName = char(java.net.InetAddress.getLocalHost.getHostName);
            catch
                warning('Could not get computer name.');
                obj.computerName = '';
            end
            obj.creationTime = datetime('now','Format','dd-MMM-yyyy HH:mm:ss');
            try
                [~,gitHashString] = system('git rev-parse HEAD');
            catch
                warning('Could not get git hash');
                gitHashString = '';
            end
            obj.gitHashString = gitHashString;
            try
                [~,gitURL] = system('git config --get remote.origin.url');
            catch
                warning('Could not get git url');
                gitURL = '';
            end
            obj.gitURL = gitURL;

            obj.performedProcessing = table(); % create empty table
            
            % check weather number of input arguments is even
            if mod(nargin, 2) ~= 0
                error('TrackingData:TrackinData', 'The number of inputs has to be even.');
            end
            
            % go through varargin and set values
            for iVar = 1 : 2 : nargin
                obj.setProperty(varargin{iVar}, varargin{iVar+1});
            end
        end
        
        %======================================================================
        %> @brief Function to set a property
        %>
        %> @details
        %> Warning: It is not checked weather the type of the value is
        %> correct!
        %>
        %> @param   obj       TrackingData class object
        %> @param   name      String: Name of the property
        %> @param   value     Not defined: Value of the property
        %======================================================================
        function setProperty(obj, name, value)
            
            if ~ischar(name)
                error('TrackingData:setParameter', ['Second input argument has to be a string and not a ' class(name)]);
            end
            
            try
                obj.(name) = value;
            catch e
                error('TrackingData:setParameter', e.message);
            end
        end
        
        %======================================================================
        %> @brief Function to save the TrackingData object
        %>
        %> @details
        %> The TrackingData object will be saved in a variable called 'data'
        %> and saved using the given filename.
        %>
        %> @param   obj          Result class object
        %> @param   filename     (optional) String: Filename with path, without extension
        %>                       (default: TrackingData.filename)
        %======================================================================
        function save(obj, filename)
            if nargin > 1
                obj.filename = filename;
            elseif isempty(obj.filename)
                error('TrackingData:save', 'The filename property is emtpy.');
            end
            
            data = obj;
            save(obj.filename, 'data');
            
        end
        
        %======================================================================
        %> @brief Function to save the TrackingData object as struct
        %>
        %> @details
        %> The TrackingData object will be saved in a struct called 'dataStruct'
        %> and saved using the given filename.
        %>
        %> @param   obj          TrackingData class object
        %> @param   filename     (optional) String: Filename with path, without extension
        %>                       (default: TrackingData.filename)
        %======================================================================
        function saveStruct(obj, filename)
            if nargin > 1
                obj.filename = filename;
            elseif isempty(obj.filename)
                error('TrackingData:save', 'The filename property is emtpy.');
            end
            
            % get struct
            dataStruct = struct(obj);
            
            % save
            save(obj.filename, 'dataStruct');
            
        end
        
        %======================================================================
        %> @brief Function to get indices of a specific row
        %>
        %> @param   obj          TrackingData class object
        %> @param   subtable     Table: Subtable with one row containing
        %>                       all columns which should be identical
        %> @retval  idxRow       Double array: Indices of the rows which are
        %>                       identically to the entries in the subtable
        %======================================================================
        function idxRow = getVariableIdx(obj, subtable)
            % Check if subtable has only one row
            if height(subtable) ~= 1
                error('TrackingData:getVariableIdx', 'subtable has to have only one row');
            end
            
            % Get column names which are used in the subtable
            variableIdentifierNames = subtable.Properties.VariableNames;
            
            % Go over all rows to find the rows which are identiacally. 
            idxRow = [];
            for iVar = 1 : obj.nVariables
               curRow = obj.variables(iVar, variableIdentifierNames);
               if isequaln(curRow, subtable)
                  idxRow = [idxRow, iVar]; 
               end
            end            
            
        end
        
        
        
        %> @cond DO_NOT_DOCUMENT
        %======================================================================
        %> @brief Function returning TrackingData.nSamples
        %>
        %> @details
        %> This function checks all variables of type 'angle', 'translation', 'moment',
        %> 'GRF', 'CoP', 'GRM', 'acc', 'gyro', or 'marker'.
        %>
        %> @param   obj         TrackingData class object
        %> @retval  nSamples    Double: Number of samples (same for each variable)
        %======================================================================
        function nSamples = get.nSamples(obj)
                        
            if isempty(obj.variables)
                nSamples = 0;
            else
                
                % get number of samples for angle, translation, moment,
                % GRF, CoP, GRM, acc, gyro, and marker variables
                types = {'angle', 'translation', 'moment', 'GRF', 'CoP', 'GRM', 'acc', 'gyro', 'marker'};
                nSampAll = (cellfun(@length,obj.variables.mean));
                nSampAll = nSampAll(ismember(obj.variables.type, types));
                
                % check whether these are equal
                if sum(diff(nSampAll)) ~= 0
                    typesString = sprintf(' %s,',types{:});
                    typesString(end) = '.'; % replace last comma
                    warning('TrackingData:nSamples', ['Number of samples is not equal for the types' typesString]);
                end

                % set number of samples
                if isempty(nSampAll )
                    nSamples = 0;
                else                    
                    nSamples = nSampAll(1);
                end
            end
        end
        
        %======================================================================
        %> @brief Function returning TrackingData.nVariables
        %>
        %> @param   obj         TrackingData class object
        %> @retval  nVariables  Double: Number of variables
        %======================================================================
        function nVariables = get.nVariables(obj)
            if isempty(obj.variables)
                nVariables = 0;
            else
                nVariables = height(obj.variables);
            end
        end
        
        %======================================================================
        %> @brief Function returning TrackingData.nTrials
        %>
        %> @param   obj         TrackingData class object
        %> @retval  nTrials  Double: Number of trials
        %======================================================================
        function nTrials = get.nTrials(obj)
            if isempty(obj.trialList)
                nTrials = 0;
            else
                nTrials = length(obj.trialList);
            end
        end
                
        %> @endcond
    end
    
    
    methods(Static)
        
        %======================================================================
        %> @brief Function to get a TrackingData object from a struct
        %>
        %> @param   dataStruct  Struct: TrackingData content
        %> @retval  obj         TrackingData class object
        %======================================================================
        function obj = struct2object(dataStruct)
            
            % create object
            obj = TrackingData();  
            
            % get dependent or constant properties since this can not be
            % overwritten
            mc = ?TrackingData;
            propNamesAll = {mc.PropertyList.Name};
            isDep = [mc.PropertyList.Dependent];
            isConst = [mc.PropertyList.Constant];
            propNamesDepConst = propNamesAll(isDep | isConst);
            
            % get fields of the struct which should be copied
            fieldNames = setdiff(fieldnames(dataStruct)', propNamesDepConst);
            
            % copy fields
            for fn =  fieldNames
                try
                    setProperty(obj, fn{1}, dataStruct.(fn{1})); % and copy
                catch 
                    warning('Could not copy field %s', fn{1});
                end
            end
            
        end
        
        % ======================================================================
        %> @brief Function to load the data into an object which was previously saved using saveStruct()
        %>
        %> @param   filename  String: Name of the file containing a struct called data
        %> @retval  obj       TrackingData class object
        % ======================================================================
        function obj = loadStruct(filename)
            if ~ischar(filename)
                error('TrackingData:loadStruct', ['First input argument has to be a string and not a ' class(filename)]);
            end
            
            % load file
            dataStruct = load(filename, 'dataStruct');
            if isfield(dataStruct, 'dataStruct')
                dataStruct = dataStruct.dataStruct;
            end
            
            % create object using the struct
            obj = TrackingData.struct2object(dataStruct);
            obj.setProperty('filename', filename);
            
        end
        
        % Declared function to compute the mean and variance of several data sets
        meanTrackingData = takeMeanAndVar(plotIt, varargin);
        
    end
    
    methods(Access = private, Static)
        
        % ======================================================================
        %> @brief Function to load the data into an object which was previously saved using saveStruct()
        %> @private
        %> @static
        %>
        %> @details
        %> Fields:
        %> - type: This has to set in all functions manipulating the data.
        %>   It's the name of the function which manipulated the data.
        %> - processingOfTrial: This is used in TrackingData.takeMeanAndVar.
        %>   It constains a table with the columns trialList and 
        %>   all the columns contained in TrackingData.performedProcessing.
        %> - entriesCorrected: This is used in TrackingData.correctGyroSum
        %>   and in TrackingData.correctVariance. It is a subtable of the
        %>   variable table describing the variables which were corrected.
        %> - nSamplesIn and nSamplesOut: These are used in TrackingData.resampleData.
        %>   They contain the number of samples of the input and of the
        %>   output, respectively.
        %> - iStart and iEnd: These are used in TrackingData.trimData. They contain
        %>   the index of start and end used to trim the data.
        %>
        %> @retval  T   Table: Default table which is used for TrackingData.performedProcessing
        % ======================================================================
        function T = createDefaultProcessingTable()
            variableNames = {'type',  'processingOfTrial', 'entriesCorrected', 'nSamplesIn',  'nSamplesOut', 'iStart', 'iEnd'};
            T = table({''}, {table()}, {table()}, {[]}, {[]}, {[]}, {[]}, 'VariableNames', variableNames);
        end
        
    end
    
end

