%======================================================================
%> @file @TrackingData/takeMeanAndVar.m
%> @brief TrackingData function to compute the mean and variance of several data sets
%> @details
%> Details: TrackingData::takeMeanAndVar()
%>
%> @author Marlies Nitschke
%> @date May, 2018
%======================================================================

% ======================================================================
%> @brief Function to compute the mean and variance of several data sets
%> @static
%>
%> @details
%> It will set the properties using the following rules:
%> - TrackingData.trialList: Combine all trial lists
%> - TrackingData.movementEvents: Use the first and last event
%>
%> We compute the mean and the variance using the data saved in the column
%> 'mean' without considering the column 'var'.
%>
%> Warning: 
%> - After calling TrackingData.takeMeanAndVar, you have to call again the
%>   preprocessing functions like TrackingData.correctGyroSum or 
%>   TrackingData.correctVariance before starting the simulation.
%> - The mean of severel means is not equal to the mean of all the single
%>   trials.
%> 
%> Examplary function call using filenames:
%> @code
%> datafiles = {'file1.mat', 'file2.mat', 'file3.mat'};
%> meanTrackingData = TrackingData.takeMeanAndVar(1, datafiles{:});
%> @endcode
%> 
%> Examplary function call using TrackingData objects:
%> @code
%> meanTrackingData = TrackingData.takeMeanAndVar(1, trackingData1, trackingData2);
%> @endcode
%> 
%> @todo We have to ensure if we can simply take the mean of the
%> translations. If we do this, we are assuming that the cycle is always at
%> the same position in the global reference frame.
%>
%>
%> @param   plotIt      Boolean: If true, plot all signals and the mean and variance
%> @param   varargin    String: Filenames of files containing TrackingData struct \n
%>                      ... or ... \n
%>                      TrackingData: TrackingData objects \n
%>                      Requirements: 
%>                      - varargin must contain at least two entries to
%>                        compute the mean and the variance out of them.
%>                      - Same studyName, subjectName, movementType, isSymmetric
%>                      - Start and end with the same movmentEvents if there are events listed
%>                      - Contain the same variables (Consitent type, name and unit)
%> @retval  meanTrackingData    TrackingData: Containing the mean and
%>                              variance of the input. 
% ======================================================================
function meanTrackingData = takeMeanAndVar(plotIt, varargin)

nData = length(varargin);

%% check if there are at least two entries in varargin
if nData < 2
    error('TrackingData:trackMeanAndVar', 'You have to input at least two TrackingData objects.');
end

%% check if all entries in varargin were of the class TrackingData
for iData = 1 : nData
    if ~isa(varargin{iData}, 'TrackingData')
        if ischar(varargin{iData})
            try
                varargin{iData} = TrackingData.loadStruct(varargin{iData});
            catch
                error('TrackingData:trackMeanAndVar', ...
                    ['Was not able to load ' num2str(iData) '. entry with the name ''' varargin{iData} '''. ' ...
                    'The file must contain a struct which can be loaded using TrackingData.loadStruct().']);
            end
        else
            error('TrackingData:trackMeanAndVar', ...
                ['The ' num2str(iData) '. entry is neither a TrackingData object nor a string containing the filename. ' ...
                'Instead it is of type ' class(varargin{iData}) '.']);
        end
    end
end

%% create a new object for the mean data
meanTrackingData = TrackingData();

%% check if studyName is the same
studyNames = cell(nData, 1);
for iData = 1 : nData
    % get studynames
    studyNames{iData} = varargin{iData}.studyName;
end
% check if there is a difference (Empty labels are not taken into account)
if sum(diff(cell2mat(studyNames))) > 0
    error('TrackingData:trackMeanAndVar', 'All objects in vargarin have to be from the same study.');
end
% Copy studyName and studyDate
meanTrackingData.setProperty('studyName', varargin{1}.studyName);
meanTrackingData.setProperty('studyDate', varargin{1}.studyDate);

%% check if subjectName is the same
subjectNames = cell(nData, 1);
for iData = 1 : nData
    % get subjectName
    subjectNames{iData} = varargin{iData}.subjectName;
end
% check if there is a difference (Empty labels are not taken into account)
if sum(diff(cell2mat(subjectNames))) > 0
    error('TrackingData:trackMeanAndVar', 'All objects in vargarin have to be from the same subject.');
end
% Copy subjectName, subjectHeight and subjectMass
meanTrackingData.setProperty('subjectName', varargin{1}.subjectName);
meanTrackingData.setProperty('subjectHeight', varargin{1}.subjectHeight);
meanTrackingData.setProperty('subjectMass', varargin{1}.subjectMass);
%% check if movementType is the same
movementTypes = cell(nData, 1);
for iData = 1 : nData
    % get movementType
    movementTypes{iData} = varargin{iData}.movementType;
end
% check if there is a difference (Empty labels are not taken into account)
if sum(diff(cell2mat(movementTypes))) > 0
    error('TrackingData:trackMeanAndVar', 'All objects in vargarin have to be from the same movement.');
end
% Copy movementType and movementDescription
meanTrackingData.setProperty('movementType', varargin{1}.movementType);
meanTrackingData.setProperty('movementDescription', varargin{1}.movementDescription);
meanTrackingData.setProperty('runningStyle', varargin{1}.runningStyle);

%% check if isSymmetric is the same
isSymmetrics = cell(nData, 1);
for iData = 1 : nData
    % get isSymmetric
    isSymmetrics{iData} = varargin{iData}.isSymmetric;
end
% check if there is a difference
if sum(diff(cell2mat(isSymmetrics))) > 0
    error('TrackingData:trackMeanAndVar', 'All objects in vargarin have to be symmetric or not.');
end
% Copy isSymmetric
meanTrackingData.setProperty('isSymmetric', varargin{1}.isSymmetric);

%% check if files start and end with the same events
startEvents = cell(nData, 1);
endEvents = cell(nData, 1);
for iData = 1 : nData

    % get all start and end events
    if ~isempty(varargin{iData}.movementEvents)
        startEventCur = varargin{iData}.movementEvents{varargin{iData}.movementEvents.index == 1, 'name'};
        startEvents{iData} = startEventCur{1};
        endEventCur = varargin{iData}.movementEvents{varargin{iData}.movementEvents.index == varargin{iData}.nSamples+1, 'name'};
        endEvents{iData} = endEventCur{1}; 
    end
   
end
% check if there is a difference (Empty labels are not taken into account)
if sum(sum(diff(cell2mat(startEvents)))) > 0 || sum(sum(diff(cell2mat(endEvents)))) > 0
    error('TrackingData:trackMeanAndVar', 'All objects in vargarin have to start and end with the same movement event.');
end

%% check if all entries contain the same variables
% get column names to identify the variables
variableIdentifierNames = intersect(varargin{1}.VARIABLEIDENTIFIER, varargin{1}.variables.Properties.VariableNames);
% get subtable of this columns
variableIdentifierEntries = varargin{1}.variables(:, variableIdentifierNames);
for iData = 2 : nData
    
    % get column names to identify the variables
    variableIdentifierNamesCur = intersect(varargin{iData}.VARIABLEIDENTIFIER, varargin{iData}.variables.Properties.VariableNames);

    % get subtable of this columns
    variableIdentifierEntriesCur = varargin{iData}.variables(:, variableIdentifierNamesCur);
    
    % check if they contain only the same variables (order can be
    % different!)
    if ~isequaln(variableIdentifierEntries, variableIdentifierEntriesCur)
        identifierStr = sprintf('%s, ', variableIdentifierNamesCur{:});
        msg = 'All objects in vargarin have to contain the same variables. %s have to be consistent.';
        error('TrackingData:trackMeanAndVar', msg, identifierStr(1:end-2));
    end
end

%% Resample data
% find the smallest data length
nSamplesAll = zeros(nData, 1);
for iData = 1 : nData
    nSamplesAll(iData) = varargin{iData}.nSamples;    
end
nSamplesMin = min(nSamplesAll);


% resample the data to nSamplesMin
disp(['Data will be resampled to ', num2str(nSamplesMin), ' samples ----']);
for iData = 1 : nData
    if varargin{iData}.nSamples ~= nSamplesMin
        varargin{iData}.resampleData(nSamplesMin);
    end
    disp(['Resampled data file number ', num2str(iData), ' of ', num2str(nData) ' files.']);
end
disp('Done');

%% Adapt events -> Keep only start and end event
% Start is at 1
iStart = find(varargin{1}.movementEvents.index == 1);
if ~isempty(iStart)
    startEvent = varargin{1}.movementEvents.name{iStart};
end
% End is at nSamples+1
iEnd = find(varargin{1}.movementEvents.index == varargin{1}.nSamples +1); % nSamples before resampling
if ~isempty(iEnd)
    endEvent = varargin{1}.movementEvents.name{iEnd};
end
% Assign start and end
meanTrackingData.movementEvents = table({startEvent, endEvent}',[1, varargin{1}.nSamples+1]'); % nNodes after resampling
meanTrackingData.movementEvents.Properties.VariableNames = {'name','index'};


%% Adapt trialList
meanTrackingData.trialList = cell(nData, 1);
for iData = 1 : nData
    meanTrackingData.trialList{iData} = varargin{iData}.trialList;
end

%% Compute the Mean
% Copy variables table
variables_all = varargin{1}.variables(:, variableIdentifierNames);
% Add a column data 
variables_all.data = cell(varargin{1}.nVariables, 1); 
for iVar = 1 : varargin{1}.nVariables
    iVarCur = varargin{1}.getVariableIdx(variables_all(iVar, variableIdentifierNames));
    variables_all.data{iVar} = NaN(nData, length(varargin{1}.variables.mean{iVarCur}));
end

% Collect all the data
for iData = 1 : nData
    for iVar = 1 : varargin{1}.nVariables
        iVarCur = varargin{iData}.getVariableIdx(variables_all(iVar, variableIdentifierNames));
        variables_all.data{iVar}(iData, :) = varargin{iData}.variables.mean{iVarCur};
    end
end

% Compute mean and variance
meanTrackingData.variables = varargin{1}.variables;
for iVar = 1 : varargin{1}.nVariables
    meanTrackingData.variables.mean{iVar} = mean(variables_all.data{iVar})';
    meanTrackingData.variables.var{iVar}  = var(variables_all.data{iVar})';
end

% set processing info
st = dbstack;
fctName = st(1).name;
T = TrackingData.createDefaultProcessingTable();
T.type{1} = fctName;
% create table with processing infos of all the data
processingOfTrial = TrackingData.createDefaultProcessingTable();
variableNames = varargin{1}.performedProcessing.Properties.VariableNames; % get variable names of performedProcessing
processingOfTrial.trialList = cell(1, 0); % add new column
warning('off','MATLAB:table:RowsAddedExistingVars'); % I was not able to find a implementation avoiding this warning.
iEntry = 0;
for iData = 1 : nData
    for iProc = 1 : height(varargin{iData}.performedProcessing)
        iEntry = iEntry +1;
        for iVarName = 1 : length(variableNames)
            curVarName = variableNames{iVarName};
            processingOfTrial{iEntry, curVarName} = varargin{iData}.performedProcessing{iProc, curVarName};
        end
        processingOfTrial.trialList(iEntry) = varargin{iData}.trialList;
    end
end
warning('on','MATLAB:table:RowsAddedExistingVars'); % Switch warning on again
% assign it to T
T.processingOfTrial{1} = processingOfTrial;
meanTrackingData.performedProcessing = [meanTrackingData.performedProcessing; T];

%% Plot the result
if plotIt
    
    nSamples = meanTrackingData.nSamples;
    types = {'angle', 'translation', 'moment', 'GRF', 'CoP', 'GRM', 'acc', 'gyro', 'marker'};
    idxToPlot= ismember(meanTrackingData.variables.type, types);
    % plot one figure for each relevant data type
    presentTypes = unique(meanTrackingData.variables.type(idxToPlot));
    for iType = 1 : length(presentTypes)
       figure('name', ['takeMeanAndVar:' presentTypes{iType}]);
       clf;
       set(gcf,'units','normalized','outerposition',[0 0 1 1]);
       
       idxCurrentType = find(ismember(meanTrackingData.variables.type, presentTypes{iType}));
       nVar = numel(idxCurrentType);
       nCol = 4;
       nRows = ceil(nVar/nCol);
       for iVar = 1 : nVar
            subplot(nRows, nCol, iVar); hold on;
                    
            % plot all trials
            for iData = 1 : nData
                k1 = plotVariable(1:nSamples, variables_all.data{idxCurrentType(iVar)}(iData, :), [], 'k');
            end
            
            % plot mean with variance
            [k2, k3] = plotVariable(1:nSamples, meanTrackingData.variables.mean{idxCurrentType(iVar)}, meanTrackingData.variables.var{idxCurrentType(iVar)} , 'r', [0, 0, 1]);
            
            % description
            xlim([1, nSamples]);
            title(meanTrackingData.variables.name{idxCurrentType(iVar)},'interpreter','none');
            ylabel(meanTrackingData.variables.unit{idxCurrentType(iVar)});
       end
       if ~isa(k3, 'matlab.graphics.primitive.Patch') && ~isnan(k3)
           legend([k1, k2, k3], {'Trials', 'Mean', 'Mean+-SD'});
       else
           legend([k1, k2], {'Trials', 'Mean'});
       end
    end
    
end

end