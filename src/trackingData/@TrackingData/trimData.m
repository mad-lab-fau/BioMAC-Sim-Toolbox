%======================================================================
%> @file @TrackingData/trimData.m
%> @brief TrackingData function to trim tracking data
%> @details
%> Details: TrackingData::trimData()
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================

% ======================================================================
%> @brief Function to trim tracking data
%>
%> @details
%> - It trims the data of the type 'angle', 'translation', 'moment','GRF', 'CoP', 'GRM', 
%>   'acc', 'gyro', 'marker', or 'time'. 
%> - Thus 'speed' and 'duration' are only single samples and are currently not trimmed, 
%>   but removed since they are not valid anymore. But there might be another better 
%>   approach.
%> - Only the events between iStart and iEnd+1 are kept. If there no events for 1 and N+1,
%>   "Start" and "End" are added since they are required in other functions.
%>
%> @param   obj             TrackingData class object which should be trimmed
%> @param   iStart          Int: Index of start sample (inclusive)
%> @param   iEnd            Int: Index of end sample (inclusive)
%> @param   plotIt          (optional) Boolean: If true, plot the signal before and after
%>                          resampling for comparison. (default: 0)
%> @param   plotSamples     (optional) Boolean: If true, plot the single sample points explecitly. (default: 0)
%> @param   plotVariance    (optional) Boolean: If true, plot the variance. (default: 0)
% ======================================================================
function trimData(obj, iStart, iEnd, plotIt, plotSamples, plotVariance)

if nargin < 4
   plotIt = 0; 
end
if nargin < 5 
   plotSamples = 0;
end
if nargin < 6
   plotVariance = 0;
end

% make a copy of the input 
nSamplesIn = obj.nSamples;
variables_in = obj.variables;

% resample data
types = {'angle', 'translation', 'moment', 'GRF', 'CoP', 'GRM', 'acc', 'gyro', 'marker', 'time'};
idxToTrim = find(ismember(obj.variables.type, types));
for iVar = idxToTrim'
    obj.variables.mean{iVar} = variables_in.mean{iVar}(iStart:iEnd);
    obj.variables.var{iVar} = variables_in.var{iVar}(iStart:iEnd); 
end

% Adapt events
% -> Keep only events between iStart and iEnd+1
idxIsIn = obj.movementEvents.index >= iStart & obj.movementEvents.index <= iEnd+1;
name = obj.movementEvents.name(idxIsIn);
index = obj.movementEvents.index(idxIsIn)-iStart+1;

% ->  If there no events for 1 and N+1, add "Start" and "End"
if ~any(index == 1)
    name = ['Start'; name];
    index = [1; index];
end
if ~any(index == obj.nSamples+1)
    name = [name; 'End'];
    index = [index; obj.nSamples+1];
end

% ->  Write events to table
obj.movementEvents = table(name,index);
obj.movementEvents.Properties.VariableNames = {'name','index'};

% set processing info
st = dbstack;
fctName = st(1).name;
T = TrackingData.createDefaultProcessingTable();
T.type{1} = fctName;
T.iStart{1} = iStart;
T.iEnd{1} = iEnd;
obj.performedProcessing = [obj.performedProcessing; T];

% plot the incoming and outgoing data
if plotIt
    % set line specifictions
    if plotSamples
        lineIn = '-r*';
        lineOut = '-bs';
    else
        lineIn = '-r';
        lineOut = '-b';
    end
    
    % plot one figure for each relevant data type
    presentTypes = unique(obj.variables.type(idxToTrim));
    for iType = 1 : length(presentTypes)
       figure('name', ['resampleData:' presentTypes{iType}]);
       clf;
       set(gcf,'units','normalized','outerposition',[0 0 1 1]);
       
       idxCurrentType = find(ismember(obj.variables.type, presentTypes{iType}));
       nVar = numel(idxCurrentType);
       nCol = 4;
       nRows = ceil(nVar/nCol);
       for iVar = 1 : nVar
            subplot(nRows, nCol, iVar); hold on;
            
            % get variance
            if plotVariance 
                varIn =  variables_in.var{idxCurrentType(iVar)};
                varOut =  obj.variables.var{idxCurrentType(iVar)};
            else
                varIn = [];
                varOut = [];
            end
            
            % plot it
            plotVariable(1:nSamplesIn, variables_in.mean{idxCurrentType(iVar)} , varIn , lineIn);
            plotVariable(iStart:iEnd, obj.variables.mean{idxCurrentType(iVar)}, varOut, lineOut);
            
            % description
            title(obj.variables.name{idxCurrentType(iVar)},'interpreter','none');
            ylabel(obj.variables.unit{idxCurrentType(iVar)});
       end
       legend('In', 'Out');
    end
end

end
