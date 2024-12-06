%======================================================================
%> @file @TrackingData/resampleData.m
%> @brief TrackingData function to resample tracking data
%> @details
%> Details: TrackingData::resampleData()
%>
%> @author Marlies Nitschke
%> @date May, 2018
%======================================================================

% ======================================================================
%> @brief Function to resample tracking data
%>
%> @details
%> - This function uses interp1().
%> - It resamples only data of the type 'angle', 'translation', 'moment','GRF', 'CoP', 'GRM', 
%>   'acc', 'gyro', 'marker', or 'time'. Thus 'speed' and 'duration' are not resampled.
%> 
%> @todo Rethink the adaption of the events. Currently we keep only the
%> start at 1 and the end at nSamples+1. However, the events are not
%> exactly at the index 1 or nSamples+1 
%>
%> @param   obj             TrackingData class object which should be resampled
%> @param   nNodes          Double: Number of nodes
%> @param   plotIt          (optional) Boolean: If true, plot the signal before and after
%>                          resampling for comparison. (default: 0)
%> @param   plotSamples     (optional) Boolean: If true, plot the single sample points explecitly. (default: 0)
%> @param   plotVariance    (optional) Boolean: If true, plot the variance. (default: 0)
% ======================================================================
function resampleData(obj, nNodes, plotIt, plotSamples, plotVariance)

if nargin < 3
   plotIt = 0; 
end
if nargin < 4 
   plotSamples = 0;
end
if nargin < 5
   plotVariance = 0;
end

% make a copy of the input 
variables_in = obj.variables;

% time vectors
nSamples = obj.nSamples;
tIn = (0:(nSamples-1))'/nSamples;
tOut = (0:(nNodes-1))'/(nNodes);
    
% resample data
types = {'angle', 'translation', 'moment', 'GRF', 'CoP', 'GRM', 'acc', 'gyro', 'marker', 'time'};
idxToResample = find(ismember(obj.variables.type, types));
for iVar = idxToResample'
    obj.variables.mean{iVar} = double(interp1(tIn, variables_in.mean{iVar},tOut,'linear','extrap'));
    obj.variables.var{iVar} = double(interp1(tIn, variables_in.var{iVar},tOut,'linear','extrap')); 
end

% Adapt events -> Keep only start and end event
% Start is at 1
iStart = find(obj.movementEvents.index == 1);
if ~isempty(iStart)
    startEvent = obj.movementEvents.name{iStart};
end
% End is at nSamples+1
iEnd = find(obj.movementEvents.index == nSamples +1); % nSamples before resampling
if ~isempty(iEnd)
    endEvent = obj.movementEvents.name{iEnd};
end
% Assign start and end
obj.movementEvents = table({startEvent, endEvent}',[1, nNodes+1]'); % nNodes after resampling
obj.movementEvents.Properties.VariableNames = {'name','index'};

% set processing info
st = dbstack;
fctName = st(1).name;
T = TrackingData.createDefaultProcessingTable();
T.type{1} = fctName;
T.nSamplesIn{1} = nSamples;
T.nSamplesOut{1} = nNodes;
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
    presentTypes = unique(obj.variables.type(idxToResample));
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
            plotVariable(1:nSamples                , variables_in.mean{idxCurrentType(iVar)} , varIn , lineIn);
            plotVariable((1:nNodes)*nSamples/nNodes, obj.variables.mean{idxCurrentType(iVar)}, varOut, lineOut);
            
            % description
            title(obj.variables.name{idxCurrentType(iVar)},'interpreter','none');
            ylabel(obj.variables.unit{idxCurrentType(iVar)});
       end
       legend('In', 'Out');
    end
end

end
