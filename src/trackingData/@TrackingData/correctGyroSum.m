%======================================================================
%> @file @TrackingData/correctGyroSum.m
%> @brief TrackingData function to set the sum of a periodic gyro signal to 0
%> @details
%> Details: TrackingData::correctGyroSum()
%>
%> @author Marlies Nitschke
%> @date May, 2018
%======================================================================

% ======================================================================
%> @brief Function to set the sum of a periodic gyro signal to 0
%>
%> @param   obj     TrackingData class object which should be corrected
%> @param   plotIt  (optional) Boolean: If true, plot the signal before and after
%>                  correction for comparison. (default: 0)
% ======================================================================
function correctGyroSum(obj, plotIt)

if nargin < 2
    plotIt = 0;
end

% is signal periodic?
isPeriodic = ~obj.isSymmetric;
if ~isPeriodic
    return;  % nothing to correct
end

% make a copy of the input
variables_in = obj.variables;

% substract to get sum(gyro) = 0 for periodic movement
idxGyr = find(ismember(obj.variables.type, 'gyro'));
nSamples = obj.nSamples;
idxVarChanged = false(obj.nVariables, 1);
for iVar = idxGyr'
    gyroSum = sum(obj.variables.mean{iVar});
    if gyroSum ~= 0
        obj.variables.mean{iVar} = obj.variables.mean{iVar} - gyroSum /nSamples;
        idxVarChanged(iVar) = 1;
        warning('TrackingData:correctGyroSum', ['The gyroscope signal with the name ''' obj.variables.name{iVar} ...
            ''' was corrected to ensure that the sum is equal to 0.']);
    end
end
entriesCorrected = obj.variables(idxVarChanged, intersect(obj.variables.Properties.VariableNames, obj.VARIABLEIDENTIFIER));

% set processing info
st = dbstack;
fctName = st(1).name;
T = TrackingData.createDefaultProcessingTable();
T.type{1} = fctName;
T.entriesCorrected{1} = entriesCorrected;
obj.performedProcessing = [obj.performedProcessing; T];

% plot the changes of the gyro signal
if plotIt
    presentTypes = unique(obj.variables.type(idxGyr));
    for iType = 1 : length(presentTypes)
        figure('name', ['correctGyroSum:' presentTypes{iType}]);
        clf;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        
        idxCurrentType = find(ismember(obj.variables.type, presentTypes{iType}));
        nVar = numel(idxCurrentType);
        nCol = 4;
        nRows = ceil(nVar/nCol);
        for iVar = 1 : nVar
            subplot(nRows, nCol, iVar); hold on;
            
            % plot it
            plotVariable(1:obj.nSamples, variables_in.mean{idxCurrentType(iVar)} , [], '-r');
            plotVariable(1:obj.nSamples, obj.variables.mean{idxCurrentType(iVar)} , [], '-b');
            
            % description
            title(obj.variables.name{idxCurrentType(iVar)},'interpreter','none');
            ylabel(obj.variables.unit{idxCurrentType(iVar)});
        end
        legend('In', 'out');
    end
    
end

end