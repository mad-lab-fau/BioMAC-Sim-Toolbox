%======================================================================
%> @file @TrackingData/useHalfGaitCycleData.m
%> @brief TrackingData cut the tracking data such that only half is used
%>
%> @details
%> Details: TrackingData::useHalfGaitCycleData()
%>
%> @author Anne Koelewijn
%> @date October, 2019
%======================================================================

% ======================================================================
%> @brief Function to cut the tracking data such that only half is used
%>
%> @details
%> Function to use only half the tracking data if the motion is symmetric but the
%> data is of the complete gait cycle.
%>
%> It is just taking the half without taking the movement events into account for cutting.
%>
%> It cuts data of the type 'angle', 'translation', 'moment','GRF', 'CoP', 'GRM',
%> 'acc', 'gyro', 'marker', 'time'. It also halves the 'duration'. The 'speed' is not changed.
%>
%> @param   obj             TrackingData class object which should be divided
% ======================================================================
function useHalfGaitCycleData(obj)    

% half time variables
types = {'angle', 'translation', 'moment', 'GRF', 'CoP', 'GRM', 'acc', 'gyro', 'marker', 'time'};
idxToResample = find(ismember(obj.variables.type, types));
for iVar = idxToResample'
    meandata = obj.variables.mean{iVar};
    vardata = obj.variables.var{iVar};
    if rem(length(meandata),2)~= 0 
        error('Data is of uneven length') 
    end
    obj.variables.mean{iVar} = meandata(1:end/2);
    obj.variables.var{iVar} = vardata(1:end/2); 
end

% half duration
iDur = find(strcmp(obj.variables.type, 'duration'));
if iDur
    obj.variables.mean{iDur} = obj.variables.mean{iDur}/2;
    obj.variables.var{iDur} = (sqrt(obj.variables.var{iDur})/2)^2;
end

% adapt events -> Keep only the onces which are before nSamples+1 (keep end)
iFirstLarger = find(obj.movementEvents.index > obj.nSamples+1, 1);
if ~isempty(iFirstLarger)
    obj.movementEvents = obj.movementEvents(1:iFirstLarger-1,:);
end

% set processing info
st = dbstack;
fctName = st(1).name;
T = TrackingData.createDefaultProcessingTable();
T.type{1} = fctName;
obj.performedProcessing = [obj.performedProcessing; T];

end
