%======================================================================
%> @file calculateMeasures.m
%> @brief Function to caluclate the error and performance measures for a simVarTable
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================

%======================================================================
%> @brief Function to caluclate the error and performance measures for a simVarTable
%>
%> @details
%> It adds the following columns for each measure to the simVarTable
%> - <measure>_sim_track comparing the columns "sim" and "mean"
%> - <measure>_sim_extra comparing the columns "sim" and "mean_extra" 
%>   if there is extra data.
%> - <measure>_track_extra comparing the columns "mean" and "mean_extra"
%>   if there is extra data.
%>
%> Supported measures and their functions:
%> - CMC (calculateCMC)
%> - Pearson (calculateCORR)
%> - rRMSE (calculateRelativeRMSE)
%> - RMSE (calculateRMSE)
%>
%> Example:
%> @code
%> simVarTable = calculateMeasures(simVarTable, {'CMC', 'RMSE'})
%> @endcode
%>
%> @param  smiVarTable  Table: Simulated, tracked and optionally reference data
%>                      created with Collocation.report or Collocation.extractData.
%> @param  measureNames Cell array: Names of the measures which should be computed
%> @retval smiVarTable  Table: Simulated, tracked and optionally reference data
%>                      including columns for the measures
%======================================================================
function simVarTable = calculateMeasures(simVarTable, measureNames)

% Definition of measures
supportedMeasures = {'CMC', 'Pearson', 'rRMSE', 'RMSE'};
supportedFcts = {'calculateCMC', 'calculateCORR', 'calculateRelativeRMSE', 'calculateRMSE'};

% Check if all requested measures are supported
idxNotSupported = find(~ismember(measureNames, supportedMeasures));
assert(isempty(idxNotSupported), 'The following measures are not supported: %s', strjoin(measureNames(idxNotSupported), ', '));
nMeas = length(measureNames);
idxRequested = nan(nMeas, 1);
for iMeas = 1 : nMeas
    idxRequested(iMeas) = find(ismember(supportedMeasures, measureNames{iMeas}));
end

% Check if there is extra data supported for some variable
hasExtraData = ~all(cellfun(@isempty, simVarTable.mean_extra));

% Add columns to the simVarTable
for iMeas = 1 : nMeas
    simVarTable.([measureNames{iMeas} '_sim_track'])(:) = nan;
    if hasExtraData
        simVarTable.([measureNames{iMeas} '_sim_extra'])(:) = nan;
        simVarTable.([measureNames{iMeas} '_track_extra'])(:) = nan;
    end
end

% Go through all rows and compute the measures
nRows = size(simVarTable, 1);
for iRow = 1 : nRows
    
    sim = simVarTable.sim{iRow};
    track =  simVarTable.mean{iRow};
    extra =  simVarTable.mean_extra{iRow};
    
    for iMeas = 1 : nMeas
        if ~isempty(sim) && ~isempty(track)
            simVarTable.([measureNames{iMeas} '_sim_track'])(iRow) = feval(supportedFcts{idxRequested(iMeas)}, sim, track);
        end
        if ~isempty(sim) && ~isempty(extra)
            simVarTable.([measureNames{iMeas} '_sim_extra'])(iRow) = feval(supportedFcts{idxRequested(iMeas)}, sim, extra);
        end
        if ~isempty(track) && ~isempty(extra)
            simVarTable.([measureNames{iMeas} '_track_extra'])(iRow) = feval(supportedFcts{idxRequested(iMeas)}, track, extra);
        end
    end
   
end

end