%======================================================================
%> @file HelperFunctions/convertUnit.m
%> @brief Function to convert unit of a variables table or simVarTable
%> @details
%> Details: convertUnit()
%>
%> @author Marlies Nitschke
%> @date March, 2022
%======================================================================

%======================================================================
%> @brief Function to convert unit of a variables table or simVarTable
%>
%> @details
%> The function converts the unit of all available data (sim, mean, var,
%> mean_extra, and var_extra). It will also save unitNew into the column
%> 'unit'.
%>
%> Usage for a simVarTable:
%> @code
%> simVarTable = convertUnit(simVarTable, 'translation', 'm', 'mm', 1000);
%> @endcode
%>
%> Usage for a TrackingData object:
%> @code
%> trackingData.setProperty('variables', convertUnit(trackingData.variables, 'angle', 'rad', 'deg', 180/pi));
%> @endcode
%>
%> @param  varTable       Table: Variables table with at least the columns 'type' and 'unit'.
%> @param  type           String: Data type for which the unit should be changed.
%> @param  unitOld        String: Old unit of the data which should be replaced.
%> @param  unitOld        String: New unit of the data which should be replaced.
%> @param  factor         Double: Factor applied to the data for conversion. 
%>                        For variance, we apply factor^2.
%> @param  varTable       Table: Variables table with converted units.
%======================================================================
function varTable = convertUnit(varTable, type, unitOld, unitNew, factor)

% Error checking
variableNames = varTable.Properties.VariableNames;
assert(any(strcmp(variableNames, 'unit')), 'A column with the unit must be given in the simVarTable.');

% Rows to convert
idxConvert = find(ismember(varTable.type, type) & ismember(varTable.unit, unitOld));

% Change unit
varTable.unit(idxConvert) = {unitNew};

% Change simulated data
if any(strcmp(variableNames, 'sim'))
    varTable.sim(idxConvert) = cellfun(@(x) x*factor, varTable.sim(idxConvert), 'UniformOutput', false);
end

% Change mean data
if any(strcmp(variableNames, 'mean'))
    varTable.mean(idxConvert) = cellfun(@(x) x*factor, varTable.mean(idxConvert), 'UniformOutput', false);
    varTable.var(idxConvert) = cellfun(@(x) x*factor^2, varTable.var(idxConvert), 'UniformOutput', false);
end

% Change extra data
if any(strcmp(variableNames, 'mean_extra'))
    varTable.mean_extra(idxConvert) = cellfun(@(x) x*factor, varTable.mean_extra(idxConvert), 'UniformOutput', false);
    varTable.var_extra(idxConvert) = cellfun(@(x) x*factor^2, varTable.var_extra(idxConvert), 'UniformOutput', false);
end

end