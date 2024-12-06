%======================================================================
%> @file @TrackingData/correctVariance.m
%> @brief TrackingData function to prevent small variances of the variables
%> @details
%> Details: TrackingData::correctVariance()
%>
%> @author Marlies Nitschke
%> @date May, 2018
%======================================================================

% ======================================================================
%> @brief Function to prevent small variances of the variables
%>
%> @details
%> - If nTrial = 1:
%>   Sets the variance to 0 if this is not already the case.
%> - If nTrial > 1: Sets all variances below 10% of the mean variance to 
%>   10% of the mean variance. If the mean variance is smaller than 10% of 
%>   the minimum variance, than minimum variance is used. The minimum variances 
%>   were chosen empirically:
%>      - angle:        var = 4 deg²        (sd = 2 deg)
%>      - translation:  var = 0.0001 m²     (sd = 0.01 m)
%>      - GRF:          var = 0.0004 BW²    (sd = 0.02 BW)
%>      - speed:        var = 0.0001 (m/s)² (sd = 0.01 m/s)
%>      - duration:     var = 0.0001 s²     (sd = 0.01 s)
%>      - acc:          var = 1 (m/s²)²     (sd = 1 m/s²)
%>      - gyro:         var = 0.1 (rad/s)²  (sd = 0.3162 rad/s)
%>
%> @todo If nTrial > 1, we set the minimum variance of speed and duration to
%> 0.0001. This was an arbitrarily chosen value. We should think about
%> this.
%>
%>
%> @param   obj     TrackingData class object which should be corrected
%> @param   plotIt  (optional) Boolean: If true, plot the signal before and after
%>                  correction of the variance for comparison. (default: 0)
% ======================================================================
function correctVariance(obj, plotIt)

if nargin < 2
    plotIt = 0;
end

% make a copy of the input
variables_in = obj.variables;

% check each variable
idxVarChanged = false(obj.nVariables, 1);
for iVar = 1 : obj.nVariables
    
    if obj.nTrials == 1
        if ~all(size(obj.variables.mean{iVar}) == size(obj.variables.var{iVar})) && ... % not same size as mean
                ~all(obj.variables.var{iVar} == 0)                                           % not all entries are 0
            % Set variance to zero (It should already be zero ideally.)
            obj.variables.var{iVar} = zeros(size(obj.variables.mean{iVar}));
        end
    else
        
        % Get minimal variance
        switch obj.variables.type{iVar}
            case 'angle'
                minVar = (2*pi/180)^2;
                unit = 'rad';
            case 'translation'
                minVar = 0.001;
                unit = 'm';
            case 'GRF'
                minVar = 0.0004;
                unit = 'BW';
            case 'CoP'
                minVar = 0.001;
                unit = 'm';
            case 'marker'
                minVar = 0.001;
                unit = 'm';
            case 'speed'
                minVar = 0.0001;
                unit = 'm/s';
            case 'duration'
                minVar = 0.0001;
                unit = 's';
            case 'acc'
                minVar = 1;
                unit = 'm/s2';
            case 'gyro'
                minVar = 0.1;
                unit = 'rad/s';
            otherwise
                warning('TrackingData:correctVariance', ['No minimal variance implemented for type ''' obj.variables.type{iVar} ...
                    '''. It will continue without ensuring a minimum variance. This can cause problems during the optimization.']);
                continue;
        end
        
        % Set minimal variance if the current variance is
        % smaller
        % average the variance over the entire gait cycle, use 10% of that
        % as the min variance
        percentage = 0.1; %10% of mean variance
        if strcmp(obj.variables.type{iVar},'CoP')
            percentage = 0;
        end
        meanVar = mean(obj.variables.var{iVar});
        if meanVar < percentage*minVar
            warning('TrackingData:correctVariance', ['The mean variance of the variable of the type ''' ...
                obj.variables.type{iVar} ''' and the name ''' obj.variables.name{iVar} ...
                ''' was too small (',num2str(meanVar),' ',unit,'). A minimum variance of ',num2str(minVar),' ',unit,' is used instead.']);
            meanVar = minVar;
       end
        idxSmaller = find(obj.variables.var{iVar} < percentage*meanVar);
        if ~isempty(idxSmaller)
            obj.variables.var{iVar}(idxSmaller) = percentage*meanVar;
            idxVarChanged(iVar) = 1;
            warning('TrackingData:correctVariance', ['The variance of the variable of the type ''' ...
                obj.variables.type{iVar} ''' and the name ''' obj.variables.name{iVar} ...
                ''' was too small at the indices ' sprintf(' %d,', idxSmaller) ...
                ' and was set to ' num2str(percentage*meanVar) ' ' unit '.']);
            
        end
        
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

% plot the changes
if plotIt
    types = {'angle', 'translation', 'moment', 'GRF', 'acc', 'gyro', 'CoP', 'marker'};
    presentTypes = unique(obj.variables.type(ismember(obj.variables.type, types)));
    for iType = 1 : length(presentTypes)
        figure('name', ['correctVariance:' presentTypes{iType}]);
        clf;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        
        idxCurrentType = find(ismember(obj.variables.type, presentTypes{iType}));
        nVar = numel(idxCurrentType);
        nCol = 4;
        nRows = ceil(nVar/nCol);
        for iVar = 1 : nVar
            subplot(nRows, nCol, iVar); hold on;
            
            % plot it
            [~, hIn] = plotVariable(1:obj.nSamples, variables_in.mean{idxCurrentType(iVar)} , variables_in.var{idxCurrentType(iVar)}, 'k', [1, 0, 0]);
            [~, hOut] = plotVariable(1:obj.nSamples, obj.variables.mean{idxCurrentType(iVar)} , obj.variables.var{idxCurrentType(iVar)}, 'k', [0, 0, 1]);
            
            % description
            title(obj.variables.name{idxCurrentType(iVar)},'interpreter','none');
            ylabel(obj.variables.unit{idxCurrentType(iVar)});
        end
        legend([hIn, hOut], {'In', 'out'});
    end
    
end


end