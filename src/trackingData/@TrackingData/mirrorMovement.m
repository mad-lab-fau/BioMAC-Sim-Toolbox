% ======================================================================
%> @file @TrackingData/mirrorMovement.m
%> @brief TrackingData function to mirror a tracking data object
%> @details
%> Details: TrackingData::mirrorMovement()
%>
%> @author Marlies Nitschke
%> @date March, 2019
% ======================================================================

% ======================================================================
%> @brief Function to mirror a tracking data object
%>
%> @details
%> This function mirrors left and right sides of a tracking data object.
%> Comparison to a circled shift:
%>  - We do not need to trust movement events.
%>  - We still have a continuous movement (we can't assume a periodic movement).
%>  - We do not have to deal with translations.
%>  - We change also the sideward direction of the movment (e.g. translation and speed).
%>
%> In the implementation we make the following assumptions (feel free to
%> extend the implementation):
%>  - Data of left and right side has to be available.
%>  - Supported variable types: translation, angle, moment, GRF, CoP, GRM, speed,
%>    traveledSpeed, duration
%>  - We use model.idxSymmetry.xsign and model.idxSymmetry.xindex for 
%>    translation, angle and moment.
%>  - For GRF, CoP, and GRM, we assume the ending '_r' and '_l'.
%>  - Movement events for left and right side are denoted with 'R_' and
%>    'L_', respectively (e.g 'L_HS').
%>
%> @param  obj      TrackingData class object which should be mirrored
%> @param  model    Model: Model associated with the tracking data. This is
%>                  required to obtain idxSymmetry
%> @param  plotIt   Boolean: If true, plot all signals and the mean and variance
%> @retval obj_out  TrackingData class object which was mirrored
% ======================================================================
function obj_out = mirrorMovement(obj, model, plotIt)

if nargin < 3
    plotIt = 0;
end

% Make deep copy
obj_out = copy(obj); % deep copy!

% Mirror variables for all supported types
variables = obj_out.variables;
for iVar = 1 : obj.nVariables
   switch variables.type{iVar}
       case {'translation', 'angle', 'moment'}
           idxCurState = model.extractState('q', variables.name{iVar});
           if ~isempty(idxCurState)
               xsign = model.idxSymmetry.xsign(idxCurState);
               xindex = model.idxSymmetry.xindex(idxCurState);
               iVarSym = extractVar(variables, variables.type{iVar}, model.dofs.Properties.RowNames{xindex});
               if isempty(iVarSym)
                  error('Data of left and right side has to be available.');
               end
               variables.mean{iVar} = xsign*obj_out.variables.mean{iVarSym};
           end
           
           % We have to set pelvis_tz_force separatly since it does not use
           % the name of the DOF
           if strcmp(variables.type{iVar}, 'moment')
               idxStates = model.extractState('q');
               idxSideward = model.idxSideward;
               if ~isempty(idxSideward)
                   nameSidewardTranslation = model.dofs.Properties.RowNames{intersect(idxStates, idxSideward)};
                   if strcmp(variables.name{iVar}, [nameSidewardTranslation '_force'])
                       xsign = model.idxSymmetry.xsign(model.extractState('q', nameSidewardTranslation));
                       variables.mean{iVar} = xsign*obj_out.variables.mean{iVar};
                   end
               end    
           end
           
       case {'GRF', 'CoP', 'GRM'}
          strR = '_r';
          strL = '_l';
          xsign = 1;
          if contains(variables.name{iVar}, '_z_')
             xsign = -1; 
          end
          if endsWith(variables.name{iVar}, strR)
            iVarSym = extractVar(variables, variables.type{iVar}, strrep(variables.name{iVar}, strR, strL));
          elseif endsWith(variables.name{iVar}, strL)
            iVarSym = extractVar(variables, variables.type{iVar}, strrep(variables.name{iVar}, strL, strR));
          else
              error('For GRF, CoP, and GRM, we assume the ending ''%s'' and ''%s''.', strR, strL);
          end
          variables.mean{iVar} = xsign*obj_out.variables.mean{iVarSym};              
          
       case 'speed'
           if strcmp(variables.name{iVar}, 'z')
               variables.mean{iVar} = -obj_out.variables.mean{iVar};              
           end
           
       case {'traveledSpeed', 'duration'}
           % Do nothing :)
          
       otherwise
           supportedTypes = {'translation', 'angle', 'moment', 'GRF', 'CoP', 'GRM', 'speed', 'traveledSpeed', 'duration'};
           supportedTypesString = sprintf('%s, ', supportedTypes{:});
           error('Type ''%s'' is not supported by this function. Please extend the function. Only the following types are supported: %s.', ...
               variables.type{iVar}, supportedTypesString(1:end-2));
   end
    
end
obj_out.setProperty('variables', variables);


% Mirror movementEvents
movementEvents = obj_out.movementEvents;
strR = 'R_';
strL = 'L_';
idxR = startsWith(movementEvents.name, strR);
idxL = startsWith(movementEvents.name, strL);
movementEvents.name(idxR) = strrep(movementEvents.name(idxR), strR, strL);
movementEvents.name(idxL) = strrep(movementEvents.name(idxL), strL, strR);
obj_out.setProperty('movementEvents', movementEvents);

% Adapt movement description
obj_out.movementDescription = sprintf('Mirrored version from ''%s''', obj.movementDescription);

% Update performedProcessing table
st = dbstack;
fctName = st(1).name;
T = TrackingData.createDefaultProcessingTable();
T.type{1} = fctName;
obj_out.performedProcessing = [obj_out.performedProcessing; T];

% plot the changes
if plotIt
    types = {'angle', 'translation', 'moment', 'GRF', 'CoP', 'GRM'};
    presentTypes = unique(obj.variables.type(ismember(obj.variables.type, types)));
    for iType = 1 : length(presentTypes)
        figure('name', ['mirrorMovement:' presentTypes{iType}]);
        clf;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        
        idxCurrentType = find(ismember(obj.variables.type, presentTypes{iType}));
        nVar = numel(idxCurrentType);
        nCol = 4;
        nRows = ceil(nVar/nCol);
        for iVar = 1 : nVar
            subplot(nRows, nCol, iVar); hold on;
            
            % plot it
            [~, hIn] = plotVariable(1:obj.nSamples, obj.variables.mean{idxCurrentType(iVar)} , obj.variables.var{idxCurrentType(iVar)}, 'r', [1, 0, 0]);
            [~, hOut] = plotVariable(1:obj.nSamples, obj_out.variables.mean{idxCurrentType(iVar)} , obj_out.variables.var{idxCurrentType(iVar)}, 'g', [0, 1, 0]);
            
            % description
            title(obj.variables.name{idxCurrentType(iVar)},'interpreter','none');
            ylabel(obj.variables.unit{idxCurrentType(iVar)});
        end
        legend([hIn, hOut], {'In', 'out'});
    end
    
end


end

%> @cond DO_NOT_DOCUMENT  
%======================================================================
%> @brief Function to extract index of variable
%> @todo This function should be moved to a generally accesible place. It 
%> could be also merged with Model.extractState() and Model.extractControl().
%======================================================================
function iVar = extractVar(variables, type, name)
if nargin == 3 && ischar(name)
    iVar = find(strcmp(variables.type,type) & strcmp(variables.name,name));
elseif nargin == 3 && iscell(name)
    idxType = strcmp(variables.type,type);
    iVar = zeros(size(name));
    for iName1 = 1 : size(name, 1)
        for iName2 = 1 : size(name, 2)
            iVar(iName1, iName2) = find(idxType & strcmp(variables.name,name{iName1, iName2}));
        end
    end
elseif nargin == 2
    iVar = find(strcmp(variables.type,type));
end
end
%> @endcond

