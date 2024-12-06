%======================================================================
%> @file periodicityConstraint.m
%> @brief Collocation function to compute periodicity constraint
%> @details
%> Details: Collocation::periodicityConstraint()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding periodic movement
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%> @param sym           Boolean if movement is symmetric (half period is optimized) or not
%======================================================================
function output = periodicityConstraint(obj,option,X,sym)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'speed') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('State vector X does not contain required states.')
end

if size(X(obj.idx.states),2) ~= (obj.nNodes+1) || size(X(obj.idx.controls),2) ~= (obj.nNodes+1) %> @todo should we add another field to the state vector instead of expanding states to obj.nNodes+1
    error('Model states and controls need to be optimized at collocation nodes + 1')
end
    
%% compute demanded output
nStates = size(obj.idx.states,1);
nControls =  size(obj.idx.controls,1);
icx = 1:nStates; % first indices are for the states
icu = (nStates+1):(nStates+nControls); % next indices are for the controls

% duration
dur = X(obj.idx.dur);
% speed
speed = X(obj.idx.speed);
% forward translation in x direction
unitdisplacementx = zeros(nStates,1);
unitdisplacementx(obj.model.idxForward) = 1;
% sideward translation in z direction
unitdisplacementz = zeros(nStates,1);
unitdisplacementz(obj.model.idxSideward) = 1;

if strcmp(option,'confun') %constraints of periodicity constraint
    output = zeros(nStates+nControls,1);
    
    % compute displacement
    displacementx = unitdisplacementx*dur*speed(1);
    if numel(speed) > 1
        % Sideward displacement in 3D model 
        displacementz = unitdisplacementz*dur*speed(2); % speed in z direction is assumed to be at speed(2)
    else
        % There is no sideward speed given
        displacementz = 0;
    end
    
    if ~sym
        % state must be periodic, with forward displacement of speed*dur
        output(icx) = X(obj.idx.states(:,end)) - X(obj.idx.states(:,1)) - displacementx - displacementz;
        % controls must be periodic
        output(icu) = X(obj.idx.controls(:,end)) - X(obj.idx.controls(:,1));
    else
        % state must be periodic, with mirroring and forward displacement of speed*dur
        output(icx) = X(obj.idx.states(:,end)) - obj.model.idxSymmetry.xsign.*X(obj.idx.states(obj.model.idxSymmetry.xindex,1)) - displacementx - displacementz;
        % controls must be periodic
        output(icu) = X(obj.idx.controls(:,end)) - obj.model.idxSymmetry.usign.*X(obj.idx.controls(obj.model.idxSymmetry.uindex,1));
    end
    
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(nStates+nControls,length(X),obj.Jnnz); %> @todo where to get Jnnz
    

    % compute derivatives of the displacement
    displacementxddur   = unitdisplacementx*speed(1);
    displacementxdspeed = unitdisplacementx*dur;
    if numel(speed) > 1
        displacementzddur   = unitdisplacementz*speed(2);
        displacementzdspeed = unitdisplacementz*dur;
    else
        displacementzddur   = 0;
        displacementzdspeed = [];
    end
    
    if ~sym
        output(icx,obj.idx.states(:,end))   = speye(nStates);
        output(icx,obj.idx.states(:,1))     = -speye(nStates);
        output(icx,obj.idx.dur)             = -displacementxddur-displacementzddur;
        output(icx,obj.idx.speed)           = [-displacementxdspeed, -displacementzdspeed];
        output(icu,obj.idx.controls(:,end)) = speye(nControls);
        output(icu,obj.idx.controls(:,1))   = -speye(nControls);
    else
        output(icx,obj.idx.states(:,end))        = speye(nStates);
        output(icx,obj.model.idxSymmetry.xindex) = -obj.model.idxSymmetry.xsign*ones(1,nStates).*speye(nStates);
        output(icx,obj.idx.dur)                  = -displacementxddur-displacementzddur;
        output(icx,obj.idx.speed)                = [-displacementxdspeed, -displacementzdspeed];
        output(icu,obj.idx.controls(:,end))      = speye(nControls);
        output(icu, obj.idx.controls(obj.model.idxSymmetry.uindex,1)) = -obj.model.idxSymmetry.usign*ones(1,nControls).*speye(nControls);
    end
    
else
    error('Unknown option');
end
end

