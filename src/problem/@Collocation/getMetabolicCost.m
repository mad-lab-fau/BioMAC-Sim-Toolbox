%======================================================================
%> @file @Collocation/getMetabolicCost.m
%> @brief Collocation function to calculate metabolic cost of a movement
%> @details
%> Details: Collocation::getMetabolicCost()
%>
%> @author Anne Koelewijn, Marlies Nitschke
%> @date March 23, 2015
%======================================================================

%======================================================================
%> @brief Function to calculate metabolic cost of a movement
%> @details
%> Function to calculate the metabolic cost of a movement using
%> Umberger's model. Ross Miller's code was also used as a reference.
%>
%> See also the metabolic cost paper:
%> A Koelewijn, E Dorschky, A van den Bogert;
%> A metabolic energy expenditure model with a continuous first derivative 
%> and its application to predictive simulations of gait. Computer Methods 
%> in Biomechanics and Biomedical Engineering 21(4):1-11, 2018.
%>
%> This function does only return the energy of the present movement.
%> => We did not multiplied by two for symmetric movements here. 
%> Since metCost is computed from mean energy expenditure, it would be false
%> to multiply the metCost by 2 for symmetric movements!!!
%>
%> To get the cost of a result, you can call:
%> @code
%> [metCost, ~, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);
%> @endcode
%>
%> @todo Agree for one method to define (compute or extract) speed!
%>
%> @param  obj            Collocation object
%> @param  X              Double matrix: State vector (i.e. result) of the problem
%> @param  name           (optional) String: name of the model to be used
%>                        default is umberger, other options are minetti, 
%>                        margaria, houdijk, bhargava, uchida, lichtwark,
%>                        kim
%> @param  getCont        (optional) Boolean: If true, get the continuous version 
%>                        which is needed if we use the output for simulation. (default: 0)
%> @param  epsilon        (optional) Double: Amount of nonlinearity in continuous model
%>                        it should be specified if a continuous model is
%>                        used
%> @param  exponent       (optional) integer: could be used to calculate square, cube
%>                        or nth power of metabolic cost. Only used for derivatives. (default: 1)


%> @retval metCost        Double: Metabolic Cost of movement in J/m/kg
%> @retval metRate        Double: Metabolic Rate of movement in W/kg
%> @retval CoT            Double: Cost of transport of movement in 1
%> @retval metCostPerMus  Double vector: Metabolic Cost of movement in J/m/kg (obj.model.nMus x 1)
%> @retval dmetCostdX     Double vector: Derivative of metCost w.r.t to X (size of X)
%> @retval dmetRatedX     Double vector: Derivative of metRate w.r.t to X (size of X)
%> @retval dCoTdX         Double vector: Derivative of metRate w.r.t to X (size of X)
%======================================================================
function [metCost, dmetCostdX, metCostPerMus, metRate, CoT, dmetRatedX, dCoTdX] = getMetabolicCost(obj, X, name, getCont, epsilon, exponent) % 

% Check whether we should return the continuous version which is needed if we use the output for simulation 
if nargin < 3
    name = 'umberger';
end
if nargin < 4
   getCont = 0; 
   epsilon = 1;
end

if getCont == 1 && nargin == 4 && ~strcmp(name,'minetti')
    error('epsilon should be specified')
end

if nargin < 6
    exponent = 1;
end

% Error checking
if ~isfield(obj.idx,'states') % check whether model states are stored in X
    error('Model states are not stored in state vector X.')
end
if ~isfield(obj.idx,'controls') % check whether controls states are stored in X
    error('Model controls are not stored in state vector X.')
end
if ~isfield(obj.idx,'dur') % check whether duration is stored in X
    error('Duration is not stored in state vector X.')
end

% Extract variables which are needed 
bodymass = obj.model.bodymass;          % Bodymass in kg
gravity  = norm(obj.model.gravity);     % Norm of gravity in m/(s^2)
nMus     = obj.model.nMus;              % Number of muscles
nNodes   = obj.nNodes;                  % Number of nodes of the colocation problem
nNodesDur= obj.nNodesDur;               % Number of nodes defining the duration
T        = X(obj.idx.dur);              % Duration of movement
h        = T/(nNodesDur-1);             % Duration of time step
if isfield(obj.idx,'speed')
    speed    = norm(X(obj.idx.speed));            % Speed in forward direction in m/s
else
    deltaX = sum(abs(diff(X(obj.idx.states(obj.model.extractState('q', 'pelvis_tx'), :)))));
    if isa(obj.model, 'Gait3d')
        deltaZ = sum(abs(diff(X(obj.idx.states(obj.model.extractState('q', 'pelvis_tz'), :)))));
    else
        deltaZ = 0;
    end
    speed = (deltaX + deltaZ) / T;
end

% Initialize paremeters
Edot = zeros(nMus, nNodesDur-1);

%Find stimulation time for the entire gait cycle
if epsilon == 1 %Only do this when postprocessing
    t_stim = obj.getStimTime(X);
else
    t_stim = zeros(obj.model.nMus,nNodes);
end


%% Compute Metabolic Cost without a gradient if nargout doesn't need it to be
if ismember(nargout,[0 1 3 4 5])


    % Compute the energy for all nodes
    states = X(obj.idx.states);
    statesd = ( states(:, 2:nNodesDur) - states(:, 1:(nNodesDur-1)) ) / h;
    controls = X(obj.idx.controls);
    for iNode = 1 : nNodesDur-1
        % Get states and controls
        % (See Collocation.dynamicConstraints as reference)
        curStates = states(:, iNode);

        % Get energy for the current node
        [Edot(:, iNode)] = obj.model.getMetabolicRate_pernode(curStates, statesd(:, iNode), controls(:,iNode), t_stim(:,iNode), name, getCont, epsilon, exponent);
    end

    % Metabolic rate in Watts for each muscle (sum up over all nodes)
    metRate = sum(Edot, 2)  / (nNodesDur-1);

    % Metabolic rate in W/kg for all muscles in total
    metRate = sum(metRate)  / bodymass;

    % Metabolic cost
    metCost    = (metRate+1)    / speed;


    % Cost of transport
    CoT    = metCost    / gravity;

    % Metabolic cost per muscle in J/kg/m
    metCostPerMus = sum(Edot, 2) / (nNodesDur-1) / bodymass / speed;
    dmetCostdX = 0;

%% Calculate gradients only if they are needed:
else
    % Compute the energy for all nodes
    dmetRatedX = zeros(size(X)); 
    states = X(obj.idx.states);
    statesd = ( states(:, 2:nNodesDur) - states(:, 1:(nNodesDur-1)) ) / h;
    controls = X(obj.idx.controls);
    for iNode = 1 : nNodesDur-1
        % Get states and controls
        % (See Collocation.dynamicConstraints as reference)
        curStates = states(:, iNode);

        % Get energy for the current node
        [Edot(:, iNode), dEdotdx, dEdotdu, dEdotdxdot] = obj.model.getMetabolicRate_pernode(curStates, statesd(:, iNode), controls(:,iNode), t_stim(:,iNode), name, getCont, epsilon, exponent);

        % Build dmetRatedX out of dmetRatedx, dmetRatedu, dmetRatedT
        dmetRatedX(obj.idx.states(:, iNode))   = dmetRatedX(obj.idx.states(:, iNode))   + dEdotdx;
        dmetRatedX(obj.idx.controls(:, iNode)) = dmetRatedX(obj.idx.controls(:, iNode)) + dEdotdu;

        %Add dEdotdxdot
        dmetRatedX(obj.idx.states(:, iNode))   = dmetRatedX(obj.idx.states(:, iNode))   - dEdotdxdot/h;
        dmetRatedX(obj.idx.states(:, iNode+1)) = dmetRatedX(obj.idx.states(:, iNode+1)) + dEdotdxdot/h;
        dmetRatedX(obj.idx.dur)                = dmetRatedX(obj.idx.dur)                - sum(dEdotdxdot.*statesd(:, iNode))/h/(nNodesDur-1);
    end

    % Metabolic rate in Watts for each muscle (sum up over all nodes)
    metRate = sum(Edot, 2)  / (nNodesDur-1);
    %To do: this should be summed over all muscles too
    dmetRatedX = dmetRatedX / (nNodesDur-1);

    % Metabolic rate in W/kg for all muscles in total
    metRate = sum(metRate)  / bodymass;
    dmetRatedX = dmetRatedX /  bodymass;

    % Metabolic cost
    metCost    = metRate    / speed;
    dmetCostdX = dmetRatedX / speed;
    if isfield(obj.idx,'speed')
        dmetCostdX(obj.idx.speed) = metRate*(-1/speed^2);
    end

    % Cost of transport
    CoT    = metCost    / gravity;
    dCoTdX = dmetCostdX / gravity;

    % Metabolic cost per muscle in J/kg/m
    metCostPerMus = sum(Edot, 2) / (nNodesDur-1) / bodymass / speed;

end


end



