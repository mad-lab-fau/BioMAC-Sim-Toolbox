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

% Determine if we are using the new dynamicPeriodicConstraints (where optimization nodes = N)
is_combined_periodic = isfield(obj.constraintTerms, 'name') && any(strcmp({obj.constraintTerms.name}, 'dynamicPeriodicConstraints'));
if is_combined_periodic
    % Find sym from dynamicPeriodicConstraints
    iCon = find(strcmp({obj.constraintTerms.name}, 'dynamicPeriodicConstraints'), 1);
    sym = obj.constraintTerms(iCon).varargin{1};
    
    if isfield(obj.idx,'speed')
        speed = norm(X(obj.idx.speed));
    else
        speed = 0; 
    end
    
    % Reconstruct states and controls with the virtual node N+1 appended
    states_orig = X(obj.idx.states);
    controls_orig = X(obj.idx.controls);
    
    x1_first = states_orig(:, 1);
    u1_first = controls_orig(:, 1);
    
    unitdisplacementx = zeros(size(x1_first));
    unitdisplacementx(obj.model.idxForward) = 1;
    unitdisplacementz = zeros(size(x1_first));
    unitdisplacementz(obj.model.idxSideward) = 1;
    
    displacementx = unitdisplacementx * T * X(obj.idx.speed(1));
    if numel(X(obj.idx.speed)) > 1
        displacementz = unitdisplacementz * T * X(obj.idx.speed(2));
    else
        displacementz = 0;
    end
    
    if sym
        x_virtual = obj.model.idxSymmetry.xsign .* x1_first(obj.model.idxSymmetry.xindex) + displacementx + displacementz;
        u_virtual = obj.model.idxSymmetry.usign .* u1_first(obj.model.idxSymmetry.uindex);
    else
        x_virtual = x1_first + displacementx + displacementz;
        u_virtual = u1_first;
    end
    
    states = [states_orig, x_virtual];
    controls = [controls_orig, u_virtual];
    nNodesDur = nNodes + 1;
    h = T / nNodes;
else
    states = X(obj.idx.states);
    controls = X(obj.idx.controls);
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
end

% Initialize parameters
statesd = ( states(:, 2:nNodesDur) - states(:, 1:(nNodesDur-1)) ) / h;
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
    for iNode = 1 : nNodesDur-1
        % Get states and controls
        % (See Collocation.dynamicConstraints as reference)
        curStates = states(:, iNode);

        % Get energy for the current node
        [Edot(:, iNode)] = obj.model.getMetabolicRate_pernode(curStates, statesd(:, iNode), controls(:,iNode), t_stim(:,iNode), name, getCont, epsilon);
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
    
    if strcmp(string(exponent), "log")
        metRate = log(metRate);
        metCost = log(metCost);
        CoT = log(CoT);
        metCostPerMus = log(metCostPerMus);
    else
        metRate = metRate.^exponent;
        metCost = metCost.^exponent;
        CoT = CoT.^exponent;
        metCostPerMus = metCostPerMus.^exponent;
    end
    
%% Calculate gradients only if they are needed:
else
    dmetRatedX = zeros(size(X)); 
    for iNode = 1 : nNodesDur-1
        curStates = states(:, iNode);
        [Edot(:, iNode), dEdotdx, dEdotdu, dEdotdxdot] = obj.model.getMetabolicRate_pernode(curStates, statesd(:, iNode), controls(:,iNode), t_stim(:,iNode), name, getCont, epsilon);
        
        dmetRatedX(obj.idx.states(:, iNode))   = dmetRatedX(obj.idx.states(:, iNode))   + dEdotdx;
        dmetRatedX(obj.idx.controls(:, iNode)) = dmetRatedX(obj.idx.controls(:, iNode)) + dEdotdu;
        dmetRatedX(obj.idx.states(:, iNode))   = dmetRatedX(obj.idx.states(:, iNode))   - dEdotdxdot/h;
        
        if is_combined_periodic && iNode == nNodesDur-1
            % Last node N+1 is virtual (reconstructed from node 1)
            dfdxdot_term = dEdotdxdot/h;
            if sym
                dfdxdot_mapped = zeros(size(dfdxdot_term));
                dfdxdot_mapped(obj.model.idxSymmetry.xindex) = dfdxdot_term .* obj.model.idxSymmetry.xsign;
            else
                dfdxdot_mapped = dfdxdot_term;
            end
            dmetRatedX(obj.idx.states(:, 1)) = dmetRatedX(obj.idx.states(:, 1)) + dfdxdot_mapped;
            
            % Gradient w.r.t duration and speed via displacement offset in virtual node
            ddisp_dT = unitdisplacementx * X(obj.idx.speed(1));
            if numel(X(obj.idx.speed)) > 1
                ddisp_dT = ddisp_dT + unitdisplacementz * X(obj.idx.speed(2));
            end
            dmetRatedX(obj.idx.dur) = dmetRatedX(obj.idx.dur) + sum((dEdotdxdot/h) .* ddisp_dT);
            
            ddisp_dspeed1 = unitdisplacementx * T;
            dmetRatedX(obj.idx.speed(1)) = dmetRatedX(obj.idx.speed(1)) + sum((dEdotdxdot/h) .* ddisp_dspeed1);
            if numel(X(obj.idx.speed)) > 1
                ddisp_dspeed2 = unitdisplacementz * T;
                dmetRatedX(obj.idx.speed(2)) = dmetRatedX(obj.idx.speed(2)) + sum((dEdotdxdot/h) .* ddisp_dspeed2);
            end
        else
            dmetRatedX(obj.idx.states(:, iNode+1)) = dmetRatedX(obj.idx.states(:, iNode+1)) + dEdotdxdot/h;
        end
        
        dmetRatedX(obj.idx.dur)                = dmetRatedX(obj.idx.dur)                - sum(dEdotdxdot.*statesd(:, iNode))/h/(nNodesDur-1);
    end
    
    % Base variables calculations
    metRate = sum(sum(Edot, 2) / (nNodesDur-1)) / bodymass;
    dmetRatedX = dmetRatedX / ((nNodesDur-1) * bodymass);
    
    metCost    = (metRate + 1) / speed;
    metCostPerMus = sum(Edot, 2) / (nNodesDur-1) / bodymass / speed;
    CoT    = metCost / gravity;
    
    % Calculate base dmetCostdX (Quotient rule: d/dx (U/V) = (U'V - UV') / V^2)
    dmetCostdX = dmetRatedX / speed; 
    
    if isfield(obj.idx,'speed')
        % speed = norm(X(obj.idx.speed))
        % d(speed)/d(speed_idx) = X(obj.idx.speed) / speed
        dmetCostdX(obj.idx.speed) = dmetCostdX(obj.idx.speed) - ((metRate + 1) / (speed^2)) * (X(obj.idx.speed) / speed);
    else
        % WARNING: If speed is calculated dynamically, you ideally need the chain rule 
        % for d(speed)/dX here to avoid discontinuous gradients. 
        % For now, updating the duration dependency part of the quotient rule:
        dmetCostdX(obj.idx.dur) = dmetCostdX(obj.idx.dur) - ((metRate + 1) / (speed^2)) * (-speed / T); 
    end
    
    % Final step: Apply Chain Rule for Exponents uniformly across functions and gradients
    if strcmp(string(exponent), "log")
        dmetRatedX = dmetRatedX / metRate;
        dmetCostdX = dmetCostdX / metCost;
        dCoTdX = dmetCostdX;
        metRate = log(metRate);
        metCost = log(metCost);
        CoT = log(CoT);
        metCostPerMus = log(metCostPerMus);
    else
        dmetRatedX = exponent * (metRate .^ (exponent-1)) .* dmetRatedX;
        dmetCostdX = exponent * (metCost .^ (exponent-1)) .* dmetCostdX;
        dCoTdX     = dmetCostdX / gravity^exponent;
        metRate = metRate.^exponent;
        metCost = metCost.^exponent;
        CoT = CoT.^exponent;
        metCostPerMus = metCostPerMus.^exponent;
    end 
end

end



