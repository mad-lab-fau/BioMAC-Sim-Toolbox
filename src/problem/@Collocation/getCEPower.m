%======================================================================
%> @file @Collocation/getCEPower.m
%> @brief Collocation function to compute the CE Power over the whole movemenent
%> @details
%> Details: Collocation::getCEPower()
%>
%> @author Antonie J. (Ton) van den Bogert, Marlies Nitschke
%> @date August, 2018
%======================================================================

%======================================================================
%> @brief Function to compute the CE Power over the whole movemenent
%>
%> @details
%> Computes the CE power according to Ton's paper:
%>
%> A. J. van den Bogert, M. Hupperets, H. Schlarb, and B. Krabbe. Predictive musculoskeletal
%> simulation using optimal control: effects of added limb mass on energy cost and kinematics
%> of walking and running. Proceedings of the Institution of Mechanical Engineers, Part P:
%> Journal of Sports Engineering and Technology, 226(2):123{133, 2012.
%>
%> His old implementation was used as basis for this implementation.
%>
%> We do only take the simulated motion into account, i.e. we do not
%> multiply the output times two for symmetric simulations.
%>
%> @todo Implement powerd. There were multiple implementations of it. See
%> comments. Add a description for the outcome powerd and rename it!
%>
%>
%> @param  obj            Collocation object
%> @param  X              Double matrix: State vector (i.e. result) of the problem
%> @retval totalCEPower   Double: Total mechanical work of the muscle fibres in W
%> @retval powers         Double array: CE power output in W of each muscle and time point
%>                        (Model.nMus x Collocation.nNodes)
%======================================================================
function [totalCEPower, powers] = getCEPower(obj, X)

    % Check if the state vector X contains all required information
    if ~isfield(obj.idx,'dur')  % check whether duration are stored in X
        error('Movement duration is not stored in state vector X.')
    end
    if ~isfield(obj.idx,'states') % check whether states are stored in X
        error('Model states are not stored in state vector X.')
    end

	% Extract information
    nNodesDur = obj.nNodesDur;
    nMus = obj.model.nMus;
    duration = X(obj.idx.dur);
    h = duration/(nNodesDur-1);
    states = X(obj.idx.states);
    statesd = ( states(:, 2:nNodesDur) - states(:, 1:(nNodesDur-1)) ) / h;
    
    % Calculate power between each pair of nodes
    powers = zeros(nMus, nNodesDur-1);
    for iNode = 1 : nNodesDur-1

        % Call mex function of the model
        powers(:, iNode) = obj.model.getMuscleCEpower(states(:,iNode), statesd(:, iNode));
        
    end

    % Calculate the mean over the movement and sum up the muscles
    % afterwards. Use therefore only positive power values.
    totalCEPower = sum( mean(max(0,powers), 2) );
   
    % Calculate powerd
%   powerd = -mean(min(0,powerdata));
%   meanpower_neg = sum(mean(min(0,powerData)));
%   powerd = meanpower/0.25 - meanpower_neg/1.2;
end
