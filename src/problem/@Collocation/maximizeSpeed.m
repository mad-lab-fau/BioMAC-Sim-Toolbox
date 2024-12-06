%======================================================================
%> @file maximizeSpeed.m
%> @brief Collocation function to maximize squared speed
%> @details
%> Details: Collocation::maximizeSpeed()
%>
%> @author Christopher Loeffelmann, Markus Gambietz
%> @date May, 2021
%======================================================================

%======================================================================
%> @brief
%> Computes the squared speed multiplied by -1
%>
%> @todo This function is not yet working if there are multiple entries for
%> the speed (i.e. in the 3D model). simVar should be norm(X(speedidx)).
%>
%> @param  obj           Collocation class object
%> @param  option        String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param  X             Double array: State vector containing at least 'speed'
%> @param  mode          'squared' (default) for -x^2; 'log' for convex function -log(speed). 
%> @retval output        Objective values for input option 'objval' or vector
%>                       with gradient
%======================================================================
function output = maximizeSpeed(obj, option, X, mode)

%% check input parameter
% => Here in maximizeSpeed, it is not really worth initializing variables
if strcmp(option,'init')
    if ~isfield(obj.idx,'speed') % check whether speed is stored in X
        error('Model speed is not stored in state vector X.')
    end

    % Return a dummy value
    output = NaN;
    return;
end
if nargin == 3
    mode = 'squared';
end

%% compute demanded output
simVar = X(obj.idx.speed);

if strcmp(option,'objval') % objective value
    if strcmp(mode,'squared')
        output = -(simVar).^2;
    elseif strcmp(mode,'log')
        output = -log(simVar);
    else
        error('Collocation.maximizeSpeed: Unknown mode')
    end
elseif strcmp(option,'gradient') % gradient
    output = zeros(size(X));
    if strcmp(mode,'squared')
        output(obj.idx.speed) = -2*(simVar);
    elseif strcmp(mode,'log')
        output(obj.idx.speed) = -(1/simVar);
    else
        error('Collocation.maximizeSpeed: Unknown mode')
    end
else
    error('Unknown option');
end
end


