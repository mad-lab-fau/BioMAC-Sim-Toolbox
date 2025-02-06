%======================================================================
%> @file cotTerm.m
%> @brief Collocation function to compute cost of transport from states
%> @details
%> Details: Collocation::cotTerm()
%>
%> @author Anne Koelewijn
%> @date October, 2019
%======================================================================

%======================================================================
%> @brief Collocation function to compute cost of transport of the muscles from states vector
%>
%> @param obj            Collocation class object
%> @param option         String parsing the demanded output: 'objval' or 'gradient'
%>                       (or 'init' for initialization)
%> @param X              Double array: State vector containing at least 'states', 'controls', 'dur',
%>                       and 'speed' of the model
%> @param name           (optional) name of the model to be used, defaults to lichtwark
%> @param epsilon        (optional) level of nonlinearity, defaults to 10^-3
%> @param exponent       (optional) Positive integer: Exponent of energy rate in optimization
%======================================================================
function output = cotTerm(obj,option,X, name, epsilon,exponent)

%% check input parameter
% => Here in the cotTerm, it is not really worth initializing variables
if strcmp(option,'init')
    if ~isfield(obj.idx,'states') ||  ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur')  || ~isfield(obj.idx,'speed')
        error('Model states, controls, duration, and speed need to be stored in state vector X.')
    end

    % Return a dummy value
    output = NaN;
    return;
end

if nargin < 4
    name = 'lichtwark';
end
if nargin < 5
    epsilon = 10^-3;
end
if nargin < 6
    exponent = 1;
end

%% compute demanded output
if strcmp(option,'objval') % objective value
    [~,~,~,~,output] = obj.getMetabolicCost(X, name, 1, epsilon, exponent); 
elseif strcmp(option,'gradient') % gradient 
    [~,~,~,~,~,~,output] = obj.getMetabolicCost(X, name, 1, epsilon, exponent);
else
    error('Unknown option');
end
end