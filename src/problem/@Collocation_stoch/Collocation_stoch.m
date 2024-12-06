% ======================================================================
%> @file Collocation_stoch.m
%> @brief Matlab class for defining an optimal control of musculoskeletal
%> model
%
%> @author Eva Dorschky, Marlies Nitschke
%> @date November, 2017
% ======================================================================

% ======================================================================
%> @brief The class for defining an optimal control of musculoskeletal
% ======================================================================
classdef Collocation_stoch < Collocation
    
    properties (SetAccess = private)
        randvals
    end
    
    methods
        %======================================================================
        %> @brief Constructor setting default Collocation options.
        %>
        %> @param   model     Model class object
        %> @param   nNodes    (optional) Double: Number of collocation nodes 
        %> @param   Euler     (optional) String: Discretization method: backward euler 'BE' or midpoint euler 'ME' 
        %> @param   logfile   (optional) String: Logfile including the path. The input can be skipped with an empty string.
        %> @param   plotLog   (optional) Bool: If true, plot log in call of Collocation.objfun() (stick figure, objective, constrains) (default: 0)
        %>
        %> @retval    obj   Collocation class object
        %======================================================================
        function [obj] = Collocation_stoch(randvals, varargin)
            % call Gait2dc constructor
            obj = obj@Collocation(varargin{:});
            obj.randvals = randvals;
            
        end
    end
end
