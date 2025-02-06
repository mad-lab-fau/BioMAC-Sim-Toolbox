% ======================================================================
%> @file Solver.m
%> @brief Matlab class (abstract) a solver inherits from
%>
%> @author 
%> @date October, 2017
% ======================================================================

% ======================================================================
%> @brief The abstract class describes the basics of a solver
% ======================================================================
classdef (Abstract) Solver < handle
    
    properties (Abstract, SetAccess = private)
        %> Struct containing the options for the solver
        solverOptions
    end
    

    methods (Abstract)
        % ======================================================================
        %> @brief Abstract function to solve the optimization problem.
        % ======================================================================
        solve(obj,problem);
        
        % ======================================================================
        %> @brief Abstract function to report the solver performance and settings.
        % ======================================================================
        [conString, texString] = report(obj,doReportSettings,solvingInfo,resultFilename, title);
    end
    
end

