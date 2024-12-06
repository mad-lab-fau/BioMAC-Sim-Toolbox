% ======================================================================
%> @file Result.m
%> @brief Matlab class for defining an optimization result
%>
%> @author Marlies Nitschke
%> @date January, 2018
% ======================================================================

% ======================================================================
%> @brief The class for defining an optimization result
%>
%> @details
%> This class was implemented to ensure a consistent data structure to save
%> results. 
%>
%> This class contains methods which are generally applicable to evaluate
%> results (e.g. writeMovie()). Project specific evaluation methods
%> should be saved in our personal scripts folder!
% ======================================================================
classdef Result < handle
    properties(SetAccess = private)       
        %> Object of the class Problem
        problem         
        %> Object of the class Solver
        solver  
        %> Double array: Solution of the problem
        X 
        %> Bool: True if optimal solution was found (default: false)
        converged = false;
        %> Struct: Further information of the solver
        info
        %> String: User name of the computer
        userName
        %> String: Name of the computer
        computerName   
        %> Datetime: Creation date with 'dd-MMM-yyyy HH:mm:ss'
        creationTime
        %> String: Git hash of the current commit
        gitHashString
        %> String: Url of git repository
        gitURL
    end
    
    properties
        %> String: Filename with path, without extension
        filename
        %> String: Personal comment this result
        comment
    end
    
    methods
        %======================================================================
        %> @brief Default constructor setting default Result object
        %>
        %> @param   problem     Object of the class Problem
        %> @param   solver      Object of the class Solver
        %> @retval  obj         Result class object
        %======================================================================
        function [obj] = Result(problem, solver)
            % add general information
            [~, obj.userName] = system('whoami');
            try
                obj.computerName = char(java.net.InetAddress.getLocalHost.getHostName);
            catch
                warning('Could not get computer name.');
                obj.computerName = '';
            end
            obj.creationTime = datetime('now','Format','dd-MMM-yyyy HH:mm:ss'); 
            try
                [~,gitHashString] = system('git rev-parse HEAD');
            catch
                warning('Could not get git hash');
                gitHashString = '';
            end
            obj.gitHashString = gitHashString;
            try
                [~,gitURL] = system('git config --get remote.origin.url');
            catch
                warning('Could not get git url');
                gitURL = '';
            end
            obj.gitURL = gitURL;
            
            % Problem
            if isa(problem, 'Problem')
                obj.problem = problem;
            else
               error('Result:Result', 'The first input must be an object of the class Problem.'); 
            end
            
            % Solver
            if isa(solver, 'Solver')
                obj.solver = solver;
            else
                error('Result:Result', 'The second input must be an object of the class Solver.');
            end

        end
        
        %======================================================================
        %> @brief Function to set the results after solving
        %>
        %> @param   obj          Result class object
        %> @param   X            Double array: Solution of the problem
        %> @param   converged    Bool: True if optimal solution was found
        %> @param   info         (optional) Struct: Further information
        %======================================================================
        function setResult(obj, X, converged, info)
            obj.X = X;
            obj.converged = converged;
            
            if nargin > 3
               obj.info = info; 
            end
            
        end
        
        %======================================================================
        %> @brief Function to save the Result object 
        %>
        %> @details
        %> The Result object will be saved in a variable called 'result'
        %> and saved using the given filename.
        %>
        %> @param   obj          Result class object
        %> @param   filename     (optional) String: Filename with path, without extension
        %>                       (default: Used Result.filename)
        %======================================================================
        function save(obj, filename)
            if nargin > 1       
                obj.filename = filename;
            elseif isempty(obj.filename)
                error('Result:save', 'The filename property is emtpy.');
            end
                
            result = obj;
            save(obj.filename, 'result');   
            
        end
    end
end

