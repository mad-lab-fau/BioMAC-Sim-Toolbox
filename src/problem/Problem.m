% ======================================================================
%> @file Problem.m
%> @brief Matlab class (abstract) for defining an optimization problem
%>
%> @author Eva Dorschky
%> @date November, 2017
% ======================================================================

% ======================================================================
%> @brief The abstract class for defining an optimization problem
% ======================================================================
classdef (Abstract) Problem < handle
    properties
        %> String: Name of result file
        name;         
        %> Struct: Intial guess for x
        initialguess;         
        %> Double array: Lower bound of state vector
        X_lb;           
        %> Double array: Upper bound of state vector
        X_ub;           
        %> Double array: Lower bound of constaints
        c_lb;      
        %> Double array: Upper bound of constraints
        c_ub;           
    end
    
    methods (Abstract)
        % ======================================================================
        %> @brief Abstract function to computing objective value at state vector X
        % ======================================================================
        f = objfun(obj,X) 
        
        % ======================================================================
        %> @brief Abstract function to computing gradient at state vector X
        % ======================================================================
        g = gradient(obj,g) %gradient function handle
        
        % ======================================================================
        %> @brief Abstract function to compute constraint violations at state vector X
        % ======================================================================
        c = confun(obj,X) %constaint function handle
        
        % ======================================================================
        %> @brief Abstract function to compute jacobian at state vector X
        % ======================================================================
        J = jacobian(obj,X) %jacobian function handle
        
        % ======================================================================
        %> @brief Abstract function providing the jacobian pattern
        % ======================================================================
        J = jacstructure(obj) %jacobian structure function handle
        
        % ======================================================================
        %> @brief Abstract function to report a current state vector of the problem
        % ======================================================================
        [conString, texString] = report(obj, X, settings, style, resultFilename)
        
    end
    
    methods
        % ======================================================================
        %> @brief Function to test the whole problem which was defined
        % ======================================================================
        function derivativetest(obj, xr)
            if ismethod(obj,'initObjectives')
                obj.initObjectives;
            end

            obj.computeJpattern;
            
            if nargin == 1
                xr = obj.X_lb + (obj.X_ub - obj.X_lb) .* rand(size(obj.X_lb));
            end
            f = obj.objfun(xr);
            g = obj.gradient(xr);
            gnum = zeros(size(xr));
            c = obj.confun(xr);
            J = obj.jacobian(xr);
            Jnum = zeros([size(c,1),numel(xr)]);
            hh = 1e-7;
            tic
            for ih = 1:numel(xr)
                if (toc > 3.0)
                    tic
                    fprintf('Checking derivatives: %4.1f%% done...\n', 100*ih/numel(xr));
                end
                xrsave = xr(ih);
                xr(ih) = xr(ih) + hh;
                fhh = obj.objfun(xr);
                gnum(ih) = (fhh-f)/hh;
                chh = obj.confun(xr);
                Jnum(:,ih) = (chh-c)/hh;
                xr(ih) = xrsave;
            end
            fprintf('Checking gradient...\n'); 	obj.matcompare(g, gnum);
            fprintf('Checking Jacobian...\n'); 	obj.matcompare(J, Jnum);
            J = full(J);
            keyboard
        end        
    end
     methods (Static)
         
        % ======================================================================
        %> @brief Helperfunction to compare matrices and mind maximum abs error
        %> @details
        %> Used for derivativetest()
        % ======================================================================
        function matcompare(a,b)
            [maxerr,irow] = max(abs(a-b));
            [maxerr,icol] = max(maxerr);
            irow = irow(icol);
            fprintf('Max. difference: %14.10f at %d %d (%14.10f vs. %14.10f)\n', ...
                maxerr, irow, icol, full(a(irow,icol)),full(b(irow,icol)));
        end
        
    end
    
end

