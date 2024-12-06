% ======================================================================
%> @file IPOPT.m
%> @brief Matlab class to call the IPOPT optimization software
%>
%> @author 
%> @date October, 2017
% ======================================================================

% ======================================================================
%> @brief Class to call the IPOPT optimization software
%>
%> @details
%> https://www.coin-or.org/Ipopt/documentation/node40.html
%>
%> You can run a derivative test before solving by:
%> @code 
%> solver.setOptionField('derivative_test', 'first-order');
%> solver.solve(problem);
%> @endcode
% ======================================================================
classdef IPOPT < Solver
    
    properties (SetAccess = private)
        %> Struct containing the options for the IPOPT.
        %> Detailed settings of the options of the IPOPT software are documented
        %> here: https://www.coin-or.org/Ipopt/documentation/node40.html.
        %> Some default values are set within the constructor.
        solverOptions;
        ipoptOptions;
    end
    
    methods
        %======================================================================
        %> @brief Default constructor setting default solver options.
        %>
        %> @details
        %> Use the setting function setOptionField() to change the setting
        %> values and do not change the default settings here!
        %>
        %> @retval    obj         IPOPT class object
        %======================================================================
        function [obj] = IPOPT()
            
            % Set default options for termination
			obj.solverOptions.tol = 0.0002; % Desired convergence tolerance.
            obj.solverOptions.max_iter = 4000; % Maximum number of iterations. 
            obj.solverOptions.constr_viol_tol = 0.001; % Absolute tolerance on the constraint violation. 
            obj.solverOptions.compl_inf_tol = 0.001; % Absolute tolerance on the complementarity
            obj.solverOptions.acceptable_tol = 10^(-6); % Determines which (scaled) overall optimality error is considered to be "acceptable." 
            
            % Set initialization options
            obj.solverOptions.bound_frac = 0.001; % Desired minimum relative distance from the initial point to bound. 
            obj.solverOptions.bound_push = 0.001; % Desired minimum absolute distance from the initial point to bound. 
            
            % Set solving options
			obj.solverOptions.hessian_approximation = 'limited-memory'; % Indicates what Hessian information is to be used.
			obj.solverOptions.check_derivatives_for_naninf = 'no'; % Indicates whether it is desired to check for Nan/Inf in derivative matrices.
            obj.solverOptions.mu_strategy = 'adaptive'; % Determines which barrier parameter update strategy is to be used.
            obj.solverOptions.linear_solver = 'mumps'; % Linear solver used for step computations. 

            % Set output options
			obj.solverOptions.print_level = 5; % Output verbosity level.
            obj.solverOptions.print_timing_statistics = 'yes'; %If selected, the program will print the CPU usage (user time) for selected tasks.
        end
        
        %======================================================================
        %> @brief Function to set a new field or change a field in the
        %> solverOptions struct.
        %>
        %> @details 
        %> Example of use: 
        %> @code
        %> obj.setOptionField('tol', 0.001)
        %> @endcode
        %>
        %> @param    obj         IPOPT class object
        %> @param    field       Fieldname to set or change
        %> @param    value       Value to set
        %======================================================================
        function setOptionField(obj, field, value)
            obj.solverOptions.(field) = value;            
        end
        
        %======================================================================
        %> @brief Function to prepare lambda, zl, and zu in case a warm
        %> start will be used. They should have the length of the problem
        %>
        %> @details 
        %> Example of use: 
        %> @code
        %> obj.setWarmstartParams(obj, 1:nCon, 1:nDof)
        %> @endcode
        %>
        %> @param    obj         IPOPT class object
        %> @param    problem     Problem: Numerical problem which should be solved
        %> @param    nConInd     (optional) Indices of the lambdas that should be used
        %> @param    nOptInd     (optional) Indices for zl and zu that should be used
        %======================================================================
        function setWarmstartParams(obj, problem, nConInd, nOptInd)
            obj.solverOptions.warm_start_init_point = 'yes';
            if nargin == 2
                nConInd = 1:problem.nConstraints;
                nOptInd = 1:problem.nVars;
            end
            obj.solverOptions.warm_start_bound_frac = 1e-16;
            obj.solverOptions.warm_start_bound_push = 1e-16;
            obj.solverOptions.warm_start_mult_bound_push = 1e-16;
            obj.solverOptions.warm_start_slack_bound_frac = 1e-16;
            obj.solverOptions.warm_start_slack_bound_push = 1e-16;       
            obj.ipoptOptions.lambda = problem.initialguess.info.lambda(nConInd);
            obj.ipoptOptions.zl     = problem.initialguess.info.zl(nOptInd);
            obj.ipoptOptions.zu     = problem.initialguess.info.zu(nOptInd);
        end
        
        %======================================================================
        %> @brief Function to solve the optimization problem.
        %>
        %> @param    obj         IPOPT class object
        %> @param    problem     Object of the class Problem
        %> @retval   result      Result class object containing the solution
        %======================================================================
        function [result] = solve(obj,problem)
            if ~isa(problem,'Problem')
                error('IPOPT:solver', '''problem'' must be an object of the Problem class.');
            end
           
            % set funcs handles
            funcs.objective         = @problem.objfun;
			funcs.gradient          = @problem.gradient;
			funcs.constraints       = @problem.confun;
			funcs.jacobian          = @problem.jacobian;
			funcs.jacobianstructure = @problem.jacstructure;
            % set bounds and options
            options = obj.ipoptOptions;
			options.lb = problem.X_lb;
			options.ub = problem.X_ub;
			options.cl = problem.c_lb;
			options.cu = problem.c_ub;
            options.ipopt = obj.solverOptions;
            if ~strcmp(problem.name , 'untitled') && ~startsWith(pwd, '/home/hpc/')
                options.ipopt.output_file = [problem.name '_ipopt.log'];
            end
            
            % decode the initial guess, do warm start if required
            % information is included
            
            % initial guess is a struct with X and info from previous result
            % this indicates that warm start is possible for this problem
            X0 = problem.initialguess.X;
            
            % intialize the objectives
            if ismethod(problem,'initObjectives')
                problem.initObjectives;
            end

            % compute Jpattern before solving
            problem.computeJpattern;
            
            % run IPOPT
            timerSolve = tic;
            [X, info] = ipopt(X0,funcs,options);   
            info.wallTime = toc(timerSolve);
            
            % try to get IPOPT version
            if isfield(options.ipopt, 'output_file')
                info.IPOPTVersion = IPOPT.getIPOPTVersion(options.ipopt.output_file);
            end
            
            % create Result object
            result = Result(problem, obj);
            converged = info.status == 0;
            result.setResult(X, converged, info);
        end
        
        
        %======================================================================
        %> @brief Function to report the solver performance and settings
        %>
        %> @param    obj               IPOPT class object
        %> @param    doReportSettings  Boolean: If true the solver settings will
        %>                             also be reported. (default if empty: 0)
        %> @param    solvingInfo       (optional) Struct: Infos of the solver which were
        %>                             saved in the result in solver.solve(). If it is not 
        %>                             given or empty, the solver status will not be reported.
        %> @param    resultFilename    (optional) String: Filename to save report in a LaTeX 
        %>                             document. If the filename is not given the report will 
        %>                             not be saved.
        %> @param    title             (optional) String: Title of the document
        %>                             (default: 'Solver Report')
        %> @retval   conString         String: Formated text to be printed in the console 
        %>                             containing the solver information and the solver options.
        %> @retval   texString         String: Formated LaTeX text containing the
        %>                             solver information and the solver options in two
        %>                             separate headers.
        %======================================================================
        function [conString, texString] = report(obj,doReportSettings,solvingInfo,resultFilename, title)
            
            % Check which reports should be made
            if nargin < 2 || isempty(doReportSettings)
               doReportSettings = 0; 
            end
            if nargin < 3 || isempty(solvingInfo)
               doReportStatus = 0;
            else
               doReportStatus = 1;
            end
            
            % 1. Solver Status
            % see https://github.com/coin-or/Ipopt/blob/master/Ipopt/src/Interfaces/IpReturnCodes_inc.h
            % and https://github.com/coin-or/Ipopt/blob/master/Ipopt/src/Interfaces/IpIpoptApplication.cpp
            if doReportStatus
                % Define dictionary for solver status
                statusDict.ID = [0:6, -1, -2, -3, -4, -10, -11, -12, -13, -100, -101, -102, -199]';
                statusDict.MsgID = {'Solve_Succeeded', ...                    % 0
                                    'Solved_To_Acceptable_Level', ...         % 1
                                    'Infeasible_Problem_Detected', ...        % 2
                                    'Search_Direction_Becomes_Too_Small', ... % 3
                                    'Diverging_Iterates', ...                 % 4
                                    'User_Requested_Stop', ...                % 5
                                    'Feasible_Point_Found', ...               % 6
                                    'Maximum_Iterations_Exceeded', ...        % -1
                                    'Restoration_Failed', ...                 % -2
                                    'Error_In_Step_Computation', ...          % -3
                                    'Maximum_CpuTime_Exceeded', ...           % -4
                                    'Not_Enough_Degrees_Of_Freedom', ...      % -10
                                    'Invalid_Problem_Definition', ...         % -11
                                    'Invalid_Option', ...                     % -12
                                    'Invalid_Number_Detected', ...            % -13
                                    'Unrecoverable_Exception', ...            % -100
                                    'NonIpopt_Exception_Thrown', ...          % -101
                                    'Insufficient_Memory', ...                % -102
                                    'Internal_Error'}';                       % -199
                statusDict.Msg =   {'Optimal Solution Found', ...
                                    'Solved To Acceptable Level', ...
                                    'Converged to a point of local infeasibility. Problem may be infeasible.', ...
                                    'Search Direction is becoming Too Small.', ...
                                    'Iterates diverging; problem might be unbounded.', ...
                                    'Stopping optimization at current point as requested by user.', ...
                                    'Feasible point for square problem found.', ...
                                    'Maximum Number of Iterations Exceeded.', ...
                                    'Restoration Failed!', ...
                                    'Error in step computation (regularization becomes too large?)!', ...
                                    'Maximum CPU time exceeded.', ...
                                    'Problem has too few degrees of freedom.', ...
                                    'Invalid_Problem_Definition', ...
                                    'Invalid option encountered.', ...
                                    'Invalid number in NLP function or derivative detected.', ...
                                    'Some uncaught Ipopt exception encountered.', ....
                                    'Unknown Exception caught in ipopt', ...
                                    'Not enough memory.', ...
                                    'INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.'}';
                statusDict = struct2table(statusDict);
                
                % Get solver status
                idxStatus = find(statusDict.ID == solvingInfo.status);
                if ~isempty(idxStatus) % message is known
                    msgIDString = statusDict.MsgID{idxStatus};
                    msgString   = statusDict.Msg{idxStatus};
                else % We do not understand the message
                    msgIDString = '';
                    msgString   = '';
                end

                % Get number of iterations
                iterString = sprintf('%u', solvingInfo.iter);

                % Get duration of solver
                cpuString = sprintf('%s (HH:MM:SS)', datestr(seconds(solvingInfo.cpu), 'HH:MM:SS'));
                if isfield(solvingInfo, 'wallTime') % Check to further support old solver objects
                    wallString = sprintf('%s (HH:MM:SS)', datestr(seconds(solvingInfo.wallTime), 'HH:MM:SS'));
                else
                    wallString = NaN;
                end

                % Fill the information in a formreated string to be printed in
                % the console
                conStringStatus = [newline, ...
                                   sprintf('** Solver Status ** \n'), ...
                                   sprintf('   - Status ID: %s \n', msgIDString), ...
                                   sprintf('   - Status Message: %s \n', msgString), ...
                                   sprintf('   - Number of iterations: %s \n', iterString), ...
                                   sprintf('   - CPU time: %s \n', cpuString), ...
                                   sprintf('   - Wall time: %s \n', wallString), ...
                                   newline]; 

                % Fill the information in a formated latex string 
                texStringStatus = [newline, ...
                                   sprintf('\\subsection{Solver Status} \n'), ...
                                   sprintf('\\begin{itemize} \n'), ...
                                   sprintf('   \\item Status ID: %s \n', strrep(msgIDString, '_', '\textunderscore ')), ...
                                   sprintf('   \\item Status Message: %s \n', msgString), ...
                                   sprintf('   \\item Number of iterations: %s \n', iterString), ...
                                   sprintf('   \\item CPU time: %s \n', cpuString), ...
                                   sprintf('   \\item Wall time: %s \n', wallString), ...
                                   sprintf('\\end{itemize} \n')]; 
            else
                conStringStatus = '';
                texStringStatus = '';
            end
                     
            % 2. Solver settings
            if doReportSettings
                % Fill the beginning of the strings
                conStringSettings = [newline, ...
                                    sprintf('** Solver Settings ** \n'), ...
                                    sprintf('   - Solver: %s \n', class(obj))];
                texStringSettings = [newline, ...
                                    sprintf('\\subsection{Solver Settings} \n'), ...
                                    sprintf('\\begin{itemize} \n'), ...
                                    sprintf('   \\item Solver: %s \n', class(obj))];

                % Iterate through the option fields and add it to the strings
                optNames = fieldnames(obj.solverOptions);
                for iOpt = 1 : length(optNames)
                    optName = optNames{iOpt};
                    optValue = obj.solverOptions.(optName);
                    if isnumeric(optValue)
                        optValue = num2str(optValue);
                    end

                    conStringSettings = [conStringSettings, ...
                                         sprintf('   - %s: %s \n', optName, optValue)];
                    texStringSettings = [texStringSettings, ...
                                         sprintf('   \\item %s: %s \n', strrep(optName, '_', '\textunderscore '), optValue)];
                end
                
                % Add the end to the strings
                optDefaut = 'For all other options, default values were used.';
                conStringSettings = [conStringSettings, ...
                                    sprintf('   - %s \n', optDefaut), ...
                                    newline];  
                texStringSettings = [texStringSettings, ...
                                    sprintf('   \\item %s \n', optDefaut), ...
                                    sprintf('\\end{itemize} \n')]; 
            
            else
                conStringSettings = '';
                texStringSettings = '';
            end
            
  
            % 3. Combine both sections
            conStringHeader = [newline, ...
                               sprintf('**** Solver **** \n')];
            texStringHeader = [newline, ...
                               sprintf('\\section{Solver} \n')];
            conString = [conStringHeader, conStringStatus, conStringSettings];
            texString = [texStringHeader, texStringStatus, texStringSettings];
            
            % 4. Save report
            if nargin > 3
                % Get title of the doucment if it was not set
                if nargin < 5
                   title = 'Solver Report';
                end
                % Save the report
                saveLaTeXDocument(texString, [resultFilename '_SolverInfo'], title);
                % Create pdf file
                [folderName,fileName]=fileparts(resultFilename);
                fileName=[fileName,'_SolverInfo'];
                % check whether the tex file has no errors so the pdf can be created
                % for sure
                try
                    callPdflatex (folderName,fileName);
                catch
                    errorMessage = sprintf('Could not make the pdf since the associated tex file has errors. Check the file:\n %s_SolverInfo.tex\n', resultFilename);
                    warning (errorMessage);
                end
            end
               
        end
        
    end
    
    methods(Static)
        %======================================================================
        %> @brief Function to get the IPOPT version of the ebertolazzi toolbox
        %>
        %> @details It is not ideal to obtain the IPOPT version from the
        %> log file, but we did not yet find a better solution.
        %>
        %> @param    obj         IPOPT class object
        %> @param    output_file String: Filename of output file
        %> @retval   version     String: Version of IPOPT
        %======================================================================
        function version = getIPOPTVersion(output_file)
            
            % Load file
            text = fileread(output_file);

            % Find version in log file
            textStart = 'This is Ipopt version ';
            idxStart = strfind(text,textStart)+length(textStart);
            textEnd = ', running with';
            idxEnd = strfind(text, textEnd)-1;
            version = text(idxStart:idxEnd);
            
        end
    end
    
end

