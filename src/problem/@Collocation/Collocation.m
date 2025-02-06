% ======================================================================
%> @file Collocation.m
%> @brief Matlab class for defining an optimal control of musculoskeletal
%> model
%
%> @author Eva Dorschky, Marlies Nitschke
%> @date November, 2017
% ======================================================================

% ======================================================================
%> @brief The class for defining an optimal control of musculoskeletal
% ======================================================================
classdef Collocation < Problem

    properties
        %> Model: Object of class Model.m
        model;
    end

    properties (SetAccess = private)
        %> Double: Number of collocation nodes
        nNodes;		
        %> String: Midpoint euler 'ME' or backward euler 'BE' or semi-implicit euler 'SIE'
        Euler;	
        %> Double: Number of variables that are optimized
        nVars;     		
        %> Struct: Indices for position of variables in state vector X
        idx
        
        %> Double: Total number of constraints
        nConstraints;     	
        %> Double: Number of non-zero entries in jacobian pattern
        Jnnz;			
        %> Double matrix: Jacobian pattern
        Jpattern;				
        
        %> Struct: Information for logging objectives and constrained violations
        log; 
        %> Bool: If true, plot log in call of Collocation.objfun() (stick figure, objective, constrains) (default: 0)
        plotLog = 0;
        
        %> Struct array: Ojective terms
        objectiveTerms
        %> Struct array: Constraints
        constraintTerms
    end
    
    properties(Dependent, SetAccess = private)
        %> Bool: True if the periodicityConstraint() was used.
        isPeriodic;
        %> Bool: True if the periodicityConstraint() was used and there sym is true.
        isSymmetric;
        %> Double: Number of collocation nodes which define the duration of the simulation.
        %> If Collocation.isPeriodic is true, it is Collocation.nNodes+1
        %> since we are simulating an additional node.
        %> If Collocation.isPeriodic is false, it is Collocation.nNodes
        %> since we do not simulate an additional node.
        %> This has to be taken into account when specifying the duration
        %> of the simulation as h=T/(nNodesDur-1).
        nNodesDur;
    end
    
    properties(Transient, SetAccess = private)
        %> Struct: Contains values extracted during initialization of objective functions.
        %> This is used to avoid the extraction in every call of the objective in each
        %> iteration thus to improve performance. Since it is a transient property it
        %> will not be saved when saving the object.
        objectiveInit;
    end

    methods
        %======================================================================
        %> @brief Constructor setting default Collocation options.
        %>
        %> @param   model     Model class object
        %> @param   nNodes    (optional) Double: Number of collocation nodes 
        %> @param   Euler     (optional) String: Discretization method: backward euler 'BE' or midpoint euler 'ME' or semi-implicit euler 'SIE
        %> @param   logfile   (optional) String: Logfile including the path. The input can be skipped with an empty string.
        %> @param   plotLog   (optional) Bool: If true, plot log in call of Collocation.objfun() (stick figure, objective, constrains) (default: 0)
        %>
        %> @retval    obj   Collocation class object
        %======================================================================
        function [obj] = Collocation(model,nNodes,Euler,logfile,plotLog)
            tic
            % initialize input parameter
            if nargin < 1
                error('Model must be defined.');
            elseif ~isa(model,'Model')
                error('The first input argument should be of type obj.model.');
            else
                obj.model = model;
            end
            
            if nargin < 2
                obj.nNodes = 50;
            elseif nNodes <= 0
                error('The number of collocation nodes is invalid.');
            elseif nNodes > 1600
                answer = input('The number of collocation nodes is very high and the simulation will likely be slow. Are you sure you want to continue? Press y/n',"s");
                if strcmp(answer,'y')
                    obj.nNodes = nNodes;
                else
                    return
                end
            else
                obj.nNodes = nNodes;
            end
            
            if nargin < 3
                obj.Euler = 'BE';
            elseif strcmp(Euler,'BE') || strcmp(Euler,'ME') || strcmp(Euler,'SIE')
                obj.Euler = Euler;
            else
                error('The discretization method is not in the valid range.');
            end
            
            if nargin < 4 || isempty(logfile)
                obj.name  = 'untitled';
                obj.log.printlog = 0;
            else
                obj.name = logfile;
                obj.log.printlog = 1;
            end
            
            if nargin < 5
                obj.plotLog = 0;
            else
               obj.plotLog = plotLog; 
            end
            
            obj.nVars = 0; % intialize number of values in state vector X
            obj.nConstraints = 0; % intialize number of constraints
            obj.initialguess.X = []; %initialize initial guess
            
            % initialize parameters for logging 
            obj.log.it = 0; % intialize iteration to zero
            if obj.log.printlog
                obj.log.file = [obj.name '.log']; % define log file saving objective values and constraint violation at iteration it
                fid = fopen(obj.log.file,'w');
                if fid == -1
                    error('Could not open %s', obj.log.file);
                end
                fclose(fid);
            end
            
        end
        
        %======================================================================
        %> @brief Function computing objective value at state vector X
        %>
        %> @param   obj Collocation class object
        %> @param   X   Double array: State vector with variables that are optimized
        %>
        %> @retval  f   Double array: Objective value at state vector X 
        %======================================================================
        function f = objfun(obj,X)
            f = 0;
            for iObj = 1:length(obj.objectiveTerms)
                W =  obj.objectiveTerms(iObj).weight; % weight of objective term
                name = obj.objectiveTerms(iObj).name; % function name of objective term
                varargin = obj.objectiveTerms(iObj).varargin; % optional input parameters of objective term
                f_tmp = obj.(name)('objval',X,varargin{:}); % objective value of objective term 
                obj.objectiveTerms(iObj).unweightedValue = f_tmp; % save current unweighted value
                obj.objectiveTerms(iObj).weightedValue = W*f_tmp; % save current weighted value
                obj.objectiveTerms(iObj).weightedValueHist(end+1, 1) = obj.objectiveTerms(iObj).weightedValue; % save current weighted value in history
                
                f = f + obj.objectiveTerms(iObj).weightedValue; % sum up objective terms
            end
            
            % logging
            drawnow('update')  % to make the IPOPT output appear immediately, use drawnow instead of pause
            obj.log.it = obj.log.it + 1;

            if toc > 100 % Every 100 seconds, write something on the log file and
                obj.printlog;
                obj.plotlog(X); 
                tic;
            end
            
        end
        
        %======================================================================
        %> @brief Function initializing the objective functions
        %>
        %> @param   obj Collocation class object
        %======================================================================
        function initObjectives(obj)
            obj.objectiveInit = struct();
            for iObj = 1:length(obj.objectiveTerms)
                name = obj.objectiveTerms(iObj).name; % function name of objective term
                varargin = obj.objectiveTerms(iObj).varargin; % optional input parameters of objective term
                obj.(name)('init',[],varargin{:}); % objective value of objective term
            end
        end
        
        %======================================================================
        %> @brief Function computing gradient at state vector X
        %>
        %> @param   obj Collocation class object
        %> @param   X   Double array: State vector with variables that are optimized
        %>
        %> @retval  g   Double array: Gradient vector at state vector X 
        %======================================================================
        function g = gradient(obj,X)
            g = zeros(size(X));
            for iObj = 1:length(obj.objectiveTerms)
                W =  obj.objectiveTerms(iObj).weight; % weight of objective term
                name = obj.objectiveTerms(iObj).name; % function name of objective term
                varargin = obj.objectiveTerms(iObj).varargin; % optional input parameters of objective term
               
                g = g + W*obj.(name)('gradient',X,varargin{:}); % sum up gradients of objective terms
            end
            pause(0)
        end
        
        %======================================================================
        %> @brief Function to compute constraint violations at state vector X
        %>
        %> @param   obj Collocation class object
        %> @param   X   Double array: State vector with variables that are optimized
        %>
        %> @retval  c   Double array: Constraint violations at state vector X 
        %======================================================================
        function c = confun(obj,X)
            c = [];
            iConStart = 1;
            for iCon = 1:length(obj.constraintTerms)
                name = obj.constraintTerms(iCon).name; % function name of constraint 
                varargin = obj.constraintTerms(iCon).varargin; % optional input parameters of constraint
                c_tmp = obj.(name)('confun',X,varargin{:}); % constraint violations
                iConStop = iConStart + length(c_tmp) - 1; 
                obj.constraintTerms(iCon).normc = norm(min(c_tmp-obj.c_lb(iConStart:iConStop),0) + max(c_tmp-obj.c_ub(iConStart:iConStop),0)); % norm of constraint violations
                obj.constraintTerms(iCon).normcHist(end+1, 1) = obj.constraintTerms(iCon).normc;
                
                c = [c;c_tmp]; % put all constraint violations in one vector
                
                iConStart = iConStop + 1;
            end
            
            if ~isempty(find(isnan(c), 1))
                warning('NaN values found in constraint function');
                keyboard
            end
            
            if ~isempty(find(isnan(c), 1)) || ~isempty(find(isinf(c), 1))
                disp('Illegal value found in constraints')
                keyboard
            end
        end
        
        %======================================================================
        %> @brief Function to compute jacobian at state vector X
        %>
        %> @param   obj Collocation class object
        %> @param   X   Double array: State vector with variables that are optimized
        %>
        %> @retval  J   Sparse matrix: Jacobian at state vector X 
        %======================================================================
        function J = jacobian(obj,X)
            J = [];
            for iCon = 1:length(obj.constraintTerms)
                name = obj.constraintTerms(iCon).name; % function name of constraint
                varargin = obj.constraintTerms(iCon).varargin; % optional input parameters of constraint
                J = [J;obj.(name)('jacobian',X,varargin{:})]; % put all jacobians of constraint violations in one sparse matrix
            end
        end
        
        %======================================================================
        %> @brief Function providing the jacobian pattern
        %>
        %> @param   obj Collocation class object
        %>
        %> @retval  J   Double matrix: Pattern of jacobian matrix
        %======================================================================
        function J = jacstructure(obj)
            if isempty(obj.Jpattern)
                obj.computeJpattern;
            end
            J = double(obj.Jpattern);
        end
        
        %======================================================================
        %> @brief Function computing the jacobian pattern
        %>
        %> @param   obj Collocation class object
        %======================================================================
        function computeJpattern(obj)
            obj.Jnnz = 1; % temporary value needed for the first memory allocation
            obj.Jpattern = sparse(obj.nConstraints,length(obj.X_lb));
            nsame = 0;
            % sparsity pattern needs to be the same 10 times in a row
            while (nsame < 10)
                X = obj.X_lb + (obj.X_ub-obj.X_lb).*rand(size(obj.X_lb)); % a random vector of unknowns (uniformly distributed between bounds)
                J = obj.jacobian(X);
                newobj.Jpattern = double(J~=0) | obj.Jpattern ;			% add any new nonzeros that were just found
                if (nnz(newobj.Jpattern - obj.Jpattern) == 0)
                    nsame = nsame + 1;
                else
                    nsame = 0;
                end
                obj.Jpattern = newobj.Jpattern ;
                obj.Jnnz = nnz(obj.Jpattern);
                
                if obj.log.printlog
                    fprintf('Jacobian sparsity: %d nonzero elements out of %d (%5.3f%%).\n', ...
                        obj.Jnnz, obj.nConstraints*length(X), 100*obj.Jnnz/(obj.nConstraints*length(X)));
                end
            end
        end
        
        
        %======================================================================
        %> @brief Function logging objective value and constraint violations
        %>
        %> @param  obj Collocation class object
        %======================================================================
        function printlog(obj)
            if ~obj.log.printlog; return; end
            fid = fopen(obj.log.file,'a');
            fprintf(fid,['%d: ',obj.log.stringConfun,' ',obj.log.stringObjval,' \n'],obj.log.it,sum([obj.constraintTerms.normc]),obj.constraintTerms.normc,sum([obj.objectiveTerms.weightedValue]),obj.objectiveTerms.weightedValue); %\todo add constraints violation
            fclose(fid);
        end
        
        %======================================================================
        %> @brief Function plotting stick figure, objective values and constraint violations
        %>
        %> @param  obj Collocation class object
        %> @param  X   Double array: State vector with variables that are optimized
        %======================================================================
        function plotlog(obj, X)
            % return if we don't wanna plot it
            if ~obj.plotLog; return; end
            
            % set figure
            figure(1);
            set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.6]);
            clf;
            
            % plot stick figure
            subplot(1, 3, 1);
            obj.model.showStick(X(obj.idx.states));
            title('Stick Figure');
            
            % plot objective
            subplot(1, 3, 2); hold on;
            plot(sum([obj.objectiveTerms(:).weightedValueHist], 2), 'DisplayName', 'Sum');
            for iObj= 1 : length(obj.objectiveTerms)
                plot(obj.objectiveTerms(iObj).weightedValueHist, 'DisplayName', obj.objectiveTerms(iObj).name);
            end
            legend(gca, 'show', 'Location', 'southoutside');
            if ~isempty([obj.objectiveTerms(:).weightedValueHist])
                setPlotLim([sum([obj.objectiveTerms(:).weightedValueHist], 2), obj.objectiveTerms(:).weightedValueHist]);
            end
            title('Objective');
            xlabel('Function Evaluations');
            grid on; box on;
            
            % plot constraint violations
            subplot(1, 3, 3); hold on;
            plot(sum([obj.constraintTerms(:).normcHist], 2), 'DisplayName', 'Sum');
            for iCon = 1 : length(obj.constraintTerms)
                plot(obj.constraintTerms(iCon).normcHist, 'DisplayName', obj.constraintTerms(iCon).name);
            end
            legend(gca, 'show', 'Location', 'southoutside');
            if ~isempty([obj.constraintTerms(:).normcHist])
                setPlotLim([sum([obj.constraintTerms(:).normcHist], 2), obj.constraintTerms(:).normcHist]);
            end
            title('Constraint Violations');
            xlabel('Function Evaluations');
            grid on; box on;
            
            %======================================================================
            %> @brief Helperfunction setting the limits for objective and constraints plot
            %>
            %> @param  values  Double matrix: Values plotted (sum is first term) (nEvaluations x nTerms) 
            %======================================================================
            function setPlotLim(values)
                nEvaluations = size(values, 1);
                start = max(1, nEvaluations -1000);
                maxy = 2*mean(values(start:end, 1));
                miny = min(min(values(start:end, 2:end)));
                maxy = max(maxy, miny+1e-6);
                set(gca,'YLim',[miny maxy]);
                set(gca, 'XLim', [start max(2, nEvaluations)]);
            end
            
        end
        
        %======================================================================
        %> @brief Function adding optimization variables to the state vector 
        %>
        %> @details Variables that should be optimized can be added to the
        %> collocation problem. The 'name' how the variables are called in
        %> the objective functions/ constraints need to be defined as first
        %> input, e.g. X(obj.idx.(name)). The lower and upper bounds xmin < x < xmax need
        %> to be defined as second and third input. Moreover, an initial
        %> guess for x can be defined. If no initial guess is provided, the
        %> mid point between the lower and upper bound is used.
        %> 
        %> The variables can be a vector or a matrix defining
        %> the variables at multiple collocation nodes (nvarspernode x 1|Collocation.nNodes|Collocation.nNodes+1).
        %> @code
        %> problem.addOptimVar('states',repmat(model.states.xmin,1,N+1),repmat(model.states.xmax,1,N+1));
        %> problem.addOptimVar('controls',repmat(model.controls.xmin,1,N+1),repmat(model.controls.xmax,1,N+1));
        %> problem.addOptimVar('dur',0.1,1);
        %> problem.addOptimVar('speed',targetspeed,targetspeed);
        %> @endcode
        %> 
        %> @param   obj     Collocation class object
        %> @param   name    String: Name of variable
        %> @param   xmin    Double, vector or matrix: Lower bounds 
        %> @param   xmax    Double, vector or matrix: Upper bounds 
        %> @param   xinit   (optional) Vector or matrix: Initial guess 
        %======================================================================
        function addOptimVar(obj,name,xmin,xmax,xinit)
            %> @todo check if xmin, xmax, xinit have the correct size
            nParam = size(xmin,1)*size(xmin,2);
            obj.idx.(name) = reshape(obj.nVars+1:obj.nVars+nParam,size(xmin));
            obj.nVars = obj.nVars+nParam;
            obj.X_lb = [obj.X_lb;reshape(xmin,nParam,1)]; 
            obj.X_ub = [obj.X_ub;reshape(xmax,nParam,1)]; 
            if nargin < 5
                obj.initialguess.X = [obj.initialguess.X;(reshape(xmin,nParam,1)+reshape(xmax,nParam,1))/2]; % use mid point as default initial guess
            else
                obj.initialguess.X = [obj.initialguess.X;reshape(xinit,nParam,1)];
            end
        end
        
        %======================================================================
        %> @brief Function adding objective term to the optimization problem
        %>
        %> @details The objective terms are summed up in objective function using
        %> the function and weights which are defined as input.
        %> 
        %> @code
        %> problem.addObjective(@regTerm,W.reg)
        %> problem.addObjective(@effortTerm,W.effort,'volumeweighted',2)
        %> problem.addObjective(@trackGRF,W.track,trackingData)
        %> problem.addObjective(@trackAngles,W.track,trackingData)
        %> @endcode
        %> 
        %> @param   obj         Collocation class object
        %> @param   fh          Function handle: Method defined in the Collocation class folder
        %> @param   weight      Double: Weight of the objective term in the objective function
        %> @param   varargin    (optional) Input parameter for function handle
        %======================================================================
        function addObjective(obj,fh,weight,varargin)
            iObj = length(obj.objectiveTerms)+1;
            s = functions(fh);
            % add objective term to struct array
            obj.objectiveTerms(iObj).name = s.function; % save function name of objective term 
            obj.objectiveTerms(iObj).unweightedValue = NaN; % initialize unweighted value of objective term for logging
            obj.objectiveTerms(iObj).weightedValue = NaN; % initialize weighted value of objective term for logging (value of last evaluation)
            obj.objectiveTerms(iObj).weightedValueHist = []; % initialize weighted value of objective term for logging (values of all evaluations)
            obj.objectiveTerms(iObj).weight = weight; % save weight of objective term 
            obj.objectiveTerms(iObj).varargin = varargin; % save optional input paramter 
            if iObj == 1
                obj.log.stringObjval = ['f = %9.4f = %9.4f (',s.function,') '];
            else
                obj.log.stringObjval = [obj.log.stringObjval,'+ %9.4f (',s.function,') '];
            end
        end
        
        %======================================================================
        %> @brief Function adding constraint to the optimization problem
        %>
        %> @details 
        %> The constraints are called by the name of the function handle. The upper and lower
        %> bounds of the constraints can be set to zero for equalitiy
        %> constraints. The constraints can be a vector or a
        %> matrix to define constraints at a single or mulitple collocation nodes. (nconstraintspernode x 1|.nNodes).
        %> 
        %> @code
        %> problem.addConstraint(@dynamicConstraints,repmat(model.init.fmin,1,N),repmat(model.init.fmax,1,N))
        %> problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),sym)
        %> @endcode
        %> 
        %> @param   obj         Collocation class object
        %> @param   fh          Function handle: Method defined in Collocation
        %> @param   cmin        Double array: Lower bound on constraint
        %> @param   cmax        Double array: Upper bound on constraint
        %> @param   varargin    (optional) Input parameter for function handle
        %======================================================================
        function addConstraint(obj,fh,cmin,cmax,varargin)
            %> @todo check if cmin and cmax have the correct size
            nCon = size(cmin,1)*size(cmin,2);
            iCon = length(obj.constraintTerms)+1;
            s = functions(fh); %give information on handle
            % add constraint to struct array
            obj.constraintTerms(iCon).name = s.function; % save function name of constraint
            obj.constraintTerms(iCon).normc = NaN; % initialize norm of constraint violations for logging (value of last evaluation)
            obj.constraintTerms(iCon).normcHist = []; % initialize norm of constraint violations for logging (values of all evaluations)
            obj.constraintTerms(iCon).varargin = varargin; % save optional input paramter 
            obj.c_lb = [obj.c_lb;reshape(cmin,nCon,1)]; % lower bound of constraint
            obj.c_ub = [obj.c_ub;reshape(cmax,nCon,1)]; % upper bound of constraint
            obj.nConstraints = obj.nConstraints+nCon; % total number of constrant
            
            if iCon == 1
                obj.log.stringConfun = ['||c|| = %9.4f  = %9.4f (',s.function,') '];
            else
                obj.log.stringConfun = [obj.log.stringConfun,'+ %9.4f (',s.function,') '];
            end
        end
        
        %======================================================================
        %> @brief Function to define an initial guess
        %>
        %> @details 
        %> The initial guess is initialized by adding the state variables. This function can be used to
        %> overwrite the initial guess, e.g. by loading an initial guess of a previous
        %> result. The input parameter 'init' can be used to define how the
        %> initial guess should be computed. If init = 'mid', the midpoint
        %> between the lower and upper bounds of the state vector X are
        %> used. If init = 'random', X is intialized by drawing random
        %> values from a uniform distribution between the lower and upper
        %> bound. If init contains the path to a result file, the result
        %> file is loaded. If neccessary, the result is resampled to the
        %> size of the current collocation problem.
        %> 
        %> @code
        %> problem.makeinitialguess('mid');
        %> @endcode
        %> 
        %> @param   obj         Collocation class object
        %> @param   init        String: Parsing the demanded initialguess
        %>                      OR
        %>                      Double vector: X which fits exactly the size of
        %>                      the current problem (no tests and no resampling)
        %======================================================================
        function makeinitialguess(obj,init)

            if isnumeric(init)
                % init is X vector
                obj.initialguess.X = init;
                obj.initialguess.info.type = 'X';

            else
                % init is string
                switch init
                    case 'mid'
                        obj.initialguess.X = (obj.X_lb + obj.X_ub) / 2; % use mid point between lower and upper bound
                    case 'random'
                        rng('shuffle'); % shuffle the generator: https://de.mathworks.com/help/matlab/math/why-do-random-numbers-repeat-after-startup.html
                        obj.initialguess.X = obj.X_lb + (obj.X_ub-obj.X_lb).*rand(size(obj.X_lb));  % use random point from uniform distribution between lower and upper bound
                    case 'random_smart'
                        load('result');
                        state_names = fieldnames(obj.idx);
                        for iname = 1:length(state_names)
                            name = state_names{iname};
                            if ~isfield(result.problem.idx,name)
                                warning(['The variable ''',name,''' does not exist in initial guess file. The initial guess of this variable could not updated.']);
                                continue
                            elseif size(obj.idx.(name),1) ~= size(result.X(result.problem.idx.(name)),1);
                                warning(['The variable ''',name,''' has a different numbers of states. The initial guess of this variable could not updated.']);
                                continue
                            end
                            %> @todo What if symmetric and not symmetric
                            %> movement is loaded?
                            if size(obj.idx.(name),2) ~= size(result.X(result.problem.idx.(name)),2);
                                obj.initialguess.X(obj.idx.(name)) = resampleX(result.X(result.problem.idx.(name)),size(obj.idx.(name),2)); % function resampleX can be overloaded and adpated by the user
                            else
                                obj.initialguess.X(obj.idx.(name)) = result.X(result.problem.idx.(name));
                            end
                            obj.initialguess.info = result.info; %save information of result which is used as initial guess, e.g. for warmstart
                        end
                        rand_nums = 0.5 + (1.5-0.5).*rand(size(obj.initialguess.X));
                        obj.initialguess.X = obj.initialguess.X.*rand_nums;
                        % reinitialize mex for current model, as mex function was
                        % initialized with result model by loading result calling load(init,'result');
                        %
                        %> @todo: current hack without changing Model Class, is there
                        %> a better way to make sure that always the right model is
                        %> initialized?
                        if ~eq(obj.model,result.problem.model)
                            obj.model.saveobj; % saveobj sets obj.init = 0
                            obj.model = obj.model.loadobj(obj.model); % loadobj calls initMex when obj.init = 0
                        end
                    otherwise
                        try
                            load(init,'result');
                        catch
                            error('Make initialguess failed. Check whether result file exists and is in the search path.');
                        end
                        state_names = fieldnames(obj.idx);
                        for iname = 1:length(state_names)
                            name = state_names{iname};
                            if ~isfield(result.problem.idx,name)
                                warning(['The variable ''',name,''' does not exist in initial guess file. The initial guess of this variable could not updated.']);
                                continue
                            elseif size(obj.idx.(name),1) ~= size(result.X(result.problem.idx.(name)),1);
                                warning(['The variable ''',name,''' has a different numbers of states. The initial guess of this variable could not updated.']);
                                continue
                            end
                            %> @todo What if symmetric and not symmetric
                            %> movement is loaded?
                            if size(obj.idx.(name),2) ~= size(result.X(result.problem.idx.(name)),2);
                                obj.initialguess.X(obj.idx.(name)) = resampleX(result.X(result.problem.idx.(name)),size(obj.idx.(name),2)); % function resampleX can be overloaded and adpated by the user
                            else
                                obj.initialguess.X(obj.idx.(name)) = result.X(result.problem.idx.(name));
                            end
                            obj.initialguess.info = result.info; %save information of result which is used as initial guess, e.g. for warmstart
                        end
                        % reinitialize mex for current model, as mex function was
                        % initialized with result model by loading result calling load(init,'result');
                        %
                        %> @todo: current hack without changing Model Class, is there
                        %> a better way to make sure that always the right model is
                        %> initialized?
                        if ~eq(obj.model,result.problem.model)
                            obj.model.saveobj; % saveobj sets obj.init = 0
                            obj.model = obj.model.loadobj(obj.model); % loadobj calls initMex when obj.init = 0
                        end
                end
                obj.initialguess.info.type = init;

            end
            
        end
        
        %======================================================================
        %> @brief Function returning isPeriodic parameter
        %>
        %> @param   obj     Collocation class object
        %> @retval  isPer   Bool: see Collocation.isPeriodic
        %======================================================================
        function isPer =  get.isPeriodic(obj)
            iperiodicityConstraint = find(strcmp({obj.constraintTerms.name},'periodicityConstraint'), 1);
            if ~isempty(iperiodicityConstraint)
                isPer = 1;
            else
                isPer = 0;
            end
        end
        
        %======================================================================
        %> @brief Function returning isSymmetric parameter
        %>
        %> @param   obj     Collocation class object
        %> @retval  isSym   Bool: see Collocation.isSymmetric
        %======================================================================
        function isSym =  get.isSymmetric(obj)
            iperiodicityConstraint = find(strcmp({obj.constraintTerms.name},'periodicityConstraint'));
            if ~isempty(iperiodicityConstraint)
                isSym = obj.constraintTerms(iperiodicityConstraint).varargin{1};
            else
                isSym = 0;
            end
        end
        
        %======================================================================
        %> @brief Function returning nNodesDur parameter
        %>
        %> @param   obj     Collocation class object
        %> @retval  isPer   Bool: see Collocation.isPeriodic
        %======================================================================
        function nNodesDur =  get.nNodesDur(obj)
            if obj.isPeriodic
                nNodesDur = obj.nNodes+1;
            else
                nNodesDur = obj.nNodes;
            end
        end

    end
    
    methods(Static)
        % Declared function to plot multiple sim var tables
        plotMultSimVarTables(tables, style, plotStance)
        
        % Declared function to plot mean and SD of multiple sim var tables
        meanTables = plotMeanSimVarTable(tables, style)
    end
    
    
end
%======================================================================
%> @brief Function for resampling state vector for makeinitialguess
%>
%> @param   x         Matrix with initial guess
%> @param   N         Number of collocation nodes
%> @retval  xrep      Matrix with resampled initial guess
%======================================================================
function xrep = resampleX(x,N)
if size(x,2) == 1
    xrep = repmat(x,1,N);
else
    xq = (0:N-1)/N * (size(x,2)-1) + 1;
    xrep = interp1(x',xq)' ;
end
end
