% ======================================================================
%> @file @Gait2dc/Gait2dc.m
%> @brief Matlab class describing the Gait2dc model
%>
%> @author Ton, Eva, Anne, Marlies
%> @date October, 2017
% ======================================================================

% ======================================================================
%> @brief The class describes the Gait2dc model
%> @details
%> - The lower body is modeled in the sagittal plane.
%> - The model is configured using an excel file. This code relies on the
%>   structure of this file (hardcoded indices for each section in the file).
%>   Please don't change the structure. Thus, the order of the table entries
%>   should not be changed!!!
%> - In the list of the CPs (Gait2dc.CPs), additional CPs can be added.
% ======================================================================
classdef Gait2dc < Model
    
    properties  (SetAccess = protected)
        %> String: Filename of parameter file including path and extension
        excelfile
        %> Table: Defining size of feet for visualization. The x and y coordinates 
        %> are defined relative to the ankle joint in m. The length of the foot is equal to
        %> the length of the foot segment defined in the table Gait2dc.segments. The size of 
        %> the foot does not influence the ground contact.
        foot
    end
    
    properties (Constant) %@todo consider to compute heel-to-ankle distance (6/180*obj.bodyheight) and ankle height (0.039*obj.bodyheight) from bodyheight
        %> Double: Distance between heel and ankle joint for default body height in m
        HEELDEFAULTOFFSET = 0.06;
        %> Double: Distance between toe and CP toe in m
        TOEOFFSET = 0.05;
        %> Double: Distance between foot sole and ankle joint in m
        FOOTSOLEOFFSET = 0.0702;
    end
    
    properties (SetObservable, AbortSet, SetAccess = protected)
        %> Double: Bodyheight in m (default: 1.80)
        bodyheight = 1.80
        %> Double: Bodymass in kg (default: 75)
        bodymass = 75
        %> Double: Speed of simulated additional Treadmill
        speed_left = 0
        speed_right = 0

    end
    
    properties (SetObservable, AbortSet)
        %> Double array: Lambda values (default: [0.01 0.1]) (1x2)
        lambda = [0.01 0.1]
        %> Double: Slope of ground in deg (positive is uphill)
        slope
        %> Double matrix: Strain engergy terms to model apparel (20 x 7)
        strainEnergyTerms
    end
    
    properties(Dependent, SetAccess = protected, SetObservable, AbortSet)
        %> Double matrix: Entries of parameter excel file. It is linked to
        %> Gait2dc.bodyheight, Gait2dc.bodymass, Gait2dc.gravity, 
        %> Gait2dc.drag_coefficient, Gait2dc.wind_speed, Gait2dc.slope, 
        %> Gait2dc.dofs, Gait2dc.segments, Gait2dc.joints, Gait2dc.muscles 
        %> and Gait2dc.CPs.
        parameter
    end
    
    properties (Dependent, SetAccess = protected)
        %> Double array: Start indices of sections in Gait2dc.parameter
        idxParamSections
        %> Double array: Indices for x position of right contact points
        idxcxCPright
        %> Double array: Indices for x position of left contact points
        idxcxCPleft
    end
    
    %> @cond DO_NOT_DOCUMENT
    properties (Access = protected, Hidden) % use hidden states to store dependent variables to save computing time
        %> Hidden double matrix: Entries of parameter excel file
        hparameter
        %> Hidden double array: Start indices of sections in Gait2dc.parameter
        hidxParamSections
        %> Hidden double array: Indices for x position of right contact points
        hidxcxCPright
        %> Hidden double array: Indices for x position of left contact points
        hidxcxCPleft
    end
    %> @endcond
    
    methods
        
        %======================================================================
        %> @brief Default constructor setting default Gait2dc object
        %>
        %> @details
        %> Initializes the model and builds and initializes the mex function.
        %>
        %> The standard model can be called using:
        %> @code
        %> Gait2dc('gait2dc_par.xls')
        %> @endcode
        %> The excelfile must be in the matlab path.
        %>
        %> @param   varargin      Variable length input with:
        %>                          - excelfile     String: Excel file with path, name and extension.
        %>                                          If the excel is given, Gait2dc.scaleParameters will be called.
        %>                            OR
        %>                          - parameter     Matrix: Model parameters like they would be loaded
        %>                                          from the excel (see Gait2dc.parameter). If the matrix
        %>                                          is given and not an excel file, Gait2dc.scaleParameters
        %>                                          will not be called.
        %>                            AND
        %>                          - bodyheight    (optional) Double: Bodyheight in m
        %>                                          (If not given, default will be used and not the entry
        %>                                          from excel or parameter matrix. Use empty to skip.)
        %>                          - bodymass      (optional) Double: Bodymass in kg
        %>                                          (If not given, default will be used and not the entry
        %>                                          from excel or parameter matrix. Use empty to skip.)
        %> @retval  obj           Gait2dc class object
        %======================================================================
        function [obj] = Gait2dc(varargin)
            
            % Set listeners
            addlistener(obj,'lambda'           ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'gravity'          ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'drag_coefficient' ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'wind_speed'       ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'slope'            ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'joints'           ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'muscles'          ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'CPs'              ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'strainEnergyTerms','PostSet',@obj.update_mexParameter);
            addlistener(obj,'bodyheight'       ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'bodymass'         ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'dofs'             ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'segments'         ,'PostSet',@obj.update_mexParameter);
            
            addlistener(obj,'parameter'        ,'PostSet',@obj.update_idxParamSections);
            
            % Set input to properties
            if ischar(varargin{1})
                obj.excelfile = varargin{1};
                % Load parameter matrix
                obj.parameter = xlsread(obj.excelfile,'','','basic'); % Only numerical values
            elseif ismatrix(varargin{1})
                obj.excelfile = '';
                % Set parameter matrix
                obj.parameter = varargin{1};
            end
            
            % Read parameters from matrix to create tables
            obj.readParameter();

            % Initialize model (Set values and defaults)
            obj.initModel(varargin{2:end});

            % Scale parameters to scale the model to the height and weight
            if ~isempty(obj.excelfile)
                obj.scaleParameters;
            end

            % Check whether (up-to-date) compilation of model exists and if not (re-)compile
            Gait2dc.getMexFiles();
            
            % Initialize mex function
            obj.initMex;
            
        end
        
        %======================================================================
        %> @brief Function to set the speed of the ground, to model a
        %> treadmill
        %>
        %> @details
        %> This function allows you to model a treadmill. It can handle two
        %> types of inputs. If it is a double, the value will be set as the
        %> treadmill speed. If it is a struct, a split-belt treadmill can be
        %> modelled
        %>
        %> @param   obj     Gait2d_osim class object
        %> @param   speed   Double or Struct: Treadmill speed in m/s.
        %>                  If it is a double, the value is used for the left and right belt. 
        %>                  If is is a struct, speed.left and speed.right are used for the 
        %>                  left and right belt, respectively.  
        %======================================================================
        function setTreadmillSpeed(obj,speed)
            if isa(speed, 'double')
                % one speed in both belts
                obj.speed_left = speed;
                obj.speed_right = speed;
            elseif isa(speed, 'struct')
                % different speeds on left and right side
                obj.speed_left = speed.left;
                obj.speed_right = speed.right;
            else
                error('The variable speed should either be a double or a struct with the fields ''left'' and ''right''.')
            end
        end
        
        %======================================================================
        %> @brief Function to compute implicit differential equation for 2D musculoskeletal model
        %>
        %> @details
        %> This function calls the mex file of gait2dc.c:
        %> [f, dfdx, dfdxdot, dfdumus, dfdMextra]  = gait2dc('Dynamics',x,xdot,umus,Mextra);
        %>
        %> with the neural excitation umus and the extra torques Mextra.
        %>
        %> The dynamic residuals will be between fmin and fmax when inputs
        %> satisfy system dynamics: fmin <= f(x,dx/dt,umus,Mextra) <= fmax
        %>
        %> The last four outputs are optional and some computation time is saved if you do
        %> not request all of them.
        %>
        %> @param   obj     Gait2dc class object
        %> @param   x       Double array: State of the model (Gait2dc.nStates x 1)
        %> @param   xdot    Double array: State derivatives (Gait2dc.nStates x 1)
        %> @param   u       Double array: Controls of the model (Gait2dc.nControls x 1)
        
        %> @retval  f       Double array: Dynamic residuals (Gait2dc.nConstraints x 1)
        %> @retval	dfdx	(optional) Double matrix: Transpose of Jacobian matrix df/dx 		(Gait2dc.nStates x Gait2dc.nConstraints)
        %> @retval	dfdxdot	(optional) Double matrix: Transpose of Jacobian matrix df/dxdot 	(Gait2dc.nStates x Gait2dc.nConstraints)
        %> @retval	dfdu	(optional) Double matrix: Transpose of Jacobian matrix df/du 		(Gait2dc.nControls x Gait3d.nConstraints)
        %======================================================================
        function [f, dfdx, dfdxdot, dfdu] = getDynamics(obj,x,xdot,u)
            
            if obj.nTor == 0
                % Get neural excitation and extra moments
                % => Hard coded it here to speed up the function call! (Use Gait2dcTest.test_speedOfMex to test this)
                idxu = 1:length(u);
                umus = u;
                idxTorque = [];
                Mextra = zeros(6, 1); % for joints but not for orientation and position
            else
                % Get neural excitation
                idxu = obj.extractControl('u');
                umus = u(idxu);

                % Get extra moments
                Mextra = zeros(obj.nDofs, 1);
                idxTorque = obj.extractControl('torque');
                Mextra(obj.idxTorqueDof) = obj.mExtraScaleFactor * u(idxTorque); % Scale them by obj.mExtraScaleFactor and assume that order is consistent.
                Mextra = Mextra(~strcmp(obj.dofs.joint, 'ground_pelvis')); % mex takes only M for the joints, but not for orientation and position of the pelvis
            end

            % Apply treadmill speed
            xdot(obj.idxcxCPleft) = xdot(obj.idxcxCPleft) + obj.speed_left;
            xdot(obj.idxcxCPright) = xdot(obj.idxcxCPright) + obj.speed_right;

            % Get dynamics    
            if nargout > 3
                
                [f, dfdx, dfdxdot, dfdumus, dfdMextra] = gait2dc('Dynamics',x,xdot,umus,Mextra);
            
                % Get dfdu from dfdumus and dfdMextra
                dfdu = zeros(obj.nControls, obj.nConstraints);
                dfdu(idxu, :) = dfdumus;
                dfdu(idxTorque, :) = dfdMextra(obj.idxTorqueDof, :) * obj.mExtraScaleFactor;  % scaling has to be considered
                
            elseif nargout > 2
                [f, dfdx, dfdxdot] = gait2dc('Dynamics',x,xdot,umus,Mextra);
            elseif nargout > 1 
                [f, dfdx] = gait2dc('Dynamics',x,xdot,umus,Mextra);  
            else
                f = gait2dc('Dynamics',x,xdot,umus,Mextra);
            end
            

        end
        
        %======================================================================
        %> @brief Function returning the ground reaction forces for the system in state x
        %>
        %> @details
        %> This function calls the mex file of gait2dc.c:
        %> [grf, dgrfdx] = gait2dc('GRF', x);
        %>
        %> @param   obj     Gait2dc class object
        %> @param   x       Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  grf     Double array (12x1) containing:
        %>                   - right Fx, Fy, Fz (in bodyweight)
        %>                   - right Mx, My, Mz (in bodyweight*m)
        %>                   - left Fx, Fy, Fz (in bodyweight)
        %>                   - left Mx, My, Mz (in bodyweight*m)
        %> @retval	dgrfdx	(optional) Double matrix: Transpose of Jacobian matrix dgrf/dx (Gait2dc.nStates x 12)
        %======================================================================
        function [grf, dgrfdx] = getGRF(obj, x)
            
            if nargout > 1
                [grf_tmp, dgrfdx_tmp] = gait2dc('GRF', x);
                dgrfdx = zeros(length(x),12);
                dgrfdx(:,[1,2,6,7,8,12]) = dgrfdx_tmp;
            else
                [grf_tmp] = gait2dc('GRF', x);
            end
            grf = zeros(12,1);
            grf([1,2,6,7,8,12],1) = grf_tmp;
            
        end
        
        %======================================================================
        %> @brief Function returns position and orientation of all segments, for the system in position q
        %> 
        %> @details 
        %> This function can be used to plot the segments or for marker tracking.
        %> It calls the mex file of gait2d.c:
        %> FK = gait2dc('Fkin', x); 
        %>
        %> The order is the following for each segment:
        %>
        %> px', 'py', 'R11', 'R12', 'R21', 'R22'
        %> 
        %> p can be obtained by: 
        %> @code
        %> idxSegment = 1; % see Gait2dc.segments for the order
        %> p = FK((1:2) + (idxSegment-1)*6);
        %> @endcode
        %> R can be obtained by: 
        %> @code
        %> R = reshape(FK((3:6) + (idxSegment-1)*6), 2, 2)';
        %> @endcode
        %>
        %> @param   obj      Gait2dc class object
        %> @param   x        Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  FK		 Double array: position and orientation of all segments. Details see above. (Gait2dc.nSegments * 6)
        %> @retval	dFKdq	 (optional) Double matrix: Jacobian of FK with respect to q (Gait2dc.nSegments * 6 x Gait2dc.nDofs)
        %> @retval  dFKdotdq (optional) Double matrix: Jacobian of dFK/dt with respect to q (Gait2dc.nSegments * 6 x Gait2dc.nDofs)
        %======================================================================
        function [FK, dFKdq, dFKdotdq] = getFkin(obj,x)
            if nargout > 2
                [FK, dFKdq,dFKdotdq] = gait2dc('Fkin', x); 
            elseif nargout > 1
                [FK, dFKdq] = gait2dc('Fkin', x); 
            else
                 FK = gait2dc('Fkin', x);  
            end          
        end
        
        %======================================================================
        %> @brief Function to extract moment arms matrix from muscle table
        %>
        %> @details
        %> It extract the moment arm matrix which is defined by the model file
        %> (excel file) and saved in model.muscles. It is equal to the variable
        %> "MA" in gait2dc.c.
        %>
        %> It is not ensured whether the order in obj.muscles fit to to the
        %> order of the joints.
        %>
        %> @param   obj        Gait2dc class object
        %> @retval  momentarms Double matrix: Muscle moment arms in m (Gait2dc.nMus x Gait2dc.nJoints)
        %======================================================================
        function momentarms = getMomentArms(obj)

            variableNames = {'dRhip', 'dRknee', 'dRankle', 'dLhip', 'dLknee', 'dLank'};
            assert(length(variableNames)==obj.nJoints, ...
                'Muscle moment arms cannot be correctly extracted since hard coded column names do not fit to the number of joints.');

            momentarms = obj.muscles{:, variableNames};

        end

        %======================================================================
        %> @brief Function returns length of entire muscle-tendon-unit, for the system in state x
        %>
        %> @details
        %> This function replicated the computation done in gait2dc.c.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  L_MTU     Double array: MTU length (m) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function L_MTU = getLMTU(obj, x)

            % Get muscle moment arm for each muscle
            momentarms = obj.getMomentArms();

            % Get q without global translation and orientation
            dofNames = obj.dofs.Properties.RowNames(ismember(obj.dofs.joint, obj.joints.Properties.RowNames));
            qJoints = x(obj.extractState('q', dofNames), :);

            % Get length
            L_MTU = obj.muscles.L0 - momentarms * qJoints;

        end

        %======================================================================
        %> @brief Function returns length of contractile element, for the system in state x
        %>
        %> @details
        %> This function simply scales s using lceopt.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  L_CE      Double array: CE length (m) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function L_CE = getLCE(obj, x)

            % Get length
            L_CE = x(obj.extractState('s')) .* obj.muscles.lceopt;

        end

        %======================================================================
        %> @brief Function returns length of serial elastic element, for the system in state x
        %>
        %> @details
        %> This function uses the difference of the length of the entire muscle-tendon-unit
        %> and of the contracile element.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  L_SEE     Double array: SEE length (m) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function L_SEE = getLSEE(obj, x)

            % Get length
            L_SEE = obj.getLMTU(x) - obj.getLCE(x);

        end

        %======================================================================
        %> @brief Function returns length change of entire muscle-tendon-unit, for the system in state x
        %>
        %> @details
        %> This function does not ensure that the order of muscle moment
        %> arms fit to the order of DoFs in qdot in the states. This must
        %> be ensured by the user. The order of muscle moment arms is defined
        %> in the excel file defining the model.
        %>
        %> The same length is also computed in MusclePath() in gait2dc.c.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  Ldot_MTU  Double array: MTU length change (m/s) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function Ldot_MTU = getLdotMTU(obj, x)

            % Get muscle moment arm for each muscle
            momentarms = obj.getMomentArms();

            % Get qdot without global translation and orientation
            dofNames = obj.dofs.Properties.RowNames(ismember(obj.dofs.joint, obj.joints.Properties.RowNames));
            qdotJoints = x(obj.extractState('qdot', dofNames), :);

            % Compute length by applying the chain rule
            % d/dt(L_MTU) = -MA * d/dt(q) with MA = -dL_MTU/dq
            Ldot_MTU = -momentarms * qdotJoints;

        end

        %======================================================================
        %> @brief Function returns length change of contractile element, for the system in state x
        %>
        %> @details
        %> This function simply scales sdot using lceopt.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   xdot      Double array: State derivatives (Gait2dc.nStates x 1)
        %>
        %> @retval  Ldot_CE   Double array: CE length change (m/s) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function Ldot_CE = getLdotCE(obj, xdot)

            % Get length change
            Ldot_CE = xdot(obj.extractState('s')) .* obj.muscles.lceopt;

        end


        %======================================================================
        %> @brief Function returns length change of serial elastic element, for the system in state x
        %>
        %> @details
        %> This function computes the length change of the SEE as Ldot_MTU-Ldot_CE.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %> @param   xdot      Double array: State derivatives (Gait2dc.nStates x 1)
        %>
        %> @retval  Ldot_SEE  Double array: SEE length change (m/s) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function Ldot_SEE = getLdotSEE(obj, x, xdot)

            % Get length change of MTU
            Ldot_MTU = obj.getLdotMTU(x);

            % Get length change of CE
            Ldot_CE = obj.getLdotCE(xdot);

            % Compute length change of SEE
            Ldot_SEE = Ldot_MTU - Ldot_CE;

        end


        %======================================================================
        %> @brief Function returning muscle forces for the system in state x
        %>
        %> @details
        %> This function calls the mex file of gait2dc.c:
        %> [forces] = gait2dc('Muscleforces', x);
        %>
        %> @param   obj        Gait2dc class object
        %> @param   x          Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  forces     Double array: Muscle forces (in N) (Gait2dc.nMus x 1)
        %> @retval  dforcesdx  Double array: Transpose of Jacobian matrix d force / d x (Gait2dc.nStates x Gait2dc.nMus)
        %======================================================================
        function [forces,dforcesdx] = getMuscleforces(obj, x)
            if nargout == 1
                [forces] = gait2dc('Muscleforces', x);
            elseif nargout == 2
                [forces, dforcesdx] = gait2dc('Muscleforces', x);
            end
        end
        
        
        %======================================================================
        %> @brief Function returning muscle forces of CE for the system in state x
        %>
        %> @details
        %> This function calls the mex file of gait2dc.c:
        %> [forces] = gait2dc('MuscleCEforces', x);
        %>
        %> @param   obj             Gait2dc class object
        %> @param   x               Double array: State of the model (Gait2dc.nStates x 1)
        %> @param   xdot            Double array: State derivatives (Gait2dc.nStates x 1)
        %>
        %> @retval  forcesCE        Double array: Muscle forces of CE (in N) (Gait2dc.nMus x 1)
        %> @retval  dforcesCEdx     Double array: Transpose of Jacobian matrix d force / d x (Gait2dc.nStates x Gait2dc.nMus)
        %> @retval  dforcesCEdxdot  Double array: Transpose of Jacobian matrix d force / d xdot (Gait2dc.nStates x Gait2dc.nMus)
        %======================================================================
        function [forcesCE,dforcesCEdx,dforcesCEdxdot] = getMuscleCEforces(obj, x, xdot)
            if nargout == 1
                [forcesCE] = gait2dc('MuscleCEforces', x, xdot);
            elseif nargout > 1
                [forcesCE, dforcesCEdx, dforcesCEdxdot] = gait2dc('MuscleCEforces', x, xdot);
            end
        end


        %======================================================================
        %> @brief Function returns power generated by entire muscle-tendon-unit, for the system in state x
        %>
        %> @details
        %> This function calls Gait2dc:getLdotMTU() which does not ensure that
        %> the order of muscle moment arms fit to the order of DoFs in qdot in
        %> the states. This must be ensured by the user. The order of muscle moment
        %> arms is defined in the excel file defining the model.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %>
        %> @retval  powers	  Double array: MTU power output (W) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function powers = getMusclePower(obj, x)

            % Get force of SEE for each muscle
            forcesSEE = obj.getMuscleforces(x);

            % Get length change of MTU
            Ldot_MTU  = obj.getLdotMTU(x);

            % Compute power
            % P_MTU = -F_SEE .* d/dt(L_MTU)
            powers = -forcesSEE .* Ldot_MTU;

        end


        %======================================================================
        %> @brief Function returns power generated by muscle contractile elements, for the system in state x
        %>
        %> @details
        %> xdot must be the state derivatives, such that the muscle balance equations are satisfied.
        %> It is up to the user to ensure that the muscle contraction balance f(x,xdot)=0 when this
        %> function is used.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %> @param   xdot      Double array: State derivatives (Gait2dc.nStates x 1)
        %>
        %> @retval  powers	  Double array: CE power output (W) of each muscle (Gait2dc.nMus x 1)
        %> @retval  conDynRes Double array: Contraction dynamics residuals (Gait2dc.nMus x 1)
        %======================================================================
        function [powers, conDynRes] = getMuscleCEpower(obj, x, xdot)
            [powers, conDynRes] = gait2dc('MuscleCEpower', x, xdot);
        end
        
        %======================================================================
        %> @brief Function returns power generated by muscle serial elastic elements, for the system in state x
        %>
        %> @details
        %> xdot must be the state derivatives, such that the muscle balance equations are satisfied.
        %> It is up to the user to ensure that the muscle contraction balance f(x,xdot)=0 when this
        %> function is used.
        %>
        %> @param   obj       Gait2dc class object
        %> @param   x         Double array: State of the model (Gait2dc.nStates x 1)
        %> @param   xdot      Double array: State derivatives (Gait2dc.nStates x 1)
        %>
        %> @retval  powersSEE Double array: SEE power output (W) of each muscle (Gait2dc.nMus x 1)
        %======================================================================
        function powersSEE = getMuscleSEEpower(obj, x, xdot)

            % Get force of SEE for each muscle
            forcesSEE = obj.getMuscleforces(x);

            % Get length change of SEE
            Ldot_SEE = getLdotSEE(obj, x, xdot);

            % Compute power
            % P_SEE = - F_SEE .* L_SEEdot
            powersSEE = - forcesSEE .* Ldot_SEE;

        end


        %======================================================================
        %> @brief Function returns joint moments, for the system in state x
        %>
        %> @details
        %> This includes passive joint moments and muscle moments and extra moments.
        %>
        %> @param   obj     Gait2dc class object
        %> @param   x       Double array: State of the model (Gait2dc.nStates x 1)
        %> @param   u       (optional) Double array: Controls of the model which are only needed if you want 
        %>                  to apply extra moments (i.e. arm torques) (Gait2dc.nControls x 1)
        %>
        %> @retval  M       Double array: Moment for each DOF (Gait2dc.nDofs x 1)
        %> @retval  dMdx	(optional) Double matrix: Transpose of Jacobian matrix dM/dx (Gait2dc.nStates x Gait2dc.nDofs)
        %> @retval  dMdu    (optional) Double matrix: Transpose of Jacobian matrix dM/du (Gait2dc.nControls x Gait2dc.nDofs)
        %======================================================================
        function  [M, dMdx, dMdu] = getJointmoments(obj, x, u)
            
            % Get extra moments
            Mextra = zeros(obj.nDofs, 1);
            if nargin > 2 % Extra moments should be applied
                % Get extra moments
                idxTorque = obj.extractControl('torque');
                Mextra(obj.idxTorqueDof) = obj.mExtraScaleFactor * u(idxTorque); % Scale them by obj.mExtraScaleFactor and assume that order is consistent.
            end            
            
            % Call the mex
            idxMToUse = find(~strcmp(obj.dofs.joint, 'ground_pelvis')); % mex gives only M for the joints, but not for orientation and position of the pelvis
            if nargout > 1
                [M_tmp, dMdx_tmp] = gait2dc('Jointmoments', x);
                dMdx = zeros(obj.nStates, obj.nDofs); 
                dMdx(:, idxMToUse) = dMdx_tmp;
            else
                [M_tmp] = gait2dc('Jointmoments', x);
            end
            M = zeros(obj.nDofs,1);
            M(idxMToUse) = M_tmp;
            
            % Add extra moments to mex output (This could also be done in
            % the mex function like in the 3D model! See gait3d.c)
            M = M + Mextra;
            
            % Get dMdu
            if nargout > 2
                % Initialize with zeros
                dMdu = zeros(obj.nControls, obj.nDofs);
                
                if nargin > 2 % Extra moments were applied 
                    % Add identity matrix for dMdMextra
                    dMdMextra = eye(obj.nDofs);
                    dMdu(idxTorque, obj.idxTorqueDof) = dMdMextra(obj.idxTorqueDof, obj.idxTorqueDof) * obj.mExtraScaleFactor; % scaling has to be considered
                end
            end            
            
        end
        
        
        %======================================================================
        %> @brief Function to show model as stick figure
        %>
        %> @details
        %> @todo If the slope ~= 0, the ranges are not perfect yet. If ranges are given, they
        %> will not be kept constant. If no ranges are given, the stick figures are not perfectly centered.
        %>
        %> @param   obj            Gait2dc class object
        %> @param   x              Double matrice: State vector of model for n time points (Gait2dc.nStates x n)
        %> @param   range          (optional) Double matrice: Defining the range of the figure
        %>                         with [xmin, xmax; ymin, ymax]. (2 x 2)
        %>                         (default: [pelvisX-1, pelvisX+1; -0.2, 2])
        %> @param   plotFeet       (optional) Bool: If true, the undeformed feet are plotted independently
        %>                         from the CPs (default: 0)
        %> @param   plotGRF        (optional) Bool: If true, it plots arrows for the GRFs (default: 0)
        %> @param   plotCPs        (optional) Bool: If true, it plots spheres for the CPs (default: 0)
        %> @param   plotJointCOSYs (optional) Bool: If true, it plots the coordinate systems of the joints (default: 0)
        %======================================================================
        function showStick(obj,x, range, plotFeet, plotGRF, plotCPs, plotJointCOSYs)
            
            % get range of figures
            if nargin > 2 && size(range, 1) == 2 && size(range, 2) == 2
                xrange = range(1, :);
                yrange = range(2, :);
            else
                idxPelvisX = obj.extractState('q', 'pelvis_tx');
                xrange = [min(x(idxPelvisX,:))-1, max(x(idxPelvisX,:))+1];
                yrange = [-0.2, 2];
            end
            
            % plot the ground
            [xrange, groundY] = obj.applySlope(xrange, [0, 0]);
            yrange = yrange + sort(groundY);
            plot(xrange, groundY,'Color', [0.5, 0.5, 0.5]);
            hold on;
            
            % plot the stick figure
            for iTime = 1 : size(x,2)
                % get locations of joints and CPs
                [R, L, R1, L1] = gait2dc('Stick', x(:, iTime));
                % Apply slope
                [Rx, Ry] = obj.applySlope(R(:, 1), R(:, 2));
                [Lx, Ly] = obj.applySlope(L(:, 1), L(:, 2));
                [R1x, R1y] = obj.applySlope(R1(:, 1), R1(:, 2));
                [L1x, L1y] = obj.applySlope(L1(:, 1), L1(:, 2));
                % get also rotations if feet or joint COSYs should be plotted
                if (nargin > 3 && plotFeet) || (nargin > 6 && plotJointCOSYs)
                    assert(obj.slope == 0, 'Gait2dc:showStick(): The slope is not yet implemented for plotting the feet or the COSYs');

                    % Get position and orientation of all segments
                    FK = obj.getFkin(x(:, iTime));
                    FK = reshape(FK, 6, obj.nSegments);
                end
                if  nargin < 4 || ~plotFeet
                    % Plot the normal stick figure
                    plot(Lx,Ly,'b',  Rx(1:2),Ry(1:2),'k', Rx(2:end),Ry(2:end),'r', 'LineWidth',2);
                    % Plot the undeformed foot
                    plot(L1x,L1y,'k', R1x,R1y,'k');
                else
                    % Plot the normal stick figure without sole
                    plot(Lx(1:3),Ly(1:3),'b',  Rx(1:2),Ry(1:2),'k', Rx(2:4),Ry(2:4),'r', 'LineWidth',2);
                    
                    % Get foot position in global coordinate system
                    % depending on position and orientation of ankle
                    p_ankle_r = FK(1:2, 4); % 4: foot_r
                    R_ankle_r = reshape(FK(3:6, 4), 2, 2)';
                    heel_r = R_ankle_r * [obj.foot{'heel_r', 'x'}; obj.foot{'heel_r', 'y'}] + p_ankle_r;
                    toe_r  = R_ankle_r * [obj.foot{'toe_r', 'x'}; obj.foot{'toe_r', 'y'}]   + p_ankle_r;
                    p_ankle_l = FK(1:2, 7); % 7: foot_l
                    R_ankle_l = reshape(FK(3:6, 7), 2, 2)';
                    heel_l = R_ankle_l * [obj.foot{'heel_l', 'x'}; obj.foot{'heel_l', 'y'}] + p_ankle_l;
                    toe_l  = R_ankle_l * [obj.foot{'toe_l', 'x'}; obj.foot{'toe_l', 'y'}]   + p_ankle_l;
                    
                    % Compose coordinates for plotting
                    foot_x_r = [Rx(end); heel_r(1); toe_r(1); Rx(end)];
                    foot_y_r = [Ry(end); heel_r(2); toe_r(2); Ry(end)];
                    foot_x_l = [Lx(end); heel_l(1); toe_l(1); Lx(end)];
                    foot_y_l = [Ly(end); heel_l(2); toe_l(2); Ly(end)];
                    plot(foot_x_l, foot_y_l, 'b', foot_x_r, foot_y_r, 'r');
                    
                    % Plot deformed CPs
                    CPs_x_r = [heel_r(1); Rx(5:end-1); toe_r(1)];
                    CPs_y_r = [heel_r(2); Ry(5:end-1); toe_r(2)];
                    CPs_x_l = [heel_l(1); Lx(4:end-1); toe_l(1)];
                    CPs_y_l = [heel_l(2); Ly(4:end-1); toe_l(2)];
                    plot(CPs_x_l, CPs_y_l, 'b', CPs_x_r, CPs_y_r, 'r', 'LineWidth', 2);
                    
                    % Plot the undeformed CPs
                    CPs_un_x_r = [heel_r(1); R1x(2:end-1); toe_r(1)];
                    CPs_un_y_r = [heel_r(2); R1y(2:end-1); toe_r(2)];
                    CPs_un_x_l = [heel_l(1); L1x(2:end-1); toe_l(1)];
                    CPs_un_y_l = [heel_l(2); L1y(2:end-1); toe_l(2)];
                    plot(CPs_un_x_l, CPs_un_y_l, 'k', CPs_un_x_r, CPs_un_y_r, 'k');
                end
                
                % get the ground reaction forces and plot them
                if nargin > 4 && plotGRF
                    grf = obj.getGRF(x(:, iTime));
                    [CoP_r, CoP_l] = obj.getCoP(grf);
                    obj.plotgrf(grf(7:9), CoP_l, 'b');  % left side GRF vector
                    obj.plotgrf(grf(1:3), CoP_r, 'r');  % right side GRF vector
                end
                
                % plot the contact points as spheres
                if nargin > 5 && plotCPs
                    for i = 1 : obj.nCPs
                       curName = obj.CPs.Properties.RowNames{i};
                       xc = x(obj.extractState('xc', curName), iTime);
                       yc = x(obj.extractState('yc', curName), iTime);
                       [xc, yc] = obj.applySlope(xc, yc);
                       if endsWith(curName, '_r') % right leg
                           cpColor = [1, 0, 0]; % see colors of segments
                       elseif endsWith(curName, '_l') % left leg
                           cpColor = [0, 0, 1]; % see colors of segments
                       else % unknown segment
                           cpColor = [0, 0, 0]; % see colors of segments for other segment
                       end
                       markerSize = 40;
                       scatter(xc,yc,markerSize,'MarkerFaceColor',cpColor, 'MarkerEdgeColor', 'none');
                    end
                end

                % plot the axes of the COSYs
                if nargin > 6 && plotJointCOSYs
                    for iSeg = 1 : obj.nSegments
                        O = FK(1:2,iSeg);				    	% position of origin
                        R = reshape(FK(3:6,iSeg)',2,2)';        % rotation matrix

                        axislength = 0.1;
                        X = O + axislength*R(:,1);				% end of X axis
                        Y = O + axislength*R(:,2);				% end of Y axis
                        plot([O(1) X(1)],[O(2) X(2)],'m','LineWidth',2)	% show X axis
                        plot([O(1) Y(1)],[O(2) Y(2)],'c','LineWidth',2)	% show Y axis
                    end
                end
            end
            axis('equal');
            grid on; box on;
            xlabel('X');
            ylabel('Y');
            xlim(xrange);
            ylim(yrange);
            hold off;

        end

        %======================================================================
        %> @brief Function to visualize treadmill as moving scatter plot
        %>
        %> @details
        %>
        %> @param   obj            Gait2dc class object
        %> @param   xpoints        Double matrix: Defining the x - positions of
        %>                         scatter plot in current frame representing the treadmill
        %> @param   ypoints        Double matrix: Defining the y - positions of
        %>                         scatter plot in current frame representing the treadmill
        %> @param   xrange         Double matrix: Defining the range of the figure
        %>                         with [xmin, xmax]. (2 x 1)
        %> @param   yrange         Double matrix: Defining the range of the figure
        %>                         with [ymin, ymax]. (2 x 1)
        %> @param   durTrial       Double: duration of gait cycle
        %> @param   nFrames        Integer: number of Nodes
        %> @retval  x_points       Double matrix: Updated x_points
        %>                         accrording to treadmill speed
        %> @retval  y_points       Double matrix: Updated y_points
        %>                         accrording to treadmill speed
        %======================================================================
        function[x_points, y_points] = showTreadmill(obj, x_points, y_points, xrange, yrange, durTrial, nFrames)  
         
            if isprop(obj, 'slope') && obj.slope ~= 0
                error('Slope is not yet implemented for treadmill visualization')
            end
            
            % Number of points representing the treadmill belts
            nDots = 200;
            
            % Initialize points for treadmill scatter plot before first frame
            if isempty(x_points) && isempty(y_points) % first call
                x_points = rand(1, nDots) * (xrange(2) - xrange(1)) + xrange(1);
                y_points = rand(1, nDots) * yrange(1);
            else
                assert(length(x_points) == nDots, 'Gait2dc:showTreadmill', 'x_points has not the correct length of %i.', nDots);
                assert(length(y_points) == nDots, 'Gait2dc:showTreadmill', 'y_points has not the correct length of %i.', nDots);
            end
            
            % Right belt
            scatter(x_points(1,1:nDots/2), y_points(1,1:nDots/2), 10, 'filled', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha',0.5);
            hold on;
            % Left belt
            scatter(x_points(1,nDots/2+1:end), y_points(1,nDots/2+1:end),10, 'filled', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha',0.5);
 
            % Update the x-coordinates of the scatter points for right and
            % left belt
            x_points(1,1:nDots/2) = x_points(1,1:nDots/2) - (obj.speed_right * durTrial) / nFrames;
            x_points(1,nDots/2+1:end) = x_points(1,nDots/2+1:end) - (obj.speed_left *durTrial) / nFrames;
            % Check if the scatter points have disappeared on the left side
            if min(x_points) < xrange(1)
                % If they have, wrap around to the right side
                x_points(x_points < xrange(1)) = x_points(x_points < xrange(1)) + (xrange(2)-xrange(1));
            end
            % Check if the scatter points have disappeared on the right side
            if max(x_points) > xrange(2)
                % If they have, wrap around to the right side
                x_points(x_points > xrange(2)) = x_points(x_points > xrange(2)) - (xrange(2)-xrange(1));
            end
        end

        %======================================================================
        %> @brief Function to plot marker positions
        %>
        %> @details
        %> The motion in x and the matrix measured have to fit to each other.
        %> If you want to plot for example only a single node,
        %> markerTable.mean must be only the data of this single node.
        %>
        %> You can use the function like this:
        %> @code
        %> figure();
        %> x = result.X(result.problem.idx.states(:, 1:result.problem.nNodes)); % get states
        %> result.problem.model.showStick(x, [], 1, 0, 1); % plot stick figure
        %> markerTable = result.problem.objectiveTerms(1).varargin{1}.variables; % assuming trackMarker is your first objective
        %> markerMean = cell2mat(markerTable.mean')/1000; % extract and convert from mm to meter
        %> result.problem.model.showMarker(x, markerTable, markerMean); % plot marker data
        %> @endcode
        %>
        %> @todo The function is not yet supporting Gait2dc.slope ~= 0.
        %>
        %> @param   obj            Gait2dc class object
        %> @param   x              Double matrix: State vector of model for n time points (Gait2dc.nStates x n)
        %> @param   markerTable    Table: Variables table specifying markers with at least the columns:
        %>                         type, name, segment, position, direction to call Gait2dc.simuMarker.
        %> @param   measuredMean   (optional) Double matrix: Measured marker data in meter(!) which will be plotted as reference.
        %>                         The matrix must have the same number of time points as the state vector x.
        %>                         Further, the matrix columns must correspond to the rows of the markerTable in the same order
        %>                         to be able to match coordinates of the markers.  (n x height(markerTable)) (default: empty)
        %======================================================================
        function showMarker(obj,x, markerTable, measuredMean)

            assert(obj.slope == 0, 'Gait2dc:showMarker(): The slope is not yet implemented for plotting markers.');

            if nargin > 3 && ~isempty(measuredMean)
                assert(size(measuredMean, 1) == size(x, 2), ...
                    'Gait2dc:showMarker(): The matrix must have the same number of time points as the state vector x.');
                assert(size(measuredMean, 2) == height(markerTable), ...
                    'Gait2dc:showMarker(): The matrix columns must correspond to the rows of the markerTable.');
                plotMean = 1;
            else
                plotMean = 0;
            end

            % get rows with markers and match rows to each other based
            % on names
            markerNames = unique(markerTable.name,'stable');
            idxMarker = nan(2, numel(markerNames)); % 2D x marker
            for iName = 1 : numel(markerNames)
                idxMarkerName = find(strcmp(markerTable.name, markerNames{iName}));
                [dims, ~] = find(markerTable.direction(idxMarkerName, :));
                idxMarker(dims, iName) = idxMarkerName;
            end

            % get q to get later simulated data
            q = x(obj.extractState('q'), :);

            % for each time point
            hold on;
            for iNode = 1 : size(x,2)
                % plot simulated markers
                markerSim = obj.simuMarker(markerTable, q(:, iNode));
                scatter(markerSim(idxMarker(1, :)), markerSim(idxMarker(2, :)), 'ok');

                % plot measured markers
                if plotMean
                    scatter(measuredMean(iNode, idxMarker(1, :)), measuredMean(iNode, idxMarker(2, :)), 'xk');
                end
            end
            hold off;
           
        end
        
        
        % Declared function to simulate acceleration and gyroscope signals
        [s, ds_dq, ds_dqd, ds_dqdd] = simuAccGyro(obj, data, q, qd, qdd, idxSegment, idxAcc, idxGyro, dlocalAll, plocalAll)
        
        % Declared function to create graphical report fom a simulation result
        reportResult(obj, problem, result, resultfilename);
        
        % Declared function to create a movie fom a simulation result
        writeMovie(obj, problem, result, resultfilename);
        
        
        
        
        %> @cond DO_NOT_DOCUMENT  
        
        % Getter and setter function:
        %======================================================================
        %> @brief Function returning matrix with entries of parameter excel file
        %>
        %> @param   obj     Gait2dc class object
        %> @retval  p       Double matrix: Entries of parameter excel file
        %======================================================================
        function p =  get.parameter(obj)
            p = obj.hparameter;
        end
        
        %======================================================================
        %> @brief Function to set matrix with entries of parameter excel file
        %>
        %> @param   obj     Gait2dc class object
        %> @param   p       Double matrix: Entries of parameter excel file
        %======================================================================
        function set.parameter(obj, p)
            obj.hparameter = p;
        end
        
        %======================================================================
        %> @brief Function returning start indices of sections in Gait2dc.parameter
        %>
        %> @param   obj               Gait2dc class object
        %> @retval  idxParamSections  Double array: Start indices of sections in Gait2dc.parameter
        %======================================================================
        function idxParamSections = get.idxParamSections(obj)
            idxParamSections = obj.hidxParamSections;
        end
        %======================================================================
        %> @brief Function returning Model.idxcxCPright
        %>
        %> @param   obj     Model class object
        %> @retval  ifor    Double vector: Indices for x position of right contact points
        %======================================================================
        function idx = get.idxcxCPright(obj)
            idx = obj.hidxcxCPright;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxcxCPleft
        %>
        %> @param   obj     Model class object
        %> @retval  isid    Double vector: Indices for x position of left contact points
        %======================================================================
        function idx = get.idxcxCPleft(obj)
            idx = obj.hidxcxCPleft;
        end
        
        %> @endcond
        
    end
    
    methods(Access = protected)
        
        %======================================================================
        %> @brief Function to initialize model with default parameters
        %>
        %> @details
        %> Sets bodyheight and mass and creates default foot table.
        %>
        %> @param   obj         Gait2dc class object
        %> @param   varargin    Cell: Can contain
        %>                      - bodyheight  Double: Bodyheight in m
        %>                      - bodymass    Double: Bodymass in kg

        %======================================================================
        function initModel(obj, varargin)
            
            % Set bodyheight and bodymass
            if nargin > 1 && ~isempty(varargin{1})
                obj.bodyheight = varargin{1};
            end
            if nargin > 2 && ~isempty(varargin{2})
                obj.bodymass = varargin{2};
            end
            
            % Update the matrix by hand since it would not update if the
            % values were chosen even though there are other values in the matrix.
            obj.update_mexParameter(struct('Name', 'bodymass'));
            obj.update_mexParameter(struct('Name', 'bodyheight'));
            
            % Set default feet size
            obj.foot = array2table(nan(4, 2), 'VariableNames', {'x', 'y'});
            obj.foot.Properties.RowNames = {'heel_r'; 'toe_r'; 'heel_l'; 'toe_l'};
            obj.foot.x(1) = -Gait2dc.HEELDEFAULTOFFSET;
            obj.foot.y(1) = -Gait2dc.FOOTSOLEOFFSET;
            obj.foot.x(2) = obj.segments{'foot_r', 'length'} - Gait2dc.HEELDEFAULTOFFSET;
            obj.foot.y(2) = -Gait2dc.FOOTSOLEOFFSET;
            obj.foot.x(3) = -Gait2dc.HEELDEFAULTOFFSET;
            obj.foot.y(3) = -Gait2dc.FOOTSOLEOFFSET;
            obj.foot.x(4) = obj.segments{'foot_l', 'length'} - Gait2dc.HEELDEFAULTOFFSET;
            obj.foot.y(4) = -Gait2dc.FOOTSOLEOFFSET;
            
        end
        
        %======================================================================
        %> @brief Initialize model mex file
        %>
        %> @details
        %> This function calls the mex file of gait2dc.c
        %> [info] = gait2dc('Initialize', model)
        %>
        %> This initializes the model.  This is required before anything is done with the model.
        %>
        %> Input:
        %>	param		Matric containing model parameters loaded from the excel file after scaling.
        %>
        %> Output:
        %>	xneutral    Model state vector (Gait2dc.nStates x 1) for a neutral state where the model is in
        %>				free fall with feet not quite touching the ground, but otherwise close to
        %>				static equilibrium.
        %>
        %> @param   obj      Gait2dc class object
        %======================================================================
        function initMex(obj)
            
            obj.init = 0;
            
            % Initialize mex function
            xneutral = gait2dc('Initialize', obj.parameter);
            
            % Set default Labmda values in mexfunction
            if ~isempty(obj.lambda)
                gait2dc('Set','Lambda',obj.lambda);
            end
            
            % Write xneutral into table of states
            obj.hstates.xneutral = xneutral;
            
            obj.init = 1;
            
        end
        
        %======================================================================
        %> @brief Function performed to update the Gait2dc.parameter, the mex and the tables
        %>
        %> @details
        %> - Called with a PostSet event after set lambda, bodyheight, bodymass, 
        %>   gravity, drag_coefficient, wind_speed, slope, joints, muscles, CPs, 
        %>   strainEnergyTerms, dofs, segments.
        %> - Reinitializes mex with new parameter matrixs.
        %> - This function calls all the other update functions. If we
        %>   would also add a listener for them, we could not ensure that they
        %>   are called in the correct order!
        %>
        %> @param   obj     Gait2dc class object
        %> @param   src     meta.property: Object describing the source of the event
        %> @param   evnt    (optional) event.Propertyevent: Object describing the event
        %======================================================================
        function update_mexParameter(obj,src,evnt)
            
            % Change obj.parameter according to changes
            if ~strcmp(src.Name, 'lambda')
                parameter_tmp = obj.parameter; % use tmp variable such that the setter is not called multiple times.
                switch src.Name
                    case 'bodyheight'
                        parameter_tmp(obj.idxParamSections(1) +5, 2) = obj.bodyheight;
                    case 'bodymass'
                        parameter_tmp(obj.idxParamSections(1) +6, 2) = obj.bodymass;
                    case 'gravity'
                        parameter_tmp(obj.idxParamSections(1) +1, 2) = obj.gravity;
                    case 'drag_coefficient'
                        parameter_tmp(obj.idxParamSections(1) +2, 2) = obj.drag_coefficient;
                    case 'wind_speed'
                        parameter_tmp(obj.idxParamSections(1) +3, 2) = obj.wind_speed;
                    case 'slope'
                        assert(~(obj.speed_right ~= 0  && obj.speed_left ~= 0 && obj.slope ~= 0), ...
                            'Treadmill speed with slope is not yet implemented');
                        parameter_tmp(obj.idxParamSections(1) +4, 2) = obj.slope;
                    case 'joints'
                        obj.hnJoints = height(obj.joints);
                        for iRow = 1 : obj.nJoints
                            idx = obj.idxParamSections(obj.joints{iRow, 'iSection'}) + obj.joints{iRow, 'iRowInSection'};
                            parameter_tmp(idx, 4) = obj.joints{iRow, 'jointK2'};
                            parameter_tmp(idx, 5) = obj.joints{iRow, 'jointD'};
                            parameter_tmp(idx, 6) = obj.joints{iRow, 'jointB'};
                            parameter_tmp(idx, 7) = obj.joints{iRow, 'jointK1'};
                            parameter_tmp(idx, 8) = obj.joints{iRow, 'jointPhi0'};
                        end
                    case 'muscles'
                        obj.hnMus = height(obj.muscles); % Has to be updated before the states
                        for iRow = 1 : obj.nMus
                            idx = obj.idxParamSections(obj.muscles{iRow, 'iSection'}) + obj.muscles{iRow, 'iRowInSection'};
                            parameter_tmp(idx, 2) = obj.muscles{iRow, 'fmax'};
                            parameter_tmp(idx, 3) = obj.muscles{iRow, 'lceopt'};
                            parameter_tmp(idx, 4) = obj.muscles{iRow, 'width'};
                            parameter_tmp(idx, 5) = obj.muscles{iRow, 'peeSlack'};
                            parameter_tmp(idx, 6) = obj.muscles{iRow, 'seeSlack'};
                            parameter_tmp(idx, 7) = obj.muscles{iRow, 'L0'};
                            parameter_tmp(idx, 8) = obj.muscles{iRow, 'dRhip'};
                            parameter_tmp(idx, 9) = obj.muscles{iRow, 'dRknee'};
                            parameter_tmp(idx,10) = obj.muscles{iRow, 'dRankle'};
                            parameter_tmp(idx,11) = obj.muscles{iRow, 'dLhip'};
                            parameter_tmp(idx,12) = obj.muscles{iRow, 'dLknee'};
                            parameter_tmp(idx,13) = obj.muscles{iRow, 'dLank'};
                            parameter_tmp(idx,14) = obj.muscles{iRow, 'kPEE'};
                            parameter_tmp(idx,15) = obj.muscles{iRow, 'umax'};
                            parameter_tmp(idx,16) = obj.muscles{iRow, 'vmax'};
                            parameter_tmp(idx,17) = obj.muscles{iRow, 'tact'};
                            parameter_tmp(idx,18) = obj.muscles{iRow, 'tdeact'};
                            parameter_tmp(idx,19) = obj.muscles{iRow, 'gmax'};
                            parameter_tmp(idx,20) = obj.muscles{iRow, 'arel'};
                            parameter_tmp(idx,21) = obj.muscles{iRow, 'FT'};
                        end
                    case 'CPs'
                        obj.hnCPs = height(obj.CPs);
                        for iRow = 1 : obj.nCPs
                            idx = obj.idxParamSections(obj.CPs{iRow, 'iSection'}) + obj.CPs{iRow, 'iRowInSection'};
                            parameter_tmp(idx, 2) = obj.CPs{iRow, 'segmentID'};
                            parameter_tmp(idx, 3) = obj.CPs{iRow, 'x'};
                            parameter_tmp(idx, 4) = obj.CPs{iRow, 'y'};
                            parameter_tmp(idx, 5) = obj.CPs{iRow, 'k1'};
                            parameter_tmp(idx, 6) = obj.CPs{iRow, 'k2'};
                            parameter_tmp(idx, 7) = obj.CPs{iRow, 'a'};
                            parameter_tmp(idx, 8) = obj.CPs{iRow, 'c'};
                            parameter_tmp(idx, 9) = obj.CPs{iRow, 'b'};
                        end
                    case 'strainEnergyTerms'
                        parameter_tmp((1:20)+obj.idxParamSections(5)+1, 2:8) = obj.strainEnergyTerms;
                    case 'dofs'
                        obj.hnDofs = height(obj.dofs); % Has to be updated before the states
                        for iRow = 4 : obj.nDofs % range of first 3 dofs is not defined in excel file
                            idx = obj.idxParamSections(obj.joints{obj.dofs(iRow, 1).joint, 'iSection'}) + obj.joints{obj.dofs(iRow, 1).joint, 'iRowInSection'};
                            parameter_tmp(idx, 2) = obj.dofs{iRow, 'range'}(1)/pi*180; %Min
                            parameter_tmp(idx, 3) = obj.dofs{iRow, 'range'}(2)/pi*180; %Max
                        end
                    case 'segments'
                        obj.hnSegments = height(obj.segments);
                        for iRow = 1 : obj.nSegments
                            idx = obj.idxParamSections(obj.segments{iRow, 'iSection'}) + obj.segments{iRow, 'iRowInSection'};
                            parameter_tmp(idx, 2) = obj.segments{iRow, 'mass'};
                            parameter_tmp(idx, 3) = obj.segments{iRow, 'inertia'};
                            parameter_tmp(idx, 4) = obj.segments{iRow, 'mass_center'}(1); % CMx
                            parameter_tmp(idx, 5) = obj.segments{iRow, 'mass_center'}(2); % CMy
                            parameter_tmp(idx, 6) = obj.segments{iRow, 'length'};
                        end
                end
                obj.parameter = parameter_tmp;
            end
            
            % call other update functions
            switch src.Name
                case 'muscles'
                    obj.update_states;
                    obj.update_idxStates;
                    obj.update_idxSymmetry;
                    obj.update_idxForward;       % States have to be updated before
                    obj.update_idxUpward;        % States have to be updated before
                    obj.update_idxSideward;      % States have to be updated before
                    obj.update_idxForwardAll;    % States and idxForward have to be updated before
                    obj.update_idxSidewardAll;   % States and idxSideward have to be updated before
                    obj.update_controls;
                    obj.update_idxControls;
                    obj.update_idxcxCP;
                    obj.update_constraints;
                case 'CPs'
                    obj.update_states;
                    obj.update_idxStates;
                    obj.update_idxSymmetry;
                    obj.update_idxForward;       % States have to be updated before
                    obj.update_idxUpward;        % States have to be updated before
                    obj.update_idxSideward;      % States have to be updated before
                    obj.update_idxForwardAll;    % States and idxForward have to be updated before
                    obj.update_idxSidewardAll;   % States and idxSideward have to be updated before
                    obj.update_idxcxCP;
                    obj.update_constraints;
                case 'dofs'
                    obj.update_states;
                    obj.update_idxStates;
                    obj.update_idxSymmetry;
                    obj.update_idxForward;       % States have to be updated before
                    obj.update_idxUpward;        % States have to be updated before
                    obj.update_idxSideward;      % States have to be updated before
                    obj.update_idxForwardAll;    % States and idxForward have to be updated before
                    obj.update_idxSidewardAll;   % States and idxSideward have to be updated before
                    obj.update_idxTorqueDof;
                    obj.update_controls; % idxTorqueDof have to be updated before
                    obj.update_idxControls;
                    obj.update_idxcxCP;
                    obj.update_constraints;
            end
            
            % Reinitialize mex with new parameter matrix
            if ~isempty(obj.init) && obj.init
                obj.initMex;
            end
        end
        
        %======================================================================
        %> @brief Function performed to update start indices of sections in Gait2dc.parameter
        %>
        %> @details
        %> - Updates Gait2dc.idxParamSections
        %> - Called with a PostSet event after set parameter
        %> - Done in a similar way as in gait2dc.c
        %>
        %> @param   obj     Gait2dc class object
        %> @param   src     (optional) meta.property: Object describing the source of the event
        %> @param   evnt    (optional) event.Propertyevent: Object describing the event
        %======================================================================
        function update_idxParamSections(obj, src, evnt)
            
            % Get indices where a new section is starting
            MAXSECTIONS = 10;
            tmp_idxParamSections = [];
            nrows = 1;
            label = obj.parameter(nrows, 1);
            
            while label ~= 999
                if label >= 1 && label < MAXSECTIONS
                    tmp_idxParamSections = [tmp_idxParamSections, nrows];
                end
                nrows = nrows + 1;
                label = obj.parameter(nrows, 1);
            end
            
            % Set hidden variable
            obj.hidxParamSections = tmp_idxParamSections;
            
        end
        
        %======================================================================
        %> @brief Function defining the table Gait2dc.states
        %>
        %> @details
        %> Table lists type, name, xmin, xmax and xneutral
        %>
        %> @todo 
        %> Variables of contact points of neutral position shouldn't be zero.
        %>
        %> @param   obj     Gait2dc class object
        %======================================================================
        function update_states(obj)
            
            type = {};
            name = {};
            xmin = [];
            xmax = [];
            xneutral = [];
            
            % generalized coordinates and generalized velocities
            if ~isempty(obj.dofs)
                dofNames = obj.dofs.Properties.RowNames;
                type(1:obj.nDofs) = {'q'}; % generalized coordinates, dofs of the model
                name(1:obj.nDofs) = dofNames;
                xmin(1:obj.nDofs) = obj.dofs.range(:,1);
                xmax(1:obj.nDofs) = obj.dofs.range(:,2);
                xneutral(1:obj.nDofs) = 0;
                type(obj.nDofs+1:obj.nDofs*2) = {'qdot'}; %generalized velocities, velocitites of dofs
                name(obj.nDofs+1:obj.nDofs*2) = dofNames;
                xmin(obj.nDofs+1:obj.nDofs*2) = -30;
                xmax(obj.nDofs+1:obj.nDofs*2) = 30;
                xneutral(obj.nDofs+1:obj.nDofs*2) = 0;
                
                xneutral(strcmp(obj.dofs.Properties.RowNames,'pelvis_ty')) = 1;
            end
            
            % length of contractil element and muscle activation
            if ~isempty(obj.muscles)
                muscleNames = obj.muscles.Properties.RowNames;
                type(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = {'s'}; % length of contractil element
                name(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = muscleNames;
                xmin(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = -1;
                xmax(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = 3;
                xneutral(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = 2;
                type(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = {'a'}; % muscle activation
                name(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = muscleNames;
                xmin(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = 0;
                xmax(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = 5;
                xneutral(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2)= 0;
            end
            
            % global position for x and y (in m) and GRF in x and y (in BW)
            if ~isempty(obj.CPs)
                types = {'xc', 'yc', 'Fx', 'Fy'}; % xc,yc(in m), Fx,Fy (in BW)
                nxc = length(types); % each contact point has this many state variables
  
                contactNames = obj.CPs.Properties.RowNames;
                range_tmp = cell2mat(obj.CPs.range);
                cp_min = range_tmp(1:2:end,:);
                cp_max = range_tmp(2:2:end,:);
                
                for iNxc = 1 : nxc
                    type(obj.nDofs*2+obj.nMus*2+iNxc : nxc : obj.nDofs*2+obj.nMus*2+obj.nCPs*nxc) = types(iNxc);
                    name(obj.nDofs*2+obj.nMus*2+iNxc : nxc : obj.nDofs*2+obj.nMus*2+obj.nCPs*nxc) = contactNames;
                    xmin(obj.nDofs*2+obj.nMus*2+iNxc : nxc : obj.nDofs*2+obj.nMus*2+obj.nCPs*nxc) = cp_min(:,iNxc);
                    xmax(obj.nDofs*2+obj.nMus*2+iNxc : nxc : obj.nDofs*2+obj.nMus*2+obj.nCPs*nxc) = cp_max(:,iNxc);
                end
                xneutral(obj.nDofs*2+obj.nMus*2+1 : obj.nDofs*2+obj.nMus*2+obj.nCPs*nxc) = 0;
            end
            
            % create table
            obj.hstates = table(type',name',xmin',xmax',xneutral');
            obj.hstates.Properties.VariableNames = {'type','name','xmin','xmax','xneutral'};
            obj.hnStates= height(obj.states);
            
        end
        
        
        %======================================================================
        %> @brief Function defining Gait2dc.idxcxCPleft and Gait2dc.idxcxCPright
        %>
        %>
        %> @param   obj     Gait2dc class object
        %======================================================================
        function update_idxcxCP(obj)
            
            if ~isempty(obj.CPs)
                obj.hidxcxCPleft = obj.extractState('xc',obj.CPs.Properties.RowNames(obj.CPs.segmentID==1));
                obj.hidxcxCPright = obj.extractState('xc',obj.CPs.Properties.RowNames(obj.CPs.segmentID==0));
            end

        end
        
        %======================================================================
        %> @brief Function defining the table Gait2dc.controls
        %>
        %> @details
        %> Table lists type, name, xmin, xmax and xneutral
        %>
        %> @param   obj     Gait2dc class object
        %======================================================================
        function update_controls(obj)
            
            type = {};
            name = {};
            xmin = [];
            xmax = [];
            xneutral = [];
            
            %  neural excitation
            if ~isempty(obj.muscles)
                muscleNames = obj.muscles.Properties.RowNames;
                type(1:obj.nMus) = {'u'};
                name(1:obj.nMus) = muscleNames;
                xmin(1:obj.nMus) = 0;
                xmax(1:obj.nMus) = 5;
                xneutral(1:obj.nMus)= 0;
            end
            
            
            % additional torques
            if ~isempty(obj.idxTorqueDof)
                armNames = obj.dofs.Properties.RowNames(obj.idxTorqueDof);
                type(obj.nMus+1:obj.nMus+length(obj.idxTorqueDof)) = {'torque'};
                name(obj.nMus+1:obj.nMus+length(obj.idxTorqueDof)) = armNames;
                xmin(obj.nMus+1:obj.nMus+length(obj.idxTorqueDof)) = -5;
                xmax(obj.nMus+1:obj.nMus+length(obj.idxTorqueDof)) = 5;
                xneutral(obj.nMus+1:obj.nMus+length(obj.idxTorqueDof))= 0;
            end
            
            % create table
            obj.hcontrols = table(type',name',xmin',xmax',xneutral');
            obj.hcontrols.Properties.VariableNames = {'type','name','xmin','xmax','xneutral'};
            obj.hnControls = height(obj.controls);
            
        end
        
        %======================================================================
        %> @brief Function defining the table Gait2dc.constraints
        %>
        %> @details
        %> Table lists type, name, equation (string), fmin and fmax
        %>
        %> The constraints are composed of Gait2dc.nConstraints = 2* Gait2dc.nDofs + 2* Gait2dc.nMus + 4* Gait2dc.nCPs
        %>   - implicite differential equations: qdot-dq/dt = 0 (Gait2dc.nDofs x 1)
        %>   - equations of motion from Autolev (Gait2dc.nDofs x 1)
        %>   - muscle contraction dynamics (Gait2dc.nMus x 1)
        %>   - muscle activation dynamics: da/dt - (u-a)(c1*u + c2) = 0 (Gait2dc.nMus x 1)
        %>   - four equations for each contact point (4* Gait2dc.nCPs x 1)
        %>
        %> @param   obj          Gait2dc class object
        %======================================================================
        function update_constraints(obj)
            
            type = {};
            name = {};
            equation = {};
            fmin = [];
            fmax = [];
           
            if ~isempty(obj.dofs)
                % first nDofs elements of the implicit differential equation are: qdot-dq/dt = 0
                dofNames = obj.dofs.Properties.RowNames;
                type(1:obj.nDofs) = {'dofs'}; 
                name(1:obj.nDofs) = dofNames;
                equation(1:obj.nDofs) = {'qdot-dq/dt = 0'}; 
                fmin(1:obj.nDofs) = zeros(obj.nDofs,1);
                fmax(1:obj.nDofs) = zeros(obj.nDofs,1);
                % next nDofs elements of the IDE are the equations of motion from Autolev (the ZERO expressions)
                type(obj.nDofs+1:obj.nDofs*2) = {'dofs'};
                name(obj.nDofs+1:obj.nDofs*2) = dofNames;
                equation(obj.nDofs+1:obj.nDofs*2) = {'equations of motion from Autolev'}; 
                fmin(obj.nDofs+1:obj.nDofs*2) = zeros(obj.nDofs,1);
                fmax(obj.nDofs+1:obj.nDofs*2) = zeros(obj.nDofs,1);
            end
            
            if ~isempty(obj.muscles)
                % next nMus elements of the IDE are the muscle contraction dynamics
                muscleNames = obj.muscles.Properties.RowNames;
                type(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = {'muscles'};
                name(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = muscleNames;
                equation(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = {'contraction dynamics: f = Fsee - Fce - Fpee'}; 
                fmin(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = zeros(obj.nMus,1);
                fmax(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = zeros(obj.nMus,1);
                % next nMus elements of the IDE are the muscle activation dynamics: da/dt - (u-a)(c1*u + c2) = 0
                type(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = {'muscles'};
                name(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = muscleNames;
                equation(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = {'activation dynamics: da/dt - (u-a)(c1*u + c2) = 0'}; 
                fmin(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = zeros(obj.nMus,1);
                fmax(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = zeros(obj.nMus,1);
            end         
            
            if ~isempty(obj.CPs)
                % next nxc * nCps equations
                contactNames = obj.CPs.Properties.RowNames;
                nxc = 4; % each contact point has this many constraints
                
                type    (obj.nDofs*2+obj.nMus*2 + (1:obj.nCPs*nxc)) = {'CP'};
                fmin    (obj.nDofs*2+obj.nMus*2 + (1:obj.nCPs*nxc)) = zeros(nxc*obj.nCPs,1);
                fmax    (obj.nDofs*2+obj.nMus*2 + (1:obj.nCPs*nxc)) = zeros(nxc*obj.nCPs,1);
                
                name    (obj.nDofs*2+obj.nMus*2 + (1:nxc:obj.nCPs*nxc)) = contactNames;
                name    (obj.nDofs*2+obj.nMus*2 + (2:nxc:obj.nCPs*nxc)) = contactNames;
                name    (obj.nDofs*2+obj.nMus*2 + (3:nxc:obj.nCPs*nxc)) = contactNames;
                name    (obj.nDofs*2+obj.nMus*2 + (4:nxc:obj.nCPs*nxc)) = contactNames;
                equation(obj.nDofs*2+obj.nMus*2 + (1:nxc:obj.nCPs*nxc)) = {'equation 1'};
                equation(obj.nDofs*2+obj.nMus*2 + (2:nxc:obj.nCPs*nxc)) = {'equation 2'};
                equation(obj.nDofs*2+obj.nMus*2 + (3:nxc:obj.nCPs*nxc)) = {'equation 3'};
                equation(obj.nDofs*2+obj.nMus*2 + (4:nxc:obj.nCPs*nxc)) = {'equation 4'};
            end
            
            % create table
            obj.hconstraints = table(type',name',equation',fmin',fmax');
            obj.hconstraints.Properties.VariableNames = {'type','name','equation','fmin','fmax'};
            obj.hnConstraints = height(obj.constraints);
            
        end
        

        
        %======================================================================
        %> @brief Function defining Gait2dc.idxSymmetry
        %>
        %> @param   obj     Gait2dc class object
        %======================================================================
        function update_idxSymmetry(obj)
            % create the symmetry operator
            % x = xsign * x(xindex) will switch left and right in the state vector x
            % (assuming that the model walks in the +X direction)
            % u = usign * u(uindex) switches left and right in the controls u
            
            % generalized coordinates
            i_dofs = 1:obj.nDofs;
            s_dof = ones(1,obj.nDofs); % sign changes are not needed
            for iDof_r = 1:obj.nDofs
                % find right dofs
                if ismember(obj.dofs.Properties.RowNames{iDof_r}(end-1:end),'_r')
                    % find corresponding left dof
                    iDof_l = find(strcmp(obj.dofs.Properties.RowNames, [obj.dofs.Properties.RowNames{iDof_r}(1:end-2),'_l']));
                    % switch right and left
                    i_dofs(iDof_r) = iDof_l;
                    i_dofs(iDof_l) = iDof_r;
                end
            end
            
            % muscles
            i_muscles = 1:obj.nMus;
            s_muscles = ones(1, obj.nMus);
            for imuscles_r = 1:obj.nMus
                %find right muscles
                if ismember(obj.muscles.Properties.RowNames{imuscles_r}(end-1:end),'_r')
                    %find corresponding left muscle
                    imuscles_l = find(strcmp(obj.muscles.Properties.RowNames, [obj.muscles.Properties.RowNames{imuscles_r}(1:end-2),'_l']));
                    %switch right and left
                    i_muscles(imuscles_r) = imuscles_l;
                    i_muscles(imuscles_l) = imuscles_r;
                end
            end
            
            % contact points (CPs)
            i_contacts_Sing = 1:obj.nCPs;
            for icontacts_r = 1:obj.nCPs
                %find right CPs
                if ismember(obj.CPs.Properties.RowNames{icontacts_r}(end-1:end),'_r')
                    %find corresponding left muscle
                    icontacts_l = find(strcmp(obj.CPs.Properties.RowNames, [obj.CPs.Properties.RowNames{icontacts_r}(1:end-2),'_l']));
                    %switch right and left
                    i_contacts_Sing(icontacts_r) = icontacts_l;
                    i_contacts_Sing(icontacts_l) = icontacts_r;
                end
            end
            % Get the index with all indices (for all state variables of each CP)
            if ~isempty(obj.CPs)
                nxc = sum(strcmp(obj.states.name, obj.CPs.Properties.RowNames{1})); % each contact point has this many state variables
            else
                nxc = 0;
            end
            i_contacts = zeros(1, obj.nCPs * nxc);
            for icontacts = 1:obj.nCPs
                i_contacts(nxc*(icontacts-1) + (1:nxc)) = nxc*(i_contacts_Sing(icontacts)-1) + (1:nxc);
            end
            
            % sign changes are not needed (would be needed for z coordinate)
            s_contact = ones(1, nxc);
            
            % Get indices for dofs with torques
            i_tordofs = 1 : obj.nTor;
            s_tordofs = ones(1,obj.nTor); % sign changes are not needed
            names = obj.dofs.Properties.RowNames(obj.idxTorqueDof);
            for itordof_r = 1 : length(obj.idxTorqueDof)
                %find right torques
                if ismember(names{itordof_r}(end-1:end),'_r')
                    %find corresponding left torques
                    itordof_l = find(strcmp(names,[names{itordof_r}(1:end-2),'_l']));
                    %switch right and left
                    i_tordofs(itordof_r) = itordof_l;
                    i_tordofs(itordof_l) = itordof_r;
                end
            end

            % the full symmetry index list
            obj.hidxSymmetry.xindex = [i_dofs , obj.nDofs + i_dofs, ...		% gen coordinates and velocities
                2*obj.nDofs + i_muscles , ...					% muscle contraction states
                2*obj.nDofs + obj.nMus + i_muscles , ...	% muscle activation states
                2*obj.nDofs + 2*obj.nMus + i_contacts]';	% contact variables
            obj.hidxSymmetry.xsign = [s_dof s_dof ones(1,2*obj.nMus) repmat(s_contact,1,obj.nCPs)]';
            obj.hidxSymmetry.uindex = [i_muscles , obj.nMus + i_tordofs]';        % controls
            obj.hidxSymmetry.usign = [s_muscles, s_tordofs]';

            
        end
        
        %> @cond DO_NOT_DOCUMENT
        %======================================================================
        %> @brief Helper function to apply the slopes to the coordinates for showStick()
        %>
        %> @param   obj    Gait2dc class object
        %> @param   x      Double (vector/matrix): x coordinate before applying the slope
        %> @param   y      Double (vector/matrix): y coordinate before applying the slope
        %> @retval  xRot   Double (vector/matrix): x coordinate after applying the slope
        %> @retval  yRot   Double (vector/matrix): y coordinate after applying the slope
        %======================================================================
        function [xRot, yRot] = applySlope(obj, x, y)
            xRot = cosd(obj.slope)*x - sind(obj.slope)*y;
            yRot = sind(obj.slope)*x + cosd(obj.slope)*y;
        end
        
        %======================================================================
        %> @brief Helper function to plot GRFs in showStick()
        %>
        %> @details
        %> Plots a ground reaction force vector represented by 2D force and
        %> CoP. 
        %> 
        %> @param   obj    Gait2dc class object
        %> @param   F      Double array: Force (2 x 1)
        %> @param   CoP    Double array: Center of pressure in x, y coordinates (2 x 1)
        %> @param   color  String: Color (e.g. color = 'r')
        %======================================================================
        function plotgrf(obj,F,CoP,color)
            
            scale = 1.0;   % scale of the force vector visualization (meters per BW)
            x = [CoP(1) CoP(1)+scale*F(1)];
            y = [CoP(2) CoP(2)+scale*F(2)];
            [x, y] = obj.applySlope(x, y);
            plot(x,y,color,'LineWidth',2);
            
        end
        %> @endcond
        
        
        % Declared function to read excelfile with parameters
        readExcelfile(obj)
        
        % Declared function to scale the parameters defined in the excelfile
        scaleParameters(obj)
    end

    methods(Static)
       
        % Declared function tomake MEX functions
        getMexFiles(mexUpdate)
        
        %> @cond DO_NOT_DOCUMENT
        %======================================================================
        %> @brief Function to load saved object of Model model
        %>
        %> @details
        %> Checks if mex function is up to date and reinitializes mex function.
        %>
        %> @param   sobj    Saved Model class object
        %> @retval  obj     Model class object
        %======================================================================
        function obj = loadobj(sobj)
            obj = sobj;
            if ~isempty(obj.init) && obj.init == 0
                Gait2dc.getMexFiles;
                obj.initMex;
            end
        end
        %> @endcond
  
    end

end


