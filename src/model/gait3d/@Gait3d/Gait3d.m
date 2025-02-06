% ======================================================================
%> @file @Gait3d/Gait3d.m
%> @brief Matlab class describing the Gait3d model
%>
%> @author Ton, Eva, Anne, Marlies
%> @date October, 2017
% ======================================================================

%======================================================================
%> @brief The class describes the Gait3d model
%>
%> @details 
%> - The model is based on Sam Hamner's running model: 
%>   https://simtk.org/home/runningsim
%> - The code describes the basic Gait3d model
%> - The model is defined by the .osim file. Scaling must be done
%>   previously on the .osim file.
%> - Currently only the following models can be used
%>    - '3D Gait Model with Simple Arms'
%>    - '3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence'
%>    - 'gait24dof22musc'
%======================================================================
classdef Gait3d < Model
    
    properties (SetAccess = protected)
        %> Struct: Information (name, file, modified, sha256) on opensim model
        osim
        %> Double: Bodymass in kg (Is computed by summing up all body
        %> mass except talus. See gait3d.c)
        %> @todo Is is not nice that the talus is missing in the
        %> computation of the bodyweight.
        bodymass
        bodyheight
        %> Handle: MEX function of the model for dynamics etc.
        hdlMEX 
    end
    
    properties (SetObservable, AbortSet)
        %> Table: Markers (Containing also the CPs. Their names start with 'CP')
        markers
    end
    
    properties (Dependent, SetAccess = protected)
        %> Double: Number of markers (height of Gait3d.markers)
        nMarkers
    end
    

    
    methods
        
        %======================================================================
        %> @brief Constructor setting default Gait3d object
        %>
        %> @details
        %> Initializes the model and builds and initializes the mex function.
        %> The mex functions can be also build with other options if
        %> Gait3d.getMexFiles(osimName) is called before the constructor is called.
        %>
        %> Sets a new default for Gait3d.mExtraScaleFactor: 10 Nm
        %>
        %> The standard model can be called using:
        %> @code
        %> Gait3d('gait3d.osim')
        %> @endcode
        %> The osim file must be in the matlab path.
        %>
        %> @param   opensimfile     String: Opensim file with path, name and extension
        %> @retval  obj             Gait3d class object
        %======================================================================
        function [obj] = Gait3d(opensimfile)
            
            % Set listeners
            addlistener(obj,'gravity'           ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'drag_coefficient'  ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'wind_speed'        ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'joints'            ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'muscles'           ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'torques'           ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'CPs'               ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'dofs'              ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'segments'          ,'PostSet',@obj.update_mexParameter);
            
            obj.loadOsimFile(opensimfile); % read opensim file
            obj.loadMomentArms; % load moment arm file
            nameMEX = Gait3d.getMexFiles(obj.osim.name); % Check whether (up-to-date) compilation of model exists and if not (re-)compile 
            obj.hdlMEX = str2func(nameMEX); % Set the handle for the mex file depending on the model name
            obj.initModel; % initialize model with default parameters
            obj.initMex; % initialize mex function
            
            % Set a new default value for mExtraScaleFactor to use this
            % value from now on
            obj.mExtraScaleFactor = 10; % in Nm
        end

        %======================================================================
        %> @brief Function computing implicit differential equation for 3D musculoskeletal model
        %>
        %> @details 
        %> This function calls the mex file of gait3d.c:
        %> [f, dfdx, dfdxdot, dfdumus, dfdMextra]  = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);
        %>
        %> with the neural excitation umus and the extra torques Mextra.
        %>
        %> The dynamic residuals will be between fmin and fmax when inputs
        %> satisfy system dynamics : fmin <= f(x,dx/dt,umus,Mextra) <= fmax
        %>
        %> The last four outputs are optional and some computation time is saved if you do
        %> not request all of them.
        %>
        %> @param   obj     Gait3d class object
        %> @param   x       Double array: State of the model (Gait3d.nStates x 1)
        %> @param   xdot    Double array: State derivatives (Gait3d.nStates x 1)
        %> @param   u       Double array: Controls of the model (Gait3d.nControls x 1)
        %>
        %> @retval  f       Double array: Dynamic residuals (Gait3d.nConstraints x 1) 
        %> @retval	dfdx	(optional) Double matrix: Transpose of Jacobian matrix df/dx 		(Gait3d.nStates x Gait3d.nConstraints)
        %> @retval	dfdxdot	(optional) Double matrix: Transpose of Jacobian matrix df/dxdot 	(Gait3d.nStates x Gait3d.nConstraints)
        %> @retval	dfdu	(optional) Double matrix: Transpose of Jacobian matrix df/du 		(Gait3d.nControls x Gait3d.nConstraints)
        %====================================================================== 
        function [f, dfdx, dfdxdot, dfdu] = getDynamics(obj,x,xdot,u)
            
            % Get neural excitation
            idxu = obj.extractControl('u');
            umus = u(idxu); 
            
            % Get extra moments (e.g., arm torques)
            Mextra = zeros(obj.nDofs, 1);
            idxTorque = obj.extractControl('torque');
            Mextra(obj.idxTorqueDof) = obj.mExtraScaleFactor * u(idxTorque); % Scale them by obj.mExtraScaleFactor and assume that order is consistent. 
            
            % Get dynamics   
            if nargout > 3
            
                [f, dfdx, dfdxdot, dfdumus, dfdMextra] = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);
                
                % Get dfdu from dfdumus and dfdMextra
                dfdu = zeros(obj.nControls, obj.nConstraints);
                dfdu(idxu, :) = dfdumus;
                dfdu(idxTorque, :) = dfdMextra(obj.idxTorqueDof, :) * obj.mExtraScaleFactor;  % scaling has to be considered
                
            elseif nargout > 1
                [f, dfdx, dfdxdot] = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);
            else
                f = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);
            end
            
        end
        
        %======================================================================
        %> @brief Function returns the ground reaction forces for the system in state x
        %>
        %> @details 
        %> This function calls the mex file of gait3d.c:
        %> [grf, dgrfdx] = obj.hdlMEX('GRF', x);
        %>
        %> @param   obj     Gait3d class object
        %> @param   x       Double array: State of the model (Gait3d.nStates x 1)
        %>
        %> @retval  grf     Double array (12x1) containing:
        %>                   - right Fx, Fy, Fz (in bodyweight)
        %>                   - right Mx, My, Mz (in bodyweight*m)
        %>                   - left Fx, Fy, Fz (in bodyweight)
        %>                   - left Mx, My, Mz (in bodyweight*m)
        %> @retval	dgrfdx	(optional) Double matrix: Transpose of Jacobian matrix dgrf/dx (Gait3d.nStates x 12)
        %====================================================================== 
        function [grf, dgrfdx] = getGRF(obj,x)
            if nargout > 1
                [grf, dgrfdx] = obj.hdlMEX('GRF', x);
            else
                grf = obj.hdlMEX('GRF', x);
            end
        end
        
         
        %======================================================================
        %> @brief Function returns position and orientation of body segments, for the system in position q
        %> 
        %> @details 
        %> Needed for tracking of marker trajectories, and for 3D visualization of the model.
        %> This function calls the mex file of gait3d.c 
        %> [FK, dFKdq,dFKdotdq] = obj.hdlMEX('Fkin', q, qdot);
        %>
        %> The order the following for each segment:
        %>
        %> 'p1', 'p2', 'p3', 'R11', 'R12', 'R13', 'R21', 'R22', 'R23', 'R31', 'R32', 'R33'
        %> 
        %> p can be obtained by: 
        %> @code
        %> p = FK((idxCurSegment-2)*12 + (1:3));
        %> @endcode
        %> R can be obtained by: 
        %> @code
        %> R = reshape(FK((idxCurSegment-2)*12 + (4:12)), 3, 3)';
        %> @endcode
        %>
        %> @param   obj      Gait3d class object
        %> @param   q        Double array: Generalized coordinates (elements of states with type 'q') (Gait3d.nDofs x 1) 
        %> @param   qdot     (optional) Double array: Generalized velocities (elements of states with type 'qdot') (Gait3d.nDofs x 1)
        %>
        %> @retval  FK		 Double array: position (3) and orientation (3x3, stored row-wise) of nSegments-1. 
        %>                   Segments are the bodies in the .osim model, not including ground. ((Gait3d.nSegments-1) * 12)
        %> @retval	dFKdq	 (optional) Double matrix: Jacobian of FK with respect to q ((Gait3d.nSegments-1) * 12 x Gait3d.nDofs)
        %> @retval  dFKdotdq (optional) Double matrix: Jacobian of dFK/dt with respect to q ((Gait3d.nSegments-1) * 12 x Gait3d.nDofs)
        %======================================================================
        function [FK, dFKdq, dFKdotdq] = getFkin(obj,q, qdot)
            if nargout > 2
                [FK, dFKdq,dFKdotdq] = obj.hdlMEX('Fkin', q, qdot);
            elseif nargout > 1
                [FK, dFKdq] = obj.hdlMEX('Fkin', q);
            else
                FK = obj.hdlMEX('Fkin', q);
            end
            
        end
        
        %======================================================================
        %> @brief Function returns muscle forces, for the system in state x
        %>
        %> @details 
        %> This function calls the mex file of gait3d.c:
        %> [forces, lengths, momentarms] = obj.hdlMEX('Muscleforces', x);
        %>
        %> @param   obj         Gait3d class object
        %> @param   x           Double array: State of the model (Gait3d.nStates x 1)
        %>
        %> @retval  forces      Double array: Muscle forces (N) (Gait3d.nMus x 1)
        %> @retval	lengths     (optional) Double array: Muscle-tendon lengths (m) (Gait3d.nMus x 1)
        %> @retval	momentarms  (optional) Double matrix: Muscle moment arms (m) (Gait3d.nMus x Gait3d.nDofs)
        %====================================================================== 
        function [forces, lengths, momentarms] = getMuscleforces(obj,x)
            if nargout > 2
                [forces, lengths, momentarms] = obj.hdlMEX('Muscleforces', x);
            elseif nargout > 1
                [forces, lengths] = obj.hdlMEX('Muscleforces', x);
            else
                forces = obj.hdlMEX('Muscleforces', x);
            end
        end
        
        %======================================================================
        %> @brief Function returning muscle forces of CE for the system in state x
        %>
        %> @details
        %> This function calls the mex file of gait3d.c:
        %> [forces] = obj.hdlMEX('MuscleCEforces', x);
        %>
        %> @param   obj             Gait3d class object
        %> @param   x               Double array: State of the model (Gait3d.nStates x 1)
        %> @param   xdot            Double array: State derivatives (Gait3d.nStates x 1)
        %>
        %> @retval  forcesCE        Double array: Muscle forces of CE (in N) (Gait3d.nMus x 1)
        %> @retval  dforcesCEdx     Double array: Transpose of Jacobian matrix d force / d x (Gait3d.nStates x Gait3d.nMus)
        %> @retval  dforcesCEdxdot  Double array: Transpose of Jacobian matrix d force / d xdot (Gait3d.nStates x Gait3d.nMus)
        %======================================================================
        function [forcesCE,dforcesCEdx,dforcesCEdxdot] = getMuscleCEforces(obj, x, xdot)
            if nargout == 1
                [forcesCE] = obj.hdlMEX('MuscleCEforces', x, xdot);
            elseif nargout > 1
                [forcesCE, dforcesCEdx, dforcesCEdxdot] = obj.hdlMEX('MuscleCEforces', x, xdot);
%                 dforcesCEdx = nan(size(dforcesCEdx));      % They are not tested yet!
%                 dforcesCEdxdot = nan(size(dforcesCEdxdot));
            end
        end
        
        %======================================================================
        %> @brief Function returns power generated by muscle contractile elements, for the system in state x
        %> 
        %> @details 
        %> xdot must be the state derivatives, such that the muscle balance equations are satisfied
        %> It is up to the user to ensure that the muscle contraction balance f(x,xdot)=0 when this
        %> function is used.
        %>
        %> @param   obj       Gait3d class object
        %> @param   x         Double array: State of the model (Gait3d.nStates x 1)
        %> @param   xdot      Double array: State derivatives (Gait3d.nStates x 1)
        %>
        %> @retval  powers	  Double array: CE power output (W) of each muscle (Gait3d.nMus x 1) 
        %> @retval  conDynRes Double array: Contraction dynamics residuals (Gait3d.nMus x 1)
        %====================================================================== 
        function [powers, conDynRes] = getMuscleCEpower(obj, x, xdot)
            [powers, conDynRes] = obj.hdlMEX('MuscleCEpower', x, xdot);
        end
        
        %======================================================================
        %> @brief Function returns joint moments, for the system in state x
        %>
        %> @details
        %> This includes passive joint moments and muscle moments and extra moments (e.g. arms).
        %>
        %> @param   obj     Gait3d class object
        %> @param   x       Double array: State of the model (Gait3d.nStates x 1)
        %> @param   u       (optional) Double array: Controls of the model which are only needed if you want 
        %>                  to apply extra moments (i.e. arm torques) (Gait3d.nControls x 1)
        %>
        %> @retval  M       Double array: Moment/force for each DOF (Gait3d.nDofs x 1)
        %> @retval  dMdx	(optional) Double matrix: Transpose of Jacobian matrix dM/dx (Gait3d.nStates x Gait3d.nDofs)
        %> @retval  dMdu    (optional) Double matrix: Transpose of Jacobian matrix dM/du (Gait3d.nControls x Gait3d.nDofs)
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
            if nargout > 1
                [M, dMdx] = obj.hdlMEX('Jointmoments', x, Mextra);
            else
                M = obj.hdlMEX('Jointmoments', x, Mextra);
            end
            
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
        %> @brief Function returns passive joint moments, for the system in state x
        %>
        %> @details
        %> This function returns the passive joint moments due to ligaments
        %>
        %> @param   obj     Gait3d class object
        %> @param   x       Double array: State of the model (Gait3d.nStates x 1)
        %>
        %> @retval  M       Double array: Passive moment for each DOF (Gait3d.nDofs x 1)
        %> @retval  dMdx	(optional) Double matrix: Transpose of Jacobian matrix dM/dx (Gait3d.nStates x Gait3d.nDofs)
        %====================================================================== 
        function  [M, dMdx] = getPassiveJointmoments(obj, x)
            
            % Call the mex
            if nargout > 1
                [M, dMdx] = obj.hdlMEX('PassiveJointmoments', x);
            else
                M = obj.hdlMEX('PassiveJointmoments', x);
            end 
        end
        
        %======================================================================
        %> @brief Function returns length of contractile element, for the system in state x
        %>
        %> @details
        %> This function obtaines L_CE from s using the equations which are
        %> also implemented in the c code.
        %>
        %> @param   obj       Gait3d class object
        %> @param   x         Double array: State of the model (Gait3d.nStates x 1)
        %>
        %> @retval  L_CE      Double array: CE length (m) of each muscle (Gait3d.nMus x 1)
        %======================================================================
        function L_CE = getLCE(obj, x)

            % Get length state s
            s = x(obj.extractState('s'));

            % Compute L_CE
            pennatopt = obj.muscles.pennatopt;
            L_CE = nan(length(pennatopt), 1);
            for iMus = 1 : obj.nMus
                if pennatopt(iMus) < 0.01
                    L_CE(iMus) = s(iMus);
                else
                    b = sin(pennatopt(iMus));
                    L_CE(iMus) = sqrt(s(iMus)^2 + b^2);
                end
            end

            % Convert from lceopt to meter
            L_CE = L_CE .* obj.muscles.lceopt;

        end

        %======================================================================
        %> @brief Function returns length change of contractile element, for the system in state x
        %>
        %> @details
        %> This function obtaines Ldot_CE from sdot using the equations which are
        %> also implemented in the c code.
        %>
        %> @param   obj       Gait3d class object
        %> @param   xdot      Double array: State derivatives (Gait3d.nStates x 1)
        %>
        %> @retval  Ldot_CE   Double array: CE length change (m/s) of each muscle (Gait3d.nMus x 1)
        %======================================================================
        function Ldot_CE = getLdotCE(obj, xdot)

            % Get length state change
            sdot = xdot(obj.extractState('s'));

            % Compute Ldot_CE
            pennatopt = obj.muscles.pennatopt;
            Ldot_CE = nan(length(pennatopt), 1);
            for iMus = 1 : obj.nMus
                if pennatopt(iMus) < 0.01
                    Ldot_CE(iMus) = sdot(iMus);
                else
                    Ldot_CE(iMus) = sdot(iMus)*cos(pennatopt(iMus)); % not exactly the approximation which is used in the c code
                end
            end

            % Convert from lceopt/s to meter/s
            Ldot_CE = Ldot_CE .* obj.muscles.lceopt;

        end

        %======================================================================
        %> @brief Function to show model as stick figure
        %>
        %> @param   obj              Gait3d class object
        %> @param   x                Double matrice: State vector of model for n time points (Gait3d.nStates x n)
        %> @param   range            (optional) Double matrice: Defining the range of the figure
        %>                           with [xmin, xmax; ymin, ymax, zmin, zmax]. (3 x 2)
        %>                           (default: [pelvisX-1, pelvisX+1; -0.2, 2, pelvisZ-1, pelvisZ+1])
        %> @param   plotFeet         (optional) Integer: Defines if and how the feet are plotted (default: 0)
        %>                           - 0: Toe segments are not plotted, 
        %>                           - 1: Toe segments are plotted using double the center of mass in x direction,
        %>                           - 2: Feet segments are plotted using the predefined position of the CPs. In future, 
        %>                                we could think of using the simulated positions saved in the states.
        %> @param   plotGRF          (optional) Bool: If true, it plots arrows for the GRFs (default: 0)
        %> @param   plotCPs          (optional) Bool: If true, it plots spheres for the CPs (default: 0)
        %> @param   plotJointCOSYs   (optional) Bool: If true, it plots the coordinate systems of the joints (default: 0)
        %> @param   az               (optional) Double: Azimuth for the view of the figure (view(az, el)). (default: 90)
        %> @param   el               (optional) Double: Elevation for the view of the figure (view(az, el)). (default: 0)
        %======================================================================
        function showStick(obj,x, range, plotFeet, plotGRF, plotCPs, plotJointCOSYs, az, el)
            if nargin < 8
                az = 90;
            end
            if nargin < 9
                el = 0;
            end
            
            
            % make coordinates for a 1 cm radius sphere
            [xs,ys,zs] = sphere(10);
            xs = xs*0.01;
            ys = ys*0.01;
            zs = zs*0.01;
            
            % define lines or polygons to draw the segments
            if ~strcmp(obj.osim.name, 'gait14dof22musc') && ~strcmp(obj.osim.name, 'gait14dof22musc and Pelvis Rotation-Obliquity-Tilt Sequence') 
                % all our usual 3D models
                polygons = {   ...
                    [obj.joints.location(12,:);	obj.joints.location(7,:); obj.joints.location(2,:); obj.joints.location(12,:)]; ...	% pelvis
                    [0 0 0; obj.joints.t1_coefs(3,5) obj.joints.t2_coefs(3,5) 0.0]; ...	% right femur
                    [0 0 0;	obj.joints.location(4,:)]; ...										% right shank
                    []; ...																		% right talus
                    [0 0 0;	obj.joints.location(6,:)]; ...										% right calcaneus
                    []; ...																		% right toes
                    [0 0 0; obj.joints.t1_coefs(8,5) obj.joints.t2_coefs(8,5) 0.0]; ...	% left femur
                    [0 0 0;	obj.joints.location(9,:)]; ...										% left shank
                    []; ...																		% left talus
                    [0 0 0;	obj.joints.location(11,:)]; ...										% left calcaneus
                    []; ...																		% left toes
                    [0 0 0; obj.joints.location(13,:); obj.joints.location(17,:); 0 0 0]; ...	% torso
                    [0 0 0;	obj.joints.location(14,:)]; ...										% right humerus
                    [0 0 0;	obj.joints.location(15,:)]; ...										% right ulna
                    [0 0 0;	obj.joints.location(16,:)]; ...										% right radius
                    []; ...																		% right hand
                    [0 0 0;	obj.joints.location(18,:)]; ...										% left humerus
                    [0 0 0;	obj.joints.location(19,:)]; ...										% left ulna
                    [0 0 0;	obj.joints.location(20,:)]; ...										% left radius
                    []; ...																		% left hand
                    };
                colors_ploygons = {'k', 'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'b', 'k', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b'};
            else
                % model used for NeuIPS challenge 2019
                polygons = {   ...
                    [obj.joints.location(12,:);	obj.joints.location(7,:); obj.joints.location(2,:); obj.joints.location(12,:)]; ...	% pelvis
                    [0 0 0; obj.joints.t1_coefs(3,5) obj.joints.t2_coefs(3,5) 0.0]; ...	        % right femur
                    [0 0 0;	obj.joints.location(4,:)]; ...										% right shank
                    [0 0 0;	obj.joints.location(5,:)]; ...		     							% right talus
                    [0 0 0;	obj.joints.location(6,:)]; ...										% right calcaneus
                    []; ...																		% right toes
                    [0 0 0; obj.joints.t1_coefs(8,5) obj.joints.t2_coefs(8,5) 0.0]; ...         % left femur
                    [0 0 0;	obj.joints.location(9,:)]; ...										% left shank
                    [0 0 0;	obj.joints.location(10,:)]; ...			     						% left talus
                    [0 0 0;	obj.joints.location(11,:)]; ...										% left calcaneus
                    []; ...																		% left toes
                    [0 0 0; obj.joints.location(13,:)]; ...	                                    % torso
                    [] ...																		% head
                    };
                colors_ploygons = {'k', 'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'b', 'k', 'k'};
                
            end

            % get range of figures
            if nargin > 2 && size(range, 1) == 3 && size(range, 2) == 2
                xrange = range(1, :);
                yrange = range(2, :);
                zrange = range(3, :);
            else
                idxPelvisX = obj.extractState('q', 'pelvis_tx');
                idxPelvisZ = obj.extractState('q', 'pelvis_tz');
                xrange = [min(x(idxPelvisX,:))-1  , max(x(idxPelvisX,:))+1];
                yrange = [-0.2, 2];
                zrange = [min(x(idxPelvisZ,:))-1  , max(x(idxPelvisZ,:))+1];
            end
            
            % plot the ground
            fill3([zrange(1), zrange(2), zrange(2), zrange(1)], ... % sidewards
                  [xrange(1), xrange(1), xrange(2), xrange(2)], ... % forwards
                  [0, 0, 0, 0], ...                                 % upwards
                  [0.5, 0.5, 0.5], 'EdgeColor', [0.5, 0.5, 0.5]); 
            
            % plot the stick figure
            hold on;
            for iTime = 1 : size(x,2)
                % run the forward kinematics
                fk = obj.getFkin(x(1:obj.nDofs, iTime));
                fk = reshape(fk, 12, obj.nSegments-1);			% forward kinematics output from MEX function does not have ground
                
                for i=1:obj.nSegments-1
                    
                    O = fk(1:3,i);							% position of origin
                    R = reshape(fk(4:12,i)',3,3)';			% rotation matrix
                    
                    % draw the XYZ axes of the 18 segments in red, green, blue
                    if nargin > 6 && plotJointCOSYs
                        axislength = 0.1;
                        X = O + axislength*R(:,1);				% end of X axis
                        Y = O + axislength*R(:,2);				% end of Y axis
                        Z = O + axislength*R(:,3);				% end of Z axis
                        % lot ZXY as XYZ so Y will be up in the Matlab Window
                        plot3([O(3) X(3)],[O(1) X(1)],[O(2) X(2)],'r','LineWidth',2)	% show X axis
                        plot3([O(3) Y(3)],[O(1) Y(1)],[O(2) Y(2)],'g','LineWidth',2)	% show Y axis
                        plot3([O(3) Z(3)],[O(1) Z(1)],[O(2) Z(2)],'b','LineWidth',2)	% show Z axis
                    end
                    
                    % Get polygon coordinates
                    p = polygons{i}';
                    
                    % Adapt polygon if feet should be plotted
                    if nargin > 3 && plotFeet == 1
                        % Toe segments are plotted using double the center of mass in x direction
                        if ismember(obj.segments.Properties.RowNames{i+1}, {'toes_r', 'toes_l'})
                            p = [0 0 0;	2*obj.segments.mass_center(i+1, 1),0,0]';
                        end
       
                    elseif nargin > 3 && plotFeet == 2  
                        % Feet segments are plotted using the CPs   
                        if ismember(obj.segments.Properties.RowNames{i+1}, {'calcn_r', 'toes_r', 'calcn_l', 'toes_l'})
                            p = obj.CPs.position(obj.CPs.segmentindex==i+1, :)'; 
                            p = p(:, [1, 3, 4, 2]); %> @todo: Order of CPs should not be hard coded!
                        end
                    end
                    
                    % Draw each segment a black line or polygon
                    % ->transform polygon coordinates to global, and plot it
                    np = size(p,2);
                    if (np>0)
                        pg = repmat(O,1,np) + R*p;
                        fill3(pg(3,:),pg(1,:),pg(2,:),'k', 'edgecolor', colors_ploygons{i});
                    end
                    
                end
                
                % get the ground reaction forces and plot them
                if nargin > 4 && plotGRF
                    grf = obj.getGRF(x(:, iTime));
                    [CoP_r, CoP_l] = obj.getCoP(grf);
                    Gait3d.plotgrf(grf(1:3), CoP_r, 'r');  % right side GRF vector
                    Gait3d.plotgrf(grf(7:9), CoP_l, 'b');  % left side GRF vector
                end
                
                % plot the contact points as spheres
                if nargin > 5 && plotCPs
                    for i = 1 : obj.nCPs
                       curName = obj.CPs.Properties.RowNames{i};
                       xc = xs + x(obj.extractState('xc', curName), iTime);
                       yc = ys + x(obj.extractState('yc', curName), iTime);
                       zc = zs + x(obj.extractState('zc', curName), iTime);
                       if endsWith(curName, '_r') % right leg
                           cpColor = [1, 0, 0]; % see colors of polygons
                       elseif endsWith(curName, '_l') % left leg
                           cpColor = [0, 0, 1]; % see colors of polygons
                       else % unknown segment
                           cpColor = [0, 0, 0]; % see colors of polygons for other segment
                       end
                       surf(zc,xc,yc,'FaceColor',cpColor, 'EdgeColor', 'none');
                    end
                end               
                
            end
            
            % finish the plot
            grid on; box on;
            xlabel('Z');
            ylabel('X');
            zlabel('Y');
            
            axis equal
            ylim(xrange);
            xlim(zrange);
            zlim(yrange);
            view(az, el);
            hold off; 
            
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
        %> @param   obj            Gait3d class object
        %> @param   x              Double matrix: State vector of model for n time points (Gait3d.nStates x n)
        %> @param   markerTable    Table: Variables table specifying markers with at least the columns:
        %>                         type, name, segment, position, direction to call Gait3d.simuMarker.
        %> @param   measuredMean   (optional) Double matrix: Measured marker data in meter(!) which will be plotted as reference.
        %>                         The matrix must have the same number of time points as the state vector x.
        %>                         Further, the matrix columns must correspond to the rows of the markerTable in the same order
        %>                         to be able to match coordinates of the markers.  (n x height(markerTable)) (default: empty)
        %======================================================================
        function showMarker(obj,x, markerTable, measuredMean)

            if nargin > 3 && ~isempty(measuredMean)
                assert(size(measuredMean, 1) == size(x, 2), ...
                    'Gait3d:showMarker(): The matrix must have the same number of time points as the state vector x.');
                assert(size(measuredMean, 2) == height(markerTable), ...
                    'Gait3d:showMarker(): The matrix columns must correspond to the rows of the markerTable.');
                plotMean = 1;
            else
                plotMean = 0;
            end

            % get rows with markers and match rows to each other based
            % on names
            markerNames = unique(markerTable.name,'stable');
            idxMarker = nan(3, numel(markerNames)); % 3D x marker
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
                scatter3(markerSim(idxMarker(3, :)), markerSim(idxMarker(1, :)), markerSim(idxMarker(2, :)), 'ok'); % order as in showStick()

                % plot measured markers
                if plotMean
                    scatter3(measuredMean(iNode, idxMarker(3, :)), measuredMean(iNode, idxMarker(1, :)), measuredMean(iNode, idxMarker(2, :)), 'xk');
                end
            end
            hold off;

        end

        
        % Declared method to track acc and gyro data
        [s, ds_dq, ds_dqd, ds_dqdd] = simuAccGyro(obj, data, q, qd, qdd, idxSegment, idxAcc, idxGyro, dlocalAll, plocalAll);
        
        % Declared function to create graphical report fom a simulation result
        reportResult(obj, problem, result, resultfilename);

        %> @cond DO_NOT_DOCUMENT  
        % Getter and setter function:
        %======================================================================
        %> @brief Function returning Gait3d.nMarkers
        %>
        %> @param   obj     Gait3d class object
        %> @retval  n       Double: Number of markers
        %======================================================================
        function n = get.nMarkers(obj)
            n = size(obj.markers,1);
        end

        %> @endcond
        
    end
    
    
    methods (Access = protected)
        %======================================================================
        %> @brief Function to initialize model with default parameters
        %>
        %> @details 
        %> Initialize contact points and stiffness of dofs
        %>
        %> @param   obj     Gait3d class object
        %======================================================================
        function initModel(obj)
            % create contact points directly from markers
            % here a more complex shoe geometry could be added
            iCP = 1;
            for imarkers=1:obj.nMarkers
                name = obj.markers.Properties.RowNames{imarkers};
                if strcmp(name(1:2),'CP')
                    name = strrep(name, 'CP_R', ''); % Remove leading R and L. This is not nice, but needed for update_isymmetry
                    name = strrep(name, 'CP_L', '');
                    CPs_tmp(iCP,1).segmentindex = obj.markers.segmentindex(imarkers);
                    CPs_tmp(iCP,1).segment = obj.markers.segment(imarkers);
                    CPs_tmp(iCP,1).name = [name, obj.markers.segment{imarkers}(end-1:end)];% This ending is needed (shoudl be _r or _l)               
                    
                    % contact variables are Fx,Fy,Fz (in BW), xc,yc,zc (in m)
                    CPs_tmp(iCP,1).range = [-2 -0.5 -2 -1 -0.5 -1; ...
                                         2   3   2  5  1.0  5];
                    
                    CPs_tmp(iCP,1).position = obj.markers.position(imarkers,:);
                    iCP = iCP + 1;
                end
            end
            
            % Make table
            CPs_tmp = struct2table(CPs_tmp);
            % Set name as row name
            CPs_tmp.Properties.RowNames = CPs_tmp.name;
            CPs_tmp.name = [];
            % Set property
            obj.CPs = CPs_tmp;
                        
            % add default stiffness parameters for all DOFs except the pelvis rotation and translation
            idxNoJoint = strcmp(obj.dofs.joint, 'ground_pelvis');
            dofs_tmp = obj.dofs; % use tmp variable such that the setter is not called multiple times.
            dofs_tmp.stiffness_K1 = ones(obj.nDofs,1);  % linear passive joint stiffness, in N/m or Nm/rad
            dofs_tmp.stiffness_K1(idxNoJoint) = NaN;
            dofs_tmp.stiffness_K2 = 5000*ones(obj.nDofs,1); % quadratic stiffness outside ROM in N/m^2 or Nm/rad^2
            dofs_tmp.stiffness_K2(idxNoJoint) = NaN;
            dofs_tmp.damping_B = ones(obj.nDofs,1); % damping in Ns/m or Nms/rad
            dofs_tmp.damping_B(idxNoJoint) = NaN;
            obj.dofs = dofs_tmp;
            
            % Add percentage of fast twitch fibers to muscle table since it
            % is needed to compute metabolic cost.
            % (Anne used a book of Ton to get them. She used the table with
            % that listed the most of the muscles that we have. And she
            % assumed 50% for the other muscles.)
            muscleNames = obj.muscles.Properties.RowNames; % names of single muscles for mapping
            FT = nan(height(obj.muscles), 1); % create new column
            FT(ismember(muscleNames, {'glut_med1_r', 'glut_med1_l'})) = 0.5; 
            FT(ismember(muscleNames, {'glut_med2_r', 'glut_med2_l'})) = 0.5; 
            FT(ismember(muscleNames, {'glut_med3_r', 'glut_med3_l'})) = 0.5; 
            FT(ismember(muscleNames, {'glut_min1_r', 'glut_min1_l'})) = 0.5; 
            FT(ismember(muscleNames, {'glut_min2_r', 'glut_min2_l'})) = 0.5; 
            FT(ismember(muscleNames, {'glut_min3_r', 'glut_min3_l'})) = 0.5; 
            FT(ismember(muscleNames, {'semimem_r'  , 'semimem_l'  })) = 0.5; 
            FT(ismember(muscleNames, {'semiten_r'  , 'semiten_l'  })) = 0.5; 
            FT(ismember(muscleNames, {'bifemlh_r'  , 'bifemlh_l'  })) = 0.35; 
            FT(ismember(muscleNames, {'bifemsh_r'  , 'bifemsh_l'  })) = 0.35; 
            FT(ismember(muscleNames, {'sar_r'      , 'sar_l'      })) = 0.5; 
            FT(ismember(muscleNames, {'add_long_r' , 'add_long_l' })) = 0.35; 
            FT(ismember(muscleNames, {'add_brev_r' , 'add_brev_l' })) = 0.55; 
            FT(ismember(muscleNames, {'add_mag1_r' , 'add_mag1_l' })) = 0.45; 
            FT(ismember(muscleNames, {'add_mag2_r' , 'add_mag2_l' })) = 0.45; 
            FT(ismember(muscleNames, {'add_mag3_r' , 'add_mag3_l' })) = 0.45; 
            FT(ismember(muscleNames, {'tfl_r'      , 'tfl_l'      })) = 0.3; 
            FT(ismember(muscleNames, {'pect_r'     , 'pect_l'     })) = 0.55; 
            FT(ismember(muscleNames, {'grac_r'     , 'grac_l'     })) = 0.45; 
            FT(ismember(muscleNames, {'glut_max1_r', 'glut_max1_l'})) = 0.5; 
            FT(ismember(muscleNames, {'glut_max2_r', 'glut_max2_l'})) = 0.5; 
            FT(ismember(muscleNames, {'glut_max3_r', 'glut_max3_l'})) = 0.5; 
            FT(ismember(muscleNames, {'iliacus_r'  , 'iliacus_l'  })) = 0.5; 
            FT(ismember(muscleNames, {'psoas_r'    , 'psoas_l'    })) = 0.5; 
            FT(ismember(muscleNames, {'quad_fem_r' , 'quad_fem_l' })) = 0.5; 
            FT(ismember(muscleNames, {'gem_r'      , 'gem_l'      })) = 0.5; 
            FT(ismember(muscleNames, {'peri_r'     , 'peri_l'     })) = 0.5; 
            FT(ismember(muscleNames, {'rect_fem_r' , 'rect_fem_l' })) = 0.55; 
            FT(ismember(muscleNames, {'vas_med_r'  , 'vas_med_l'  })) = 0.5; 
            FT(ismember(muscleNames, {'vas_int_r'  , 'vas_int_l'  })) = 0.5; 
            FT(ismember(muscleNames, {'vas_lat_r'  , 'vas_lat_l'  })) = 0.55; 
            FT(ismember(muscleNames, {'med_gas_r'  , 'med_gas_l'  })) = 0.45; 
            FT(ismember(muscleNames, {'lat_gas_r'  , 'lat_gas_l'  })) = 0.45; 
            FT(ismember(muscleNames, {'soleus_r'   , 'soleus_l'   })) = 0.25; 
            FT(ismember(muscleNames, {'tib_post_r' , 'tib_post_l' })) = 0.45; 
            FT(ismember(muscleNames, {'flex_dig_r' , 'flex_dig_l' })) = 0.6; 
            FT(ismember(muscleNames, {'flex_hal_r' , 'flex_hal_l' })) = 0.5; 
            FT(ismember(muscleNames, {'tib_ant_r'  , 'tib_ant_l'  })) = 0.3; 
            FT(ismember(muscleNames, {'per_brev_r' , 'per_brev_l' })) = 0.55; 
            FT(ismember(muscleNames, {'per_long_r' , 'per_long_l' })) = 0.4; 
            FT(ismember(muscleNames, {'per_tert_r' , 'per_tert_l' })) = 0.65; 
            FT(ismember(muscleNames, {'ext_dig_r'  , 'ext_dig_l'  })) = 0.6; 
            FT(ismember(muscleNames, {'ext_hal_r'  , 'ext_hal_l'  })) = 0.5; 
            FT(ismember(muscleNames, {'ercspn_r'   , 'ercspn_l'   })) = 0.4335; 
            FT(ismember(muscleNames, {'intobl_r'   , 'intobl_l'   })) = 0.5;
            FT(ismember(muscleNames, {'extobl_r'   , 'extobl_l'   })) = 0.5;
            obj.muscles.FT = FT;
            
            if strcmp(obj.osim.name, 'gait14dof22musc') 
               warning('Gait3d:initModel', 'Percentage of fast twitch fibers is not defined for all muscles of this model. This will be needed to compute metabolic cost.');
            end

        end
        
        %======================================================================
        %> @brief Initialize model mex file
        %>
        %> @details 
        %> This function calls the mex file of gait3d.c 
        %> [info] = obj.hdlMEX('Initialize', model)
        %>
        %> This initializes the model.  This is required before anything is done with the model.
        %>
        %> Input:
        %>	model		Struct containing model parameters. 
        %>
        %> Output:
        %>	info		Struct with information about the size of the model, fields:
        %>   	Nx				Number of state variables: 2* Gait3d.nDofs + 2* Gait3d.nMuscles + 9* Gait3d.nCPs
        %>		Nf				Number of functions returned by the Dynamics function
        %>		fmin,fmax		Lower and upper bounds for dynamics residuals f
        %>      Bodyweight      Sum of weights of all segments. This is used to set bodymass
        %>
        %> @param   obj     Gait3d class object
        %======================================================================
        function initMex(obj)
            obj.init = 0;
 
            % make struct for mex input
            % (Extracting parameters explicitly makes more clear, what is
            % required by gait3d.c)
            % After changing this parameters initMex must be called again.
            mexModel = struct();
            mexModel.gravity = obj.gravity;
            mexModel.drag_coefficient = obj.drag_coefficient;
            mexModel.wind_speed = obj.wind_speed;
            mexModel.nDofs = obj.nDofs;
            mexModel.nSegments = obj.nSegments;
            mexModel.nMus = obj.nMus;
            mexModel.nCPs = obj.nCPs;
            for iDof = 1:obj.nDofs
                mexModel.dofs{iDof} = table2struct(obj.dofs(iDof,:));
                mexModel.dofs{iDof}.name = obj.dofs.Properties.RowNames{iDof};
            end
            for iSeg = 1:obj.nSegments
                mexModel.segments{iSeg} = table2struct(obj.segments(iSeg,:));
                mexModel.segments{iSeg}.name = obj.segments.Properties.RowNames{iSeg};
            end
            for iMus = 1:obj.nMus
                mexModel.muscles{iMus} = table2struct(obj.muscles(iMus,:));
                mexModel.muscles{iMus}.name = obj.muscles.Properties.RowNames{iMus};
            end
            for iJoin = 1:obj.nJoints
                mexModel.joints{iJoin} = table2struct(obj.joints(iJoin,:));
                mexModel.joints{iJoin}.name = obj.joints.Properties.RowNames{iJoin};
            end
            mexModel.CPs = {};
            for iCP = 1:obj.nCPs
                mexModel.CPs{iCP} = table2struct(obj.CPs(iCP,:));
                mexModel.CPs{iCP}.name = obj.CPs.Properties.RowNames{iCP};
            end
            
            % Call mex
            init_tmp = obj.hdlMEX('Initialize',mexModel);
                
            % Set bodymass
            obj.bodymass = init_tmp.Bodyweight / norm(obj.gravity);
            if isempty(obj.bodyheight)
                obj.bodyheight = 1.8; %from opensim
            end
            % Write fmin and fmax into table of constraints
            obj.hconstraints.fmin(1:length(init_tmp.fmin)) = init_tmp.fmin;
            obj.hconstraints.fmax(1:length(init_tmp.fmax)) = init_tmp.fmax;

            obj.init = 1;
        end
        
        %======================================================================
        %> @brief Function to load the content from the opensim file
        %>
        %> @details
        %> Calls readOsim() if there is no .mat file corresponding to the
        %> requested osimfile. Therefore it searches for a file called 
        %> strrep(osimfile, '.osim', '.mat') and checks the sha code to see 
        %> weather the .mat file is still up to date with the opensim file.
        %>
        %> @param   obj         Gait3d class object
        %> @param   osimfile    String: Filename of opensim file including path and extension
        %====================================================================== 
        function loadOsimFile(obj,osimfile)
            % Compute and store hash
            fileID = fopen(osimfile);
            fileID_read = fread(fileID);
            md = java.security.MessageDigest.getInstance('SHA-256');
            osim_sha256_new = typecast(md.digest(uint8(fileID_read))', 'uint8');
            matName = strrep(osimfile, '.osim', '.mat');
            if ~exist(matName,'file') 
                Gait3d.readOsim(osimfile);
            end
            load(matName,'model','osim_sha256')
            if (osim_sha256 ~= osim_sha256_new) & ~strcmp(computer,'MACA64') % OpenSim is not supported on Apple Silicon Mac
                model = Gait3d.readOsim(osimfile);
            end
            obj.dofs = model.dofs;
            obj.joints = model.joints;
            obj.segments = model.segments;
            obj.markers = model.markers;
            obj.muscles = model.muscles;
            obj.gravity = model.gravity;
            obj.osim = model.osim;
            obj.osim.sha256 = osim_sha256_new;
            
            % Set torque table
            % => This has to be done early to set passive moment arms
            % correctly
            switch obj.osim.name
                case {'3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence with Torques', ...
                    '3DGaitModelwithSimpleArmsandPelvisRotation-Obliquity-TiltSequencewithTorques', ...
                    'Running Model For Motions In All Directions With Torques', ...
                    'RunningModelForMotionsInAllDirectionsWithTorques'}
                    % torque driven model has torque actuators at each DOF (except for the global pelvis DOFs)
                    % (there is probably a better way to code this)
                    torques = struct();
                    dofNames = obj.dofs.Properties.RowNames(~ismember(obj.dofs.joint, 'ground_pelvis'));
                    for idof = 1:length(dofNames)
                        torques(idof,1).name = dofNames(idof);
                        torques(idof,1).dof = dofNames(idof);
                    end
                    torques = struct2table(torques);
                    torques.Properties.RowNames = torques.name;
                    torques.name = [];
                    obj.torques = torques;
                case {'3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence', ...
                        '3DGaitModelwithSimpleArmsandPelvisRotation-Obliquity-TiltSequence', ...
                        'Running Model For Motions In All Directions', ...
                        'RunningModelForMotionsInAllDirections', 'scaled_model_1dof_s01', 'Catelli_high_hip_flexion'}
                    % gait3d_pelvis213 model
                    % (there are only torque actuators at the arms)
                    armdof_names = {'arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_r','pro_sup_r','arm_flex_l','arm_add_l','arm_rot_l','elbow_flex_l','pro_sup_l'};
                    torques = struct();
                    for iarm = 1:length(armdof_names)
                        torques(iarm,1).name = armdof_names(iarm);
                        torques(iarm,1).dof = armdof_names(iarm);
                    end
                    torques = struct2table(torques);
                    torques.Properties.RowNames = torques.name;
                    torques.name = [];
                    obj.torques = torques;
                otherwise
                    error('Model name in obj.sim.name not recognized.')
            end
        end 
        
        %======================================================================
        %> @brief Function to load or compute moment arms
        %>
        %> @details
        %> Loads moment arms save in the file strrep(obj.osim.file,'.osim','_momentarms.mat').
        %> If there is not such a file, the moment arms will be computed and saved.
        %>
        %> We define here also the ranges of the dofs used to get the
        %> muscle moments and the ones used to get the passive moments:
        %>  - range_muscleMoment: They were manually tune to result in good
        %>    polynomials.
        %>  - range_passiveMoment: By default, they are equal to the ranges 
        %>    used for the muscle moments. But for the arms, we smaller the
        %>    range of the osim file by 2 degrees at each boundary 
        %>    (e.g. range_osim = [-90, 90], range_passiveMoment = [-88, 88])
        %>    
        %>
        %>
        %>
        %> @param   obj         Gait3d class object
        %====================================================================== 
        function loadMomentArms(obj)
            
            % define ranges which are used to compute muscle arms
            range_muscleMoment = table(NaN(height(obj.dofs), 2), 'VariableNames', {'range_muscleMoment'}, 'RowNames', obj.dofs.Properties.RowNames);
            if ~ismember(obj.osim.name, {'gait14dof22musc' , 'gait14dof22musc and Pelvis Rotation-Obliquity-Tilt Sequence', 'sIMUlate', ...
                    '3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence with Torques', ...
                    '3DGaitModelwithSimpleArmsandPelvisRotation-Obliquity-TiltSequencewithTorques', ...
                    'Running Model For Motions In All Directions With Torques', ...
                    'RunningModelForMotionsInAllDirectionsWithTorques'})
                % all our usual 3D models
                range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'} = [- 28, 78] / 180 * pi;
                range_muscleMoment{'hip_adduction_r' , 'range_muscleMoment'} = [- 13, 13] / 180 * pi;
                range_muscleMoment{'hip_rotation_r'  , 'range_muscleMoment'} = [-  3, 18] / 180 * pi;
                range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'} = [-118,  8] / 180 * pi;
                range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'} = [- 38, 38] / 180 * pi;
                range_muscleMoment{'subtalar_angle_r', 'range_muscleMoment'} = [- 13, 13] / 180 * pi;
                range_muscleMoment{'mtp_angle_r'     , 'range_muscleMoment'} = [-  8, 48] / 180 * pi;
                range_muscleMoment{'hip_flexion_l'   , 'range_muscleMoment'} = range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'};
                range_muscleMoment{'hip_adduction_l' , 'range_muscleMoment'} = range_muscleMoment{'hip_adduction_r' , 'range_muscleMoment'};
                range_muscleMoment{'hip_rotation_l'  , 'range_muscleMoment'} = range_muscleMoment{'hip_rotation_r'  , 'range_muscleMoment'};
                range_muscleMoment{'knee_angle_l'    , 'range_muscleMoment'} = range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'};
                range_muscleMoment{'ankle_angle_l'   , 'range_muscleMoment'} = range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'};
                range_muscleMoment{'subtalar_angle_l', 'range_muscleMoment'} = range_muscleMoment{'subtalar_angle_r', 'range_muscleMoment'};
                range_muscleMoment{'mtp_angle_l'     , 'range_muscleMoment'} = range_muscleMoment{'mtp_angle_r'     , 'range_muscleMoment'};
                range_muscleMoment{'lumbar_extension', 'range_muscleMoment'} = [- 38,  3] / 180 * pi;
                range_muscleMoment{'lumbar_bending'  , 'range_muscleMoment'} = [-  8,  8] / 180 * pi;
                range_muscleMoment{'lumbar_rotation' , 'range_muscleMoment'} = [- 18, 18] / 180 * pi;
            elseif strcmp(obj.osim.name, 'sIMUlate')
                % all our usual 3D models
                range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'} = [- 28, 78] / 180 * pi;
                range_muscleMoment{'hip_adduction_r' , 'range_muscleMoment'} = [- 13, 13] / 180 * pi;
                range_muscleMoment{'hip_rotation_r'  , 'range_muscleMoment'} = [-  3, 18] / 180 * pi;
                range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'} = [-118,  8] / 180 * pi;
                range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'} = [- 38, 38] / 180 * pi;
                range_muscleMoment{'subtalar_angle_r', 'range_muscleMoment'} = [- 13, 13] / 180 * pi;
                range_muscleMoment{'mtp_angle_r'     , 'range_muscleMoment'} = [-  8, 48] / 180 * pi;
                range_muscleMoment{'hip_flexion_l'   , 'range_muscleMoment'} = range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'};
                range_muscleMoment{'hip_adduction_l' , 'range_muscleMoment'} = range_muscleMoment{'hip_adduction_r' , 'range_muscleMoment'};
                range_muscleMoment{'hip_rotation_l'  , 'range_muscleMoment'} = range_muscleMoment{'hip_rotation_r'  , 'range_muscleMoment'};
                range_muscleMoment{'knee_angle_l'    , 'range_muscleMoment'} = range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'};
                range_muscleMoment{'ankle_angle_l'   , 'range_muscleMoment'} = range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'};
                range_muscleMoment{'subtalar_angle_l', 'range_muscleMoment'} = range_muscleMoment{'subtalar_angle_r', 'range_muscleMoment'};
                range_muscleMoment{'mtp_angle_l'     , 'range_muscleMoment'} = range_muscleMoment{'mtp_angle_r'     , 'range_muscleMoment'};
            elseif ismember(obj.osim.name, {'3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence with Torques', ...
                    '3DGaitModelwithSimpleArmsandPelvisRotation-Obliquity-TiltSequencewithTorques', ...
                    'Running Model For Motions In All Directions With Torques', ...
                    'RunningModelForMotionsInAllDirectionsWithTorques'})
                % Do nothing
            else
                % model used for NeuIPS challenge 2019
                range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'} = [- 28, 78] / 180 * pi;
                range_muscleMoment{'hip_adduction_r' , 'range_muscleMoment'} = [- 13, 13] / 180 * pi;
                range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'} = [-118,  8] / 180 * pi;
                range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'} = [- 38, 38] / 180 * pi;
                range_muscleMoment{'hip_flexion_l'   , 'range_muscleMoment'} = range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'};
                range_muscleMoment{'hip_adduction_l' , 'range_muscleMoment'} = range_muscleMoment{'hip_adduction_r' , 'range_muscleMoment'};
                range_muscleMoment{'knee_angle_l'    , 'range_muscleMoment'} = range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'};
                range_muscleMoment{'ankle_angle_l'   , 'range_muscleMoment'} = range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'};
                range_muscleMoment{'lumbar_extension', 'range_muscleMoment'} = [- 38,  3] / 180 * pi; % For lumbar_extension we use the one used for gait3d.osim. It's locked anyway.
            end
            
            % load or compute moment arms
            filename = strrep(obj.osim.file,'.osim','_momentarms.mat');
            filename = strsplit(filename,{'/','\'});
            filename = which(filename{end});
            if ~exist(filename, 'file')
                disp('Momentarm file could not be found. Momentarms will be computed.')
                obj.computeMomentArms(range_muscleMoment); % compute moment arms
            end
            load(which(filename),'momentarm_model','osim_sha256');
            if (osim_sha256 ~= obj.osim.sha256) & ~strcmp(computer,'MACA64') % OpenSim is not supported on Apple Silicon Mac
                b = input('The momentarm file is not up-to-date. Do you want new polynomials? (Y or N)','s');
                if (b(1) == 'Y') || (b(1) == 'y')
                    [momentarm_model] = obj.computeMomentArms(range_muscleMoment);
                end
            end
            
            % get polynomials
            lparam_count = zeros(obj.nMus,1);
            lparams = cell(obj.nMus,1); % the exponents for each coeficient - size is lparam_count x dof_count
            lcoefs = cell(obj.nMus,1);
            for imus = 1:obj.nMus
                % obj.muscles{imus}.dof_rom = momentarm_model{imus}.rom; % this is the range of motion for each dof where the momentarm model is valid
                lparam_count(imus) = momentarm_model{imus}.num_lparams; % the number of coeficients in the poly
                lparams{imus} = momentarm_model{imus}.lparams; % the exponents for each coeficient - size is lparam_count x dof_count
                lcoefs{imus} = momentarm_model{imus}.lcoef; % the array of coeficients - size is lparam_count x 1
                
            end
            obj.muscles.lparam_count = lparam_count;
            obj.muscles.lparams = lparams;
            obj.muscles.lcoefs = lcoefs;
            
            % add ranges which were used to compute the moment arms to the dofs table
            dofs_tmp = obj.dofs;
            dofs_tmp = [range_muscleMoment dofs_tmp];
            
            % Add ranges of motion in radians which are used to compute passive
            % moments (Are used in gait3d.c)
            if ~ismember(obj.osim.name, {'3D Gait Model with Simple Arms and Pelvis Rotation-Obliquity-Tilt Sequence with Torques', ...
                    '3DGaitModelwithSimpleArmsandPelvisRotation-Obliquity-TiltSequencewithTorques', ...
                    'Running Model For Motions In All Directions With Torques', ...
                    'RunningModelForMotionsInAllDirectionsWithTorques'})
                % By default equal to range_muscleMoment
                range_passiveMoment = range_muscleMoment; % default equal to range_muscleMoment
                range_passiveMoment.Properties.VariableNames{strcmp(range_passiveMoment.Properties.VariableNames, 'range_muscleMoment')} = 'range_passiveMoment';
                % For torque actuated dofs, use the range_osim and smaller the bounds
                amountSmallerBound = [-2, 2] / 180 * pi; % in radian
                for iTorqueDof = 1 : length(obj.idxTorqueDof)
                    range_passiveMoment{obj.idxTorqueDof(iTorqueDof) , 'range_passiveMoment'} = obj.dofs{obj.idxTorqueDof(iTorqueDof) , 'range_osim'} - amountSmallerBound;
                end
            else
                % The same as for the muscle model (see range_muscleModel for the standard model)
                range_passiveMoment = table(NaN(height(obj.dofs), 2), 'VariableNames', {'range_passiveMoment'}, 'RowNames', obj.dofs.Properties.RowNames);
                range_passiveMoment{'hip_flexion_r'   , 'range_passiveMoment'} = [- 28, 78] / 180 * pi;
                range_passiveMoment{'hip_adduction_r' , 'range_passiveMoment'} = [- 13, 13] / 180 * pi;
                range_passiveMoment{'hip_rotation_r'  , 'range_passiveMoment'} = [-  3, 18] / 180 * pi;
                range_passiveMoment{'knee_angle_r'    , 'range_passiveMoment'} = [-118,  8] / 180 * pi;
                range_passiveMoment{'ankle_angle_r'   , 'range_passiveMoment'} = [- 38, 38] / 180 * pi;
                range_passiveMoment{'subtalar_angle_r', 'range_passiveMoment'} = [- 13, 13] / 180 * pi;
                range_passiveMoment{'mtp_angle_r'     , 'range_passiveMoment'} = [-  8, 48] / 180 * pi;
                range_passiveMoment{'hip_flexion_l'   , 'range_passiveMoment'} = range_passiveMoment{'hip_flexion_r'   , 'range_passiveMoment'};
                range_passiveMoment{'hip_adduction_l' , 'range_passiveMoment'} = range_passiveMoment{'hip_adduction_r' , 'range_passiveMoment'};
                range_passiveMoment{'hip_rotation_l'  , 'range_passiveMoment'} = range_passiveMoment{'hip_rotation_r'  , 'range_passiveMoment'};
                range_passiveMoment{'knee_angle_l'    , 'range_passiveMoment'} = range_passiveMoment{'knee_angle_r'    , 'range_passiveMoment'};
                range_passiveMoment{'ankle_angle_l'   , 'range_passiveMoment'} = range_passiveMoment{'ankle_angle_r'   , 'range_passiveMoment'};
                range_passiveMoment{'subtalar_angle_l', 'range_passiveMoment'} = range_passiveMoment{'subtalar_angle_r', 'range_passiveMoment'};
                range_passiveMoment{'mtp_angle_l'     , 'range_passiveMoment'} = range_passiveMoment{'mtp_angle_r'     , 'range_passiveMoment'};
                range_passiveMoment{'lumbar_extension', 'range_passiveMoment'} = [- 38,  3] / 180 * pi;
                range_passiveMoment{'lumbar_bending'  , 'range_passiveMoment'} = [-  8,  8] / 180 * pi;
                range_passiveMoment{'lumbar_rotation' , 'range_passiveMoment'} = [- 18, 18] / 180 * pi;
                % For torque actuated dofs, use the range_osim and smaller the bounds
                amountSmallerBound = [-2, 2] / 180 * pi; % in radian
                for iDof = 24:33
                    range_passiveMoment{iDof, 'range_passiveMoment'} = obj.dofs{iDof, 'range_osim'} - amountSmallerBound;
                end
            end
            dofs_tmp = [range_passiveMoment, dofs_tmp];
            
            % reasign it to the dofs table
            obj.dofs = dofs_tmp;
        end
        
        %======================================================================
        %> @brief Function performed to update the parameter of the mex and the tables
        %>
        %> @details
        %> - Called with a PostSet event after set Gait3d.gravity, Gait3d.drag_coefficient, 
        %>   Gait3d.wind_speed, Gait3d.joints, Gait3d.muscles, Gait3d.CPs, Gait3d.dofs, Gait3d.segments.
        %> - Reinitializes mex with new parameters.
        %> - This function calls all the other update functions. If we
        %>   would also add a listener for them, we could not ensure that they
        %>   are called in the correct order!
        %>
        %> @param   obj     Gait3d class object
        %> @param   src     meta.property: Object describing the source of the event
        %> @param   evnt    (optional) event.Propertyevent: Object describing the event
        %======================================================================
        function update_mexParameter(obj,src,evnt)
            
            % call other update functions
            switch src.Name
                case 'joints'
                    obj.hnJoints = height(obj.joints);
                case 'segments'
                    obj.hnSegments = height(obj.segments);
                case 'muscles'
                    obj.hnMus = height(obj.muscles); % Has to be updated before the states
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
                    obj.update_constraints;
                case 'torques'
                    obj.hnTor = height(obj.torques);
                    obj.update_idxTorqueDof;     % Has to be updated first
                    obj.update_idxSymmetry;
                    obj.update_controls;
                    obj.update_idxControls;
                case 'CPs'
                    obj.hnCPs = height(obj.CPs);
                    obj.update_states;
                    obj.update_idxStates;
                    obj.update_idxSymmetry;
                    obj.update_idxForward;       % States have to be updated before
                    obj.update_idxUpward;        % States have to be updated before
                    obj.update_idxSideward;      % States have to be updated before
                    obj.update_idxForwardAll;    % States and idxForward have to be updated before
                    obj.update_idxSidewardAll;   % States and idxSideward have to be updated before
                    obj.update_constraints;
                case 'dofs'
                    obj.hnDofs = height(obj.dofs); % Has to be updated before the states
                    obj.update_states;
                    obj.update_idxStates;
                    obj.update_idxSymmetry;
                    obj.update_idxForward;       % States have to be updated before
                    obj.update_idxUpward;        % States have to be updated before
                    obj.update_idxSideward;      % States have to be updated before
                    obj.update_idxForwardAll;    % States and idxForward have to be updated before
                    obj.update_idxSidewardAll;   % States and idxSideward have to be updated before
                    obj.update_idxTorqueDof;  
                    obj.update_controls;         % idxTorqueDof have to be updated before
                    obj.update_idxControls;
                    obj.update_constraints;
            end
            
            % Reinitialize mex with new parameter matrix
            if ~isempty(obj.init) && obj.init
                obj.initMex;
            end
        end
        
        
        %======================================================================
        %> @brief Function defining the table Gait3d.states
        %>
        %> @details
        %> Table lists type, name, xmin, xmax and xneutral
        %>
        %> @todo 
        %> Variables of contact points of neutral position shouldn't be zero.
        %>
        %> @param   obj     Gait3d class object
        %======================================================================
        function update_states(obj)
            type = {};
            name = {};
            xmin = [];   
            xmax = [];  
            xneutral = [];
            if ~isempty(obj.dofs)
                dofNames = obj.dofs.Properties.RowNames;
                type(1:obj.nDofs) = {'q'}; % generalized coordinates, dofs of the model
                name(1:obj.nDofs) = dofNames;
                xmin(1:obj.nDofs) = obj.dofs.range_osim(:,1);
                xmax(1:obj.nDofs) = obj.dofs.range_osim(:,2);
                xneutral(1:obj.nDofs) = obj.dofs.neutral_position; % It is defined in gait3d.osim as default_value: Changing the default_value will change the passive joint moment!
                type(obj.nDofs+1:obj.nDofs*2) = {'qdot'}; %generalized velocities, velocitites of dofs
                name(obj.nDofs+1:obj.nDofs*2) = dofNames;
                xmin(obj.nDofs+1:obj.nDofs*2) = -30;
                xmax(obj.nDofs+1:obj.nDofs*2) = 30;
                xneutral(obj.nDofs+1:obj.nDofs*2) = 0;
            end
            
            
            if ~isempty(obj.muscles)
                muscleNames = obj.muscles.Properties.RowNames;
                type(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = {'s'}; % length of contractil element
                name(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = muscleNames;
                xmin(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = -1;
                xmax(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = 3;
                xneutral(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = 2; % muscle state variable s is set to 2.0 which should result in slack SEE and low muscle force
                type(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = {'a'};
                name(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = muscleNames;
                xmin(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = 0;
                xmax(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = 5;
                xneutral(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2)= 0; % muscle active state a is zero
            end
            
            if ~isempty(obj.CPs)
                types = {'Fx', 'Fy', 'Fz', 'xc', 'yc', 'zc'}; %Fx,Fy,Fz (in BW), xc,yc,zc (in m)
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
            
  
            obj.hstates = table(type',name',xmin',xmax',xneutral');
            obj.hstates.Properties.VariableNames = {'type','name','xmin','xmax','xneutral'};
            obj.hnStates= height(obj.states);

        end
        
        %======================================================================
        %> @brief Function defining the table Gait3d.controls
        %>
        %> @details
        %> Table lists type, name, xmin, xmax and xneutral
        %>
        %> @param   obj     Gait3d class object
        %======================================================================
        function update_controls(obj)
            type = {};
            name = {};
            xmin = [];   
            xmax = [];  
            xneutral = [];
            
            if ~isempty(obj.muscles)
                muscleNames = obj.muscles.Properties.RowNames;
                type(1:obj.nMus) = {'u'};
                name(1:obj.nMus) = muscleNames;
                xmin(1:obj.nMus) = 0;
                xmax(1:obj.nMus) = 5;
                xneutral(1:obj.nMus)= 0;
            end
            
            itors = obj.idxTorqueDof;
            if ~isempty(itors)
                torNames = obj.dofs.Properties.RowNames(itors);
                type(obj.nMus+1:obj.nMus+length(itors)) = {'torque'};
                name(obj.nMus+1:obj.nMus+length(itors)) = torNames;
                xmin(obj.nMus+1:obj.nMus+length(itors)) = -5;
                xmax(obj.nMus+1:obj.nMus+length(itors)) = 5;
                xneutral(obj.nMus+1:obj.nMus+length(itors))= 0;
            end
            obj.hcontrols = table(type',name',xmin',xmax',xneutral');
            obj.hcontrols.Properties.VariableNames = {'type','name','xmin','xmax','xneutral'};
            obj.hnControls = height(obj.controls);

        end	
        
        %======================================================================
        %> @brief Function defining the table Gait3d.constraints
        %>
        %> @details
        %> Table lists type, name, equation (string), fmin and fmax
        %>
        %> The constraints are composed of Gait3d.nConstraints = 2* Gait3d.nDofs + 2* Gait3d.nMus + 13* Gait3d.nCPs)
        %>   - implicite differential equations: qdot-dq/dt = 0 (Gait3d.nDofs x 1)
        %>   - equations of motion from Autolev (Gait3d.nDofs x 1)
        %>   - muscle contraction dynamics (Gait3d.nMus x 1)
        %>   - muscle activation dynamics: da/dt - rate * (u-a) = 0 (Gait3d.nMus x 1)
        %>   - six equality constraints for contact points (Gait3d.nCPs*6 x 1)
        %>
        %> @param   obj          Gait3d class object
        %======================================================================
        function update_constraints(obj)
            type = {};
            name = {};
            equation = {};
           
            if ~isempty(obj.dofs)
                % first nDofs elements of the implicit differential equation are: qdot-dq/dt = 0
                dofNames = obj.dofs.Properties.RowNames;
                type(1:obj.nDofs) = {'dofs'}; 
                name(1:obj.nDofs) = dofNames;
                equation(1:obj.nDofs) = {'qdot-dq/dt = 0'}; 
                % next nDofs elements of the IDE are the equations of motion from Autolev (the ZERO expressions)
                type(obj.nDofs+1:obj.nDofs*2) = {'dofs'};
                name(obj.nDofs+1:obj.nDofs*2) = dofNames;
                equation(obj.nDofs+1:obj.nDofs*2) = {'equations of motion from Autolev'}; 
            end
            
            if ~isempty(obj.muscles)
                % next nMus elements of the IDE are the muscle contraction dynamics
                muscleNames = obj.muscles.Properties.RowNames;
                type(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = {'muscles'};
                name(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = muscleNames;
                equation(obj.nDofs*2+1:obj.nDofs*2+obj.nMus) = {'contraction dynamics: f = Fsee - Fce - Fpee'}; 
                % next nMus elements of the IDE are the muscle activation dynamics: da/dt - rate * (u-a)= 0
                type(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = {'muscles'};
                name(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = muscleNames;
                equation(obj.nDofs*2+obj.nMus+1:obj.nDofs*2+obj.nMus*2) = {'activation dynamics: da/dt - rate * (u-a) = 0'}; 
            end         
            
            if ~isempty(obj.CPs)
                % next nxc * nCps equations
                contactNames = obj.CPs.Properties.RowNames;
                nEqu = 6; % each contact point has this many equality constraints
                type    (obj.nDofs*2+obj.nMus*2 + (1:obj.nCPs*nEqu)) = {'CP'};
                for iEqu = 1 : nEqu
                    name    (obj.nDofs*2+obj.nMus*2 + (iEqu:nEqu:obj.nCPs*nEqu)) = contactNames;
                    equation(obj.nDofs*2+obj.nMus*2 + (iEqu:nEqu:obj.nCPs*nEqu)) = {['equality ' num2str(iEqu)]};
                end
            end
            
            % create table
            fmin = NaN(size(type)); % Will be replaced in initMex() using output of mex
            fmax = NaN(size(type)); % Will be replaced in initMex() using output of mex
            obj.hconstraints = table(type',name',equation',fmin',fmax');
            obj.hconstraints.Properties.VariableNames = {'type','name','equation','fmin','fmax'};
            obj.hnConstraints = height(obj.constraints);
            
        end
        
        %======================================================================
        %> @brief Function defining Gait3d.idxSymmetry
        %>
        %> @param   obj     Gait3d class object
        %======================================================================
        function update_idxSymmetry(obj)
            % create the symmetry operator
            % x = xsign * x(xindex) will switch left and right in the state vector x
            % (assuming that the model walks in the +X direction)
            % u = usign * u(uindex) switches left and right in the controls u

            names_dof_signChange = {'pelvis_list','pelvis_rotation','pelvis_obliquity','pelvis_tz','lumbar_bending','lumbar_rotation'};

            % generalized coordinates
            i_dofs = 1:obj.nDofs;
            s_dof = ones(1,obj.nDofs);
            for iDof_r = 1:obj.nDofs
                %find right dofs
                if ismember(obj.dofs.Properties.RowNames{iDof_r}(end-1:end),'_r')
                    %find corresponding left dof
                    %[~,iDof_l] = Gait3d.getElementbyName(obj.dofs,[obj.dofs{iDof_r}.Properties.RowNames(1:end-2),'_l']);
                    iDof_l = find(strcmp(obj.dofs.Properties.RowNames,[obj.dofs.Properties.RowNames{iDof_r}(1:end-2),'_l']));
                    %switch right and left
                    i_dofs(iDof_r) = iDof_l;
                    i_dofs(iDof_l) = iDof_r;
                end
                % sign changes are needed in pelvis list, pelvis rotation, pelvis z,
                % lumbar bending and lumbar rotation
                if nnz(strcmp(obj.dofs.Properties.RowNames{iDof_r},names_dof_signChange))
                    s_dof(iDof_r) = -1;
                end
            end 
            
            % muscles: 1-43 are right leg, 44-86 are left leg, 87-92 are back muscles (RLRLRL)
            i_muscles = 1:obj.nMus;
            s_muscles = ones(1, obj.nMus);
            for imuscles_r = 1:obj.nMus
                %find right muscles
                if ismember(obj.muscles.Properties.RowNames{imuscles_r}(end-1:end),'_r')
                    %find corresponding left muscle
                    %[~,imuscles_l] = Gait3d.getElementbyName(obj.muscles,[obj.muscles{imuscles_r}.Properties.RowNames(1:end-2),'_l']);
                    imuscles_l = find(strcmp(obj.muscles.Properties.RowNames,[obj.muscles.Properties.RowNames{imuscles_r}(1:end-2),'_l']));
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
            nxc = 0;
            types = {''};
            if ~isempty(obj.CPs)
                types = obj.states{strcmp(obj.states.name, obj.CPs.Properties.RowNames{1}), 'type'}; % get types of CP variables
                nxc = length(types); % each contact point has this many state variables
            end
            i_contacts = zeros(1, obj.nCPs * nxc);
            for icontacts = 1:obj.nCPs
                i_contacts(nxc*(icontacts-1) + (1:nxc)) = nxc*(i_contacts_Sing(icontacts)-1) + (1:nxc);
            end
            % sign changes are needed in global Fz, global z coordinate of contact
            % point, and local z coordinate of foot point
            s_contact = ones(1, nxc);
            s_contact(strcmp(types, 'Fz')) = -1;
            s_contact(strcmp(types, 'zc')) = -1;
            
            % Get indices for dofs with torques
            i_tordofs = 1 : obj.nTor;
            s_tordofs = ones(1,obj.nTor);
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
                % sign changes are needed in pelvis list, pelvis rotation, pelvis z,
                % lumbar bending and lumbar rotation
                if nnz(strcmp(names{itordof_r},names_dof_signChange))
                    s_tordofs(itordof_r) = -1;
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

        % Declared function to compute moment arms using Opensim model
        [momentarms] = computeMomentArms(obj,range_muscleMoment,examine,name);
    end
    
    methods(Static, Access = protected)
        % Declared function to read Opensim file
        model = readOsim(osimfile)

        %======================================================================
        %> @brief Function used in showStick() to plot GRFs
        %>
        %> @details
        %> Plots a ground reaction force vector represented by 3D force and
        %> CoP. 
        %> 
        %> @param   F      Double array: Force (3 x 1)
        %> @param   CoP    Double array: Center of pressure in x, y, z coordinates (3 x 1)
        %> @param   color  String: Color (e.g. color = 'r')
        %======================================================================
        function plotgrf(F,CoP,color)
            
            scale = 1.0;   % scale of the force vector visualization (meters per BW)
            x = [CoP(1) CoP(1)+scale*F(1)];
            y = [CoP(2) CoP(2)+scale*F(2)];
            z = [CoP(3) CoP(3)+scale*F(3)];
            plot3(z,x,y,color,'LineWidth',2);
            
        end

        
    end
    
    methods(Static)
       
        % Declared function tomake MEX functions
        name_MEX = getMexFiles(osimName,clean, rebuild_contact, optimize)
        
        % Declared function to randomize model
        obj = createRandomPerson(obj, symmusrand, massratio,lengthratio)
  
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
                Gait3d.getMexFiles(obj.osim.name);
                obj.initMex;
            end
        end
        %> @endcond
        
    end    
end



