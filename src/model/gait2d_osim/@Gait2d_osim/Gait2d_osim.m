% ======================================================================
%> @file @Gait2d_osim/Gait2d_osim.m
%> @brief Matlab class describing the Gait2d_osim model. This is based on
%> the Gait3d model and OpenSim's gait10dof18musc
%>
%> @author Ton, Eva, Anne, Marlies, Markus
%> @date May, 2022
% ======================================================================

%======================================================================
%> @brief The class describes the Gait2d_osim model
%>
%> @details
%> - The model is exactly OpenSim's gait10dof18musc model
%> - The code is copied from gait3d entirely, with one dimension removed
%======================================================================
classdef Gait2d_osim < Model

    properties
        bodyheight; % This doesnt get set automatically, so no need to do protect it;
    end

    properties (SetAccess = protected)
        %> Struct: Information (name, file, modified, sha256) on opensim model
        osim
        %> Double: Bodymass in kg (Is computed by summing up all body
        %> mass except talus. See gait2d_osim.m)
        %> @todo Is is not nice that the talus is missing in the
        %> computation of the bodyweight.
        bodymass
        %> Handle: MEX function of the model for dynamics etc.
        hdlMEX
        %> hash: Hash to link the .tmp osim model
        hash
        %> Double: Speed of simulated additional Treadmill
        speed_left = 0
        speed_right = 0
    end

    properties (SetObservable, AbortSet)
        %> Table: Markers (Containing also the CPs. Their names start with 'CP')
        markers
        %> Copy of the segments Table from Model.m, to allow easy access for
        % changing the segments in the superclass
        segmentsTable = table();
    end

    properties (Dependent, SetAccess = protected)
        %> Double: Number of markers (height of Gait2d_osim.markers)
        nMarkers
        %> Double array: Indices for x position of right contact points
        idxcxCPright
        %> Double array: Indices for x position of left contact points
        idxcxCPleft
    end
    
    %> @cond DO_NOT_DOCUMENT
    properties (Access = protected, Hidden) % use hidden states to store dependent variables to save computing time
        %> Hidden double array: Indices for x position of right contact points
        hidxcxCPright
        %> Hidden double array: Indices for x position of left contact points
        hidxcxCPleft
    end
    %> @endcond

    methods

        %======================================================================
        %> @brief Constructor setting default Gait2d_osim object
        %>
        %> @details
        %> Initializes the model and builds and initializes the mex function.
        %> The mex functions can be also build with other options if
        %> Gait2d_osim.getMexFiles(osimName) is called before the constructor is called.
        %>
        %> Sets a new default for Gait2d_osim.mExtraScaleFactor: 10 Nm
        %>
        %> The standard model can be called using:
        %> @code
        %> Gait2d_osim('gait10dof18musc.osim')
        %> @endcode
        %> The osim file must be in the matlab path.
        %> @param   scale_factors   String, List or Float: .trc marker
        %scalefile, list of scalefactors: [nSegments] or [nSegments * 3],
        %(3d case); or float: bodyheight
        %> @param   bodyweight  Float: bodyweight in kg
        %> @param   opensimfile     String: Opensim file with path, name and extension

        %> @retval  obj             Gait2d_osim class object
        %======================================================================
        function [obj] = Gait2d_osim(opensimfile, scale_factors, bodymass, varargin)

            % Set listeners
            obj.add_listeners();
            % Add scaling:
            if nargin > 1
                if nargin < 3 
                    bodymass = 72.6;
                end
                % Scaling writes to a _tmp opensimfile, therefore
                [opensimfile, obj.hash] = Gait2d_osim.scaleOsim(opensimfile, scale_factors, bodymass, varargin{:});  
            end

            obj.loadOsimFile(opensimfile); % read opensim file
            obj.loadMomentArms; % load moment arm file
            nameMEX = Gait2d_osim.getMexFiles(obj.osim.name); % Check whether (up-to-date) compilation of model exists and if not (re-)compile
            obj.hdlMEX = str2func(nameMEX); % Set the handle for the mex file depending on the model name
            obj.initModel; % initialize model with default parameters
            obj.initMex; % initialize mex function

            % Set a new default value for mExtraScaleFactor to use this
            % value from now on
            obj.mExtraScaleFactor = 10; % in Nm

            % For 
            obj.segmentsTable = obj.segments(:,:);
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
        %> @brief Function computing implicit differential equation for 3D musculoskeletal model
        %>
        %> @details
        %> This function calls the mex file of gait2d_osim.c:
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
        %> @param   obj     Gait2d_osim class object
        %> @param   x       Double array: State of the model (Gait2d_osim.nStates x 1)
        %> @param   xdot    Double array: State derivatives (Gait2d_osim.nStates x 1)
        %> @param   u       Double array: Controls of the model (Gait2d_osim.nControls x 1)
        %>
        %> @retval  f       Double array: Dynamic residuals (Gait2d_osim.nConstraints x 1)
        %> @retval	dfdx	(optional) Double matrix: Transpose of Jacobian matrix df/dx 		(Gait2d_osim.nStates x Gait2d_osim.nConstraints)
        %> @retval	dfdxdot	(optional) Double matrix: Transpose of Jacobian matrix df/dxdot 	(Gait2d_osim.nStates x Gait2d_osim.nConstraints)
        %> @retval	dfdu	(optional) Double matrix: Transpose of Jacobian matrix df/du 		(Gait2d_osim.nControls x Gait2d_osim.nConstraints)
        %======================================================================
        function [f, dfdx, dfdxdot, dfdu] = getDynamics(obj,x,xdot,u)

            % Get neural excitation
            idxu = obj.extractControl('u');
            umus = u(idxu);

            % Get extra moments (arm torques)
            Mextra = zeros(obj.nDofs, 1);
            idxTorque = obj.extractControl('torque');
            Mextra(obj.hidxTorqueDof) = obj.mExtraScaleFactor * u(idxTorque); % Scale them by obj.mExtraScaleFactor and assume that order is consistent.

            % Apply treadmill speed
            xdot(obj.idxcxCPleft) = xdot(obj.idxcxCPleft) + obj.speed_left;
            xdot(obj.idxcxCPright) = xdot(obj.idxcxCPright) + obj.speed_right;

            % Get dynamics
            if nargout > 3

                [f, dfdx, dfdxdot, dfdumus, dfdMextra] = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);

                % Get dfdu from dfdumus and dfdMextra
                dfdu = zeros(obj.nControls, obj.nConstraints);
                dfdu(idxu, :) = dfdumus;
                dfdu(idxTorque, :) = dfdMextra(obj.hidxTorqueDof, :) * obj.mExtraScaleFactor;  % scaling has to be considered

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
        %> This function calls the mex file of gait2d_osim.c:
        %> [grf, dgrfdx] = obj.hdlMEX('GRF', x);
        %>
        %> @param   obj     Gait2d_osim class object
        %> @param   x       Double array: State of the model (Gait2d_osim.nStates x 1)
        %>
        %> @retval  grf     Double array (12x1) containing:
        %>                   - right Fx, Fy, Fz (in bodyweight)
        %>                   - right Mx, My, Mz (in bodyweight*m)
        %>                   - left Fx, Fy, Fz (in bodyweight)
        %>                   - left Mx, My, Mz (in bodyweight*m)
        %> @retval	dgrfdx	(optional) Double matrix: Transpose of Jacobian matrix dgrf/dx (Gait2d_osim.nStates x 12)
        %======================================================================
        function [grf, dgrfdx] = getGRF(obj,x)
            if nargout > 1
                [grf_tmp, dgrfdx_tmp] = obj.hdlMEX('GRF', x);
                dgrfdx = zeros(length(x),12);
                dgrfdx(:,[1,2,6,7,8,12]) = dgrfdx_tmp;
            else
                grf_tmp = obj.hdlMEX('GRF', x);
            end
            grf = zeros(12,1);
            grf([1,2,6,7,8,12],1) = grf_tmp;

        end


        %======================================================================
        %> @brief Function returns position and orientation of body segments, for the system in position q
        %>
        %> @details
        %> Needed for tracking of marker trajectories, and for 3D visualization of the model.
        %> This function calls the mex file of gait2d_osim.c
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
        %> @param   obj      Gait2d_osim class object
        %> @param   q        Double array: Generalized coordinates (elements of states with type 'q') (Gait2d_osim.nDofs x 1)
        %> @param   qdot     (optional) Double array: Generalized velocities (elements of states with type 'qdot') (Gait2d_osim.nDofs x 1)
        %>
        %> @retval  FK		 Double array: position (2) and orientation (2x2, stored row-wise) of nSegments-1.
        %>                   Segments are the bodies in the .osim model, not including ground. ((Gait2d_osim.nSegments-1) * 6)
        %> @retval	dFKdq	 (optional) Double matrix: Jacobian of FK with respect to q ((Gait2d_osim.nSegments-1) * 6 x Gait2d_osim.nDofs)
        %> @retval  dFKdotdq (optional) Double matrix: Jacobian of dFK/dt with respect to q ((Gait2d_osim.nSegments-1) * 6 x Gait2d_osim.nDofs)
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
        %> This function calls the mex file of gait2d_osim.c:
        %> [forces, lengths, momentarms] = obj.hdlMEX('Muscleforces', x);
        %>
        %> @param   obj         Gait2d_osim class object
        %> @param   x           Double array: State of the model (Gait2d_osim.nStates x 1)
        %>
        %> @retval  forces      Double array: Muscle forces (N) (Gait2d_osim.nMus x 1)
        %> @retval	lengths     (optional) Double array: Muscle-tendon lengths (m) (Gait2d_osim.nMus x 1)
        %> @retval	momentarms  (optional) Double matrix: Muscle moment arms (m) (Gait2d_osim.nMus x Gait2d_osim.nDofs)
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
        %> This function calls the mex file of gait2d_osim.c:
        %> [forces] = obj.hdlMEX('MuscleCEforces', x);
        %>
        %> @param   obj             Gait2d_osim class object
        %> @param   x               Double array: State of the model (Gait2d_osim.nStates x 1)
        %> @param   xdot            Double array: State derivatives (Gait2d_osim.nStates x 1)
        %>
        %> @retval  forcesCE        Double array: Muscle forces of CE (in N) (Gait2d_osim.nMus x 1)
        %> @retval  dforcesCEdx     Double array: Transpose of Jacobian matrix d force / d x (Gait2d_osim.nStates x Gait2d_osim.nMus)
        %> @retval  dforcesCEdxdot  Double array: Transpose of Jacobian matrix d force / d xdot (Gait2d_osim.nStates x Gait2d_osim.nMus)
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
        %> @param   obj       Gait2d_osim class object
        %> @param   x         Double array: State of the model (Gait2d_osim.nStates x 1)
        %> @param   xdot      Double array: State derivatives (Gait2d_osim.nStates x 1)
        %>
        %> @retval  powers	  Double array: CE power output (W) of each muscle (Gait2d_osim.nMus x 1)
        %> @retval  conDynRes Double array: Contraction dynamics residuals (Gait2d_osim.nMus x 1)
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
        %> @param   obj     Gait2d_osim class object
        %> @param   x       Double array: State of the model (Gait2d_osim.nStates x 1)
        %> @param   u       (optional) Double array: Controls of the model which are only needed if you want
        %>                  to apply extra moments (i.e. arm torques) (Gait2d_osim.nControls x 1)
        %>
        %> @retval  M       Double array: Moment/force for each DOF (Gait2d_osim.nDofs x 1)
        %> @retval  dMdx	(optional) Double matrix: Transpose of Jacobian matrix dM/dx (Gait2d_osim.nStates x Gait2d_osim.nDofs)
        %> @retval  dMdu    (optional) Double matrix: Transpose of Jacobian matrix dM/du (Gait2d_osim.nControls x Gait2d_osim.nDofs)
        %======================================================================
        function  [M, dMdx, dMdu] = getJointmoments(obj, x, u)

            % Get extra moments
            Mextra = zeros(obj.nDofs, 1);
            if nargin > 2 % Extra moments should be applied
                % Get extra moments (arm torques)
                idxTorque = obj.extractControl('torque');
                Mextra(obj.hidxTorqueDof) = obj.mExtraScaleFactor * u(idxTorque); % Scale them by obj.mExtraScaleFactor and assume that order is consistent.
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
                    dMdu(idxTorque, obj.hidxTorqueDof) = dMdMextra(obj.hidxTorqueDof, obj.hidxTorqueDof) * obj.mExtraScaleFactor; % scaling has to be considered
                end
            end

        end

        %======================================================================
        %> @brief Function to show model as stick figure
        %>
        %> @param   obj              Gait2d_osim class object
        %> @param   x                Double matrice: State vector of model for n time points (Gait2d_osim.nStates x n)
        %> @param   range            (optional) Double matrice: Defining the range of the figure
        %>                           with [xmin, xmax; ymin, ymax, zmin, zmax]. (3 x 2)
        %>                           (default: [pelvisX-1, pelvisX+1; -0.2, 2, pelvisZ-1, pelvisZ+1])
        %> @param   plotGRF          (optional) Bool: If true, it plots arrows for the GRFs (default: 0)
        %> @param   plotCPs          (optional) Bool: If true, it plots spheres for the CPs (default: 0)
        %> @param   plotJointCOSYs   (optional) Bool: If true, it plots the coordinate systems of the joints (default: 0)
        %======================================================================
        function showStick(obj, x, range, plotGRF, plotCPs, plotJointCOSYs)
            % make coordinates for a 1 cm radius sphere
            [xs,ys] = sphere(10);
            xs = xs*0.01;
            ys = ys*0.01;

            % Define positions for left and right leg
            l_pos = obj.joints.location(5,3);
            r_pos = obj.joints.location(2,3);

            % get range of figures
            if nargin > 2 && size(range, 1) == 3 && size(range, 2) == 2
                xrange = range(1, :);
                yrange = range(2, :);
            else
                idxPelvisX = obj.extractState('q', 'pelvis_tx');
                xrange = [min(x(idxPelvisX,:))  , max(x(idxPelvisX,:))];
                yrange = [-0.2, 2];
            end

            % plot the ground
            hold on;
            fill([xrange(1), xrange(1), xrange(2), xrange(2)], ... % forwards
                [0, 0, 0, 0], ...                                 % upwards
                [0.5, 0.5, 0.5], 'EdgeColor', [0.5, 0.5, 0.5]);

            % plot the stick figure
            for iTime = 1 : size(x,2)
                % run the forward kinematics
                fk = obj.getFkin(x(1:obj.nDofs, iTime));
                fk = reshape(fk, [6, obj.nSegments-1]);			% forward kinematics output from MEX function does not have ground

                for i=1:obj.nSegments-1

                    R = zeros(3,3);
                    O = zeros(3,1);

                    O_in = fk(1:2,i);                           % position of origin
                    R_in = reshape(fk(3:6,i)', 2,2)';			% rotation matrix

                    R(1:2,1:2) = R_in;
                    O(1:2) = O_in;
                    R(3,3) = 1;
                    % Draw Pelvis
                    if i == 1
                        fktemp = fk;
                        fktemp(3,:) = 0;
                        fktemp(3,[2 5]) = [r_pos l_pos];
                        idx = [1 2 5;...
                            2 5 8;...
                            5 8 1];
                        for j = 1:3
                            pg = fktemp(1:2,idx(j,:));
                            p = fill(pg(1,:),pg(2,:),[128 128 128]/255, 'edgecolor', 'k');
                            %p.FaceAlpha = 0.5;
                        end
                    elseif i == 3 || i == 4
                        line = [O(1:2)' r_pos; O_old(1:2)' r_pos]';
                        plot(line(1,:),line(2,:),'r')
                    elseif i == 6 || i == 7
                        line = [O(1:2)' l_pos; O_old(1:2)' l_pos]';
                        plot(line(1,:),line(2,:),'b')
                    elseif i == 8
                        foot_poly = [0 0 0; 2*obj.segments.mass_center(i+1,1:2) obj.joints.location(2,3);2*obj.segments.mass_center(i+1,1:2) obj.joints.location(5,3)];

                        pg = repmat(O,1,size(foot_poly,2)) + R*foot_poly';
                        fill(pg(1,:),pg(2,:),[128 128 128]/255, 'edgecolor', 'k');
                    end
                    % Plot a polygon for the foot, using the default cp locations
                    if i == 4 || i == 7
                        j = mod(i,2)+1; %zip([4 7],[1 2])
                        colors = ['r';'b'];
                        foot_poly = [0 0 obj.joints.location(i-2,3); obj.CPs.position(j,1:2) obj.joints.location(i-2,3);obj.CPs.position(j+2,1:2) obj.joints.location(i-2,3)];
                        pg = repmat(O,1,size(foot_poly,2)) + R*foot_poly';
                        fill(pg(1,:),pg(2,:),colors(j), 'edgecolor', colors(j));
                    end

                    O_old = O;

                    % draw the XY axes of the 18 segments in red, green, blue
                    if nargin > 6 && plotJointCOSYs
                        axislength = 0.1;
                        X = O + axislength*R(:,1);				% end of X axis
                        Y = O + axislength*R(:,2);				% end of Y axis
                        % plot ZXY as XYZ so Y will be up in the Matlab Window
                        plot([O(1) X(1)],[O(2) X(2)],'r','LineWidth',2)	% show X axis
                        plot([O(1) Y(1)],[O(2) Y(2)],'g','LineWidth',2)	% show Y axis
                    end



                end

                % get the ground reaction forces and plot them
                if nargin > 4 && plotGRF
                    grf = obj.getGRF(x(:, iTime));
                    [CoP_r, CoP_l] = obj.getCoP(grf);
                    Gait2d_osim.plotgrf(grf(1:3), CoP_r, 'r');  % right side GRF vector
                    Gait2d_osim.plotgrf(grf(7:9), CoP_l, 'b');  % left side GRF vector
                end

                % plot the contact points as spheres
                if nargin > 5 && plotCPs
                    xc_idx = obj.extractState('xc');
                    yc_idx = obj.extractState('yc');
                    for i = 1 : obj.nCPs
                        xc = xs + x(xc_idx(i), iTime);
                        yc = ys + x(yc_idx(i), iTime);
                        surf(xc,yc,'FaceColor',cpColor, 'EdgeColor', 'none');
                    end
                end

            end

            % finish the plot
            grid on; box on;
            xlabel('X');
            ylabel('Y');
            if abs(xrange(1) - xrange(2)) < 0.1
                xrange(2) = xrange(2)+2.2;
                xrange = xrange - 1.1;
            end
            if yrange(1)==yrange(2)
                yrange(2) = yrange(2)+2.2;
                yrange = yrange - 1.1;
            end
            axis equal
            ylim(yrange);
            xlim(xrange);
            hold off;

        end

        %======================================================================
        %> @brief Function to visualize treadmill as moving scatter plot
        %>
        %> @details
        %>
        %> @param   obj            Gait2d_osim class object
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
                assert(length(x_points) == nDots, 'Gait2d_osim:showTreadmill', 'x_points does not have the correct length of %i.', nDots);
                assert(length(y_points) == nDots, 'Gait2d_osim:showTreadmill', 'y_points does not have the correct length of %i.', nDots);
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
        %> @param   obj            Gait2d_osim class object
        %> @param   x              Double matrix: State vector of model for n time points (Gait2d_osim.nStates x n)
        %> @param   markerTable    Table: Variables table specifying markers with at least the columns:
        %>                         type, name, segment, position, direction to call Gait2d_osim.simuMarker.
        %> @param   measuredMean   (optional) Double matrix: Measured marker data in meter(!) which will be plotted as reference.
        %>                         The matrix must have the same number of time points as the state vector x.
        %>                         Further, the matrix columns must correspond to the rows of the markerTable in the same order
        %>                         to be able to match coordinates of the markers.  (n x height(markerTable)) (default: empty)
        %======================================================================
        function showMarker(obj,x, markerTable, measuredMean)

            if nargin > 3 && ~isempty(measuredMean)
                assert(size(measuredMean, 1) == size(x, 2), ...
                    'Gait2d_osim:showMarker(): The matrix must have the same number of time points as the state vector x.');
                assert(size(measuredMean, 2) == height(markerTable), ...
                    'Gait2d_osim:showMarker(): The matrix columns must correspond to the rows of the markerTable.');
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
                scatter3(0, markerSim(idxMarker(1, :)), markerSim(idxMarker(2, :)), 'ok'); % order as in showStick()

                % plot measured markers
                if plotMean
                    scatter3(0, measuredMean(iNode, idxMarker(1, :)), measuredMean(iNode, idxMarker(2, :)), 'xk');
                end
            end
            hold off;

        end


        % Declared method to track acc and gyro data
        [s, ds_dq, ds_dqd, ds_dqdd] = simuAccGyro(obj, data, q, qd, qdd);

        % Declared function to create graphical report fom a simulation result
        reportResult(obj, problem, result, resultfilename);

        %> @cond DO_NOT_DOCUMENT
        % Getter and setter function:
        %======================================================================
        %> @brief Function returning Gait2d_osim.nMarkers
        %>
        %> @param   obj     Gait2d_osim class object
        %> @retval  n       Double: Number of markers
        %======================================================================
        function n = get.nMarkers(obj)
            n = size(obj.markers,1);
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


    methods (Access = protected)
        %======================================================================
        %> @brief Function to initialize model with default parameters
        %>
        %> @details
        %> Initialialize contact points, alignment and stiffness of dofs
        %>
        %> @param   obj     Gait2d_osim class object
        %======================================================================
        function initModel(obj)
            % create contact points directly from markers
            % here a more complex shoe geometry could be added
            % set default alignment

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
            FT(ismember(muscleNames, {'hamstrings_r'  , 'hamstrings_l'  })) = 0.35;
            FT(ismember(muscleNames, {'bifemsh_r'  , 'bifemsh_l'  })) = 0.1;
            FT(ismember(muscleNames, {'glut_max_r', 'glut_max_l'})) = 0.45;
            FT(ismember(muscleNames, {'iliopsoas_r'    , 'iliopsoas_l'    })) = 0.5;
            FT(ismember(muscleNames, {'gastroc_r'      , 'gastroc_l'      })) = 0.5;
            FT(ismember(muscleNames, {'rect_fem_r' , 'rect_fem_l' })) = 0.65;
            FT(ismember(muscleNames, {'vasti_r'  , 'vasti_l'  })) = 0.5;
            FT(ismember(muscleNames, {'soleus_r'   , 'soleus_l'   })) = 0.2;
            FT(ismember(muscleNames, {'tib_ant_r'  , 'tib_ant_l'  })) = 0.25;
            obj.muscles.FT = FT;

            if ~obj.nMus==18
                warning('Gait2d_osim:initModel', 'Percentage of fast twitch fibers is not defined for all muscles of this model. This will be needed to compute metabolic cost.');
            end

            % Simplify the OpenSim foot model for the simulation
            obj.simplifyFoot();

            % Constrain the knee angle to 0 degrees
            obj.dofs{5,3}(2)=0;
            obj.dofs{8,3}(2)=0;

            obj.muscles.flen = obj.muscles.flen.*0+1.4; % Counter the millard muscles
            obj.muscles.kactive = obj.muscles.kactive.*0+0.5; % Millard kactive
        end


        %======================================================================
        %> @brief Initialize model mex file
        %>
        %> @details
        %> This function calls the mex file of gait2d_osim.c
        %> [info] = obj.hdlMEX('Initialize', model)
        %>
        %> This initializes the model.  This is required before anything is done with the model.
        %>
        %> Input:
        %>	model		Struct containing model parameters.
        %>
        %> Output:
        %>	info		Struct with information about the size of the model, fields:
        %>   	Nx				Number of state variables: 2* Gait2d_osim.nDofs + 2* Gait2d_osim.nMuscles + 9* Gait2d_osim.nCPs
        %>		Nf				Number of functions returned by the Dynamics function
        %>		fmin,fmax		Lower and upper bounds for dynamics residuals f
        %>      Bodyweight      Sum of weights of all segments. This is used to set bodymass
        %>
        %>
        %>
        %> @param   obj     Gait2d_osim class object
        %======================================================================
        function initMex(obj)
            obj.init = 0;

            % make struct for mex input
            % (Extracting parameters explicitly makes more clear, what is
            % required by gait2d_osim.c)
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
        %> @param   obj         Gait2d_osim class object
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
                Gait2d_osim.readOsim(osimfile);
            end
            load(matName,'model','osim_sha256')
            if (osim_sha256 ~= osim_sha256_new) & ~strcmp(computer,'MACA64') % OpenSim is not supported on Apple Silicon Mac
                model = Gait2d_osim.readOsim(osimfile);
            end
            obj.dofs = model.dofs;
            obj.joints = model.joints;
            obj.segments = model.segments;
            obj.markers = model.markers;
            obj.muscles = model.muscles;
            obj.gravity = model.gravity;
            obj.CPs = model.CPs;
            obj.osim = model.osim;
            obj.osim.sha256 = osim_sha256_new;
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
        %> @param   obj         Gait2d_osim class object
        %======================================================================
        function loadMomentArms(obj)

            % define ranges which are used to compute muscle arms
            range_muscleMoment = table(NaN(height(obj.dofs), 2), 'VariableNames', {'range_muscleMoment'}, 'RowNames', obj.dofs.Properties.RowNames);
            if strcmp(obj.osim.name, 'gait10dof18musc_wrap.osim') || strcmp(obj.osim.name, 'gait2d.osim') || strcmp(obj.osim.name, 'gait10dof18musc_wrap_ll.osim')
                % all our usual 2D_osim model
                range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'} = [- 60, 110] / 180 * pi;
                range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'} = [-160, 0] / 180 * pi;
                range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'} = [- 60, 60] / 180 * pi;
                range_muscleMoment{'hip_flexion_l'   , 'range_muscleMoment'} = range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'};
                range_muscleMoment{'knee_angle_l'    , 'range_muscleMoment'} = range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'};
                range_muscleMoment{'ankle_angle_l'   , 'range_muscleMoment'} = range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'};
                % No lumbar extending/flexing muscles
                %range_muscleMoment{'lumbar_extension', 'range_muscleMoment'} = [- 38,  3] / 180 * pi;
            else % default ranges
                %display('Standard gait2d_osim model')
                range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'} = [- 30, 160] / 180 * pi;
                range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'} = [-160, 0] / 180 * pi;
                range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'} = [- 60, 60] / 180 * pi;
                range_muscleMoment{'hip_flexion_l'   , 'range_muscleMoment'} = range_muscleMoment{'hip_flexion_r'   , 'range_muscleMoment'};
                range_muscleMoment{'knee_angle_l'    , 'range_muscleMoment'} = range_muscleMoment{'knee_angle_r'    , 'range_muscleMoment'};
                range_muscleMoment{'ankle_angle_l'   , 'range_muscleMoment'} = range_muscleMoment{'ankle_angle_r'   , 'range_muscleMoment'};
            end

            % load or compute moment arms
            filename = strrep(obj.osim.file,'.osim','_momentarms.mat');
            filename = strsplit(filename,{'/','\'});
            filename = which(filename{end});
            if ~exist(filename, 'file')
                disp('Momentarm file could not be found. Momentarms will be computed.')
                [momentarm_model] = obj.computeMomentArms(range_muscleMoment); % compute moment arms
            else
                load(which(filename),'momentarm_model','osim_sha256');
                if (osim_sha256 ~= obj.osim.sha256) & ~strcmp(computer,'MACA64') % OpenSim is not supported on Apple Silicon Mac
                    %b = input('The momentarm file is not up-to-date. Do you want new polynomials? (Y or N)','s');
                    display('Calculating new polynomials')
                    b = 'y';
                    if (b(1) == 'Y') || (b(1) == 'y')
                        [momentarm_model] = obj.computeMomentArms(range_muscleMoment);
                    end
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
            % moments (by default equal to range_muscleMoment; Are used in gait2d_osim.c)
            range_passiveMoment = range_muscleMoment; % default equal to range_muscleMoment; arm dofs are also defined in computeMomentArms()
            range_passiveMoment.Properties.VariableNames{strcmp(range_passiveMoment.Properties.VariableNames, 'range_muscleMoment')} = 'range_passiveMoment';

            dofs_tmp = [range_passiveMoment, dofs_tmp];

            if obj.nDofs == 10
                dofs_tmp.range_passiveMoment(10,:) = obj.dofs.range_osim(10,:);
            end
            % reasign it to the dofs table
            obj.dofs = dofs_tmp;
        end

        %======================================================================
        %> @brief Function to add Listeners to the obj
        %>
        %> @details
        %> - PostSet event listeners enable reinitialition for mex with new parameters.
        %>
        %> @param   obj     Gait2d_osim class object
        %======================================================================
        function add_listeners(obj)
            addlistener(obj,'gravity'           ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'drag_coefficient'  ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'wind_speed'        ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'joints'            ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'muscles'           ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'CPs'               ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'dofs'              ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'segments'          ,'PostSet',@obj.update_mexParameter);
            addlistener(obj,'segmentsTable'     ,'PostSet',@obj.update_mexParameter);
        end

        %======================================================================
        %> @brief Function performed to update the parameter of the mex and the tables
        %>
        %> @details
        %> - Called with a PostSet event after set Gait2d_osim.gravity, Gait2d_osim.drag_coefficient,
        %>   Gait2d_osim.wind_speed, Gait2d_osim.joints,
        %>   Gait2d_osim.muscles, Gait2d_osim.CPs, Gait2d_osim.dofs, Gait2d_osim.segments.
        %> - Reinitializes mex with new parameters.
        %> - This function calls all the other update functions. If we
        %>   would also add a listener for them, we could not ensure that they
        %>   are called in the correct order!
        %>
        %> @param   obj     Gait2d_osim class object
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
                case 'segmentsTable'
                    obj.segments = obj.segmentsTable;
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
                    obj.update_idxcxCP;
                    obj.update_constraints;
                case 'torques'
                    obj.hnTor = height(obj.torques);
                    obj.update_idxTorqueDof;     % Has to be updated first
                    obj.update_idxSymmetry;
                    obj.update_controls;
                    obj.update_idxcxCP;
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
                    obj.update_idxcxCP;
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
                    obj.update_controls;         % hidxArmdof have to be updated before
                    obj.update_idxControls;
                    obj.update_constraints;
                    obj.update_idxcxCP;
                case 'slope'
                    assert(~(obj.speed_right ~= 0  && obj.speed_left ~= 0 && obj.slope ~= 0), ...
                        'Treadmill speed with slope is not yet implemented');
            end

            % Reinitialize mex with new parameter matrix
            if ~isempty(obj.init) && obj.init
                obj.initMex;
            end
        end

        %======================================================================
        %> @brief Function defining the table Gait2d_osim.states
        %>
        %> @details
        %> Table lists type, name, xmin, xmax and xneutral
        %>
        %> @todo
        %> Variables of contact points of neutral position shouldn't be zero.
        %>
        %> @param   obj     Gait2d_osim class object
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
                xneutral(1:obj.nDofs) = obj.dofs.neutral_position; % It is defined in gait2d_osim.osim as default_value: Changing the default_value will change the passive joint moment!
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
                types = {'Fx', 'Fy', 'xc', 'yc'};%, 'xf', 'yf'}; %Fx,Fy,Fz (in BW), xc,yc,zc,xf,yf,zf (in m)
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
        %> @brief Function defining Gait2dc.idxcxCPleft and Gait2dc.idxcxCPright
        %>
        %>
        %> @param   obj     Gait2dc class object
        %======================================================================
        function update_idxcxCP(obj)
            
            if ~isempty(obj.CPs)
                obj.hidxcxCPleft = obj.extractState('xc',obj.CPs.Properties.RowNames(obj.CPs.segmentindex==8));
                obj.hidxcxCPright = obj.extractState('xc',obj.CPs.Properties.RowNames(obj.CPs.segmentindex==5));
            end

        end

        %======================================================================
        %> @brief Function defining the table Gait2d_osim.controls
        %>
        %> @details
        %> Table lists type, name, xmin, xmax and xneutral
        %>
        %> @param   obj     Gait2d_osim class object
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

            iarms = obj.hidxTorqueDof;
            if ~isempty(iarms)
                armNames = obj.dofs.Properties.RowNames(iarms);
                type(obj.nMus+1:obj.nMus+length(iarms)) = {'torque'};
                name(obj.nMus+1:obj.nMus+length(iarms)) = armNames;
                xmin(obj.nMus+1:obj.nMus+length(iarms)) = -5;
                xmax(obj.nMus+1:obj.nMus+length(iarms)) = 5;
                xneutral(obj.nMus+1:obj.nMus+length(iarms))= 0;
            end
            obj.hcontrols = table(type',name',xmin',xmax',xneutral');
            obj.hcontrols.Properties.VariableNames = {'type','name','xmin','xmax','xneutral'};
            obj.hnControls = height(obj.controls);

        end

        %======================================================================
        %> @brief Function defining the table Gait2d_osim.constraints
        %>
        %> @details
        %> Table lists type, name, equation (string), fmin and fmax
        %>
        %> The constraints are composed of Gait2d_osim.nConstraints = 2* Gait2d_osim.nDofs + 2* Gait2d_osim.nMus + 13* Gait2d_osim.nCPs)
        %>   - implicite differential equations: qdot-dq/dt = 0 (Gait2d_osim.nDofs x 1)
        %>   - equations of motion from Autolev (Gait2d_osim.nDofs x 1)
        %>   - muscle contraction dynamics (Gait2d_osim.nMus x 1)
        %>   - muscle activation dynamics: da/dt - rate * (u-a) = 0 (Gait2d_osim.nMus x 1)
        %>   - nine equality and four inequality constraints for contact points (Gait2d_osim.nCPs*13 x 1)
        %>
        %> @param   obj          Gait2d_osim class object
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
                nEqu = 4; % each contact point has this many equality constraints
                nInEqu = 0; % each contact point has this many inequality constraints
                nTotal = nEqu + nInEqu;
                type    (obj.nDofs*2+obj.nMus*2 + (1:obj.nCPs*nTotal)) = {'CP'};
                for iEqu = 1 : nEqu
                    name    (obj.nDofs*2+obj.nMus*2 + (iEqu:nTotal:obj.nCPs*nTotal)) = contactNames;
                    equation(obj.nDofs*2+obj.nMus*2 + (iEqu:nTotal:obj.nCPs*nTotal)) = {['equality ' num2str(iEqu)]};
                end
                for iInEqu = 1 : nInEqu
                    idx = iInEqu + nEqu;
                    name    (obj.nDofs*2+obj.nMus*2 + (idx:nTotal:obj.nCPs*nTotal)) = contactNames;
                    equation(obj.nDofs*2+obj.nMus*2 + (idx:nTotal:obj.nCPs*nTotal)) = {['inequality ' num2str(iInEqu)]};
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
        %> @brief Function defining Gait2d_osim.idxSymmetry
        %>
        %> @param   obj     Gait2d_osim class object
        %======================================================================
        function update_idxSymmetry(obj)
            % create the symmetry operator
            % x = xsign * x(xindex) will switch left and right in the state vector x
            % (assuming that the model walks in the +X direction)
            % u = u(uindex) switches left and right in the controls u
            % M = Msign * M(Mindex) switches left and right in the applied forces and moments M
            % generalized coordinates
            i_dofs = 1:obj.nDofs;
            s_dof = ones(1,obj.nDofs);
            for iDof_r = 1:obj.nDofs
                %find right dofs
                if ismember(obj.dofs.Properties.RowNames{iDof_r}(end-1:end),'_r')
                    %find corresponding left dof
                    %[~,iDof_l] = Gait2d_osim.getElementbyName(obj.dofs,[obj.dofs{iDof_r}.Properties.RowNames(1:end-2),'_l']);
                    iDof_l = find(strcmp(obj.dofs.Properties.RowNames,[obj.dofs.Properties.RowNames{iDof_r}(1:end-2),'_l']));
                    %switch right and left
                    i_dofs(iDof_r) = iDof_l;
                    i_dofs(iDof_l) = iDof_r;
                end
                % sign changes are needed in pelvis list, pelvis rotation, pelvis z,
                % lumbar bending and lumbar rotation
                if nnz(strcmp(obj.dofs.Properties.RowNames{iDof_r},{'pelvis_list','pelvis_rotation','pelvis_obliquity','pelvis_tz','lumbar_bending','lumbar_rotation'}))
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
                    %[~,imuscles_l] = Gait2d_osim.getElementbyName(obj.muscles,[obj.muscles{imuscles_r}.Properties.RowNames(1:end-2),'_l']);
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
            s_contact(strcmp(types, 'zf')) = -1;

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

        %======================================================================
        %> @brief Function to convert OpenSim's complicated foot into a singular bone
        %>
        %> @details
        %> Adds the mass, inertial properties of the talus, calcn, toes
        %>
        %> @param   F      Double array: Force (3 x 1)
        %> @param   CoP    Double array: Center of pressure in x, y, z coordinates (3 x 1)
        %> @param   color  String: Color (e.g. color = 'r')
        %======================================================================
        function simplifyFoot(obj)
            %
            seg = obj.segments; % Short name to refer to (lazy typer)
            %% Some naughty coded adjustments where the com is with respect to the ankle joint
            % ToDo: ask someone who knows if that is correct
            seg{'calcn_r','mass_center'} = seg{'calcn_r','mass_center'} + obj.joints{'subtalar_r','location'};
            seg{'toes_r','mass_center'} = seg{'toes_r','mass_center'} + obj.joints{'subtalar_r','location'} + obj.joints{'mtp_r','location'};
            seg{'calcn_l','mass_center'} = seg{'calcn_l','mass_center'} + obj.joints{'subtalar_l','location'};
            seg{'toes_l','mass_center'} = seg{'toes_l','mass_center'} + obj.joints{'subtalar_l','location'} + obj.joints{'mtp_l','location'};

            %% ToDo: The center of mass is not in a global coordinate system, fix that
            %% Right Foot
            mass_r = seg{'talus_r','mass'}+seg{'calcn_r','mass'}+seg{'toes_r','mass'};
            % CoM calculations:
            CoM_r = [];
            for i = 1:3
                com_i = (seg{'talus_r','mass'}*seg{'talus_r','mass_center'}(1,i)+...
                    seg{'calcn_r','mass'}*seg{'calcn_r','mass_center'}(1,i)+...
                    seg{'toes_r','mass'}*seg{'toes_r','mass_center'}(1,i))/mass_r;
                CoM_r(i) = com_i;
            end
            % Inertia calculations:
            inertia_r = zeros(3,3);
            for i = 1:3
                j = mod(i+1,3);
                k = mod(i+2,3);
                j(j==0)=3;
                k(k==0)=3;
                inertia_r(i,i) = seg{'talus_r','inertia'}{1}(i,i) + ...
                    seg{'calcn_r','inertia'}{1}(i,i)+seg{'toes_r','inertia'}{1}(i,i) + ... %Sum of inertial masses
                    +seg{'talus_r','mass'}*((seg{'talus_r','mass_center'}(1,k)-CoM_r(k))^2) + ...
                    +seg{'calcn_r','mass'}*((seg{'calcn_r','mass_center'}(1,k)-CoM_r(k))^2) + ...
                    +seg{'toes_r','mass'}*((seg{'toes_r','mass_center'}(1,k)-CoM_r(k))^2)+ ... % Steiner Anteile
                    +seg{'talus_r','mass'}*((seg{'talus_r','mass_center'}(1,j)-CoM_r(j))^2) + ...
                    +seg{'calcn_r','mass'}*((seg{'calcn_r','mass_center'}(1,j)-CoM_r(j))^2) + ...
                    +seg{'toes_r','mass'}*((seg{'toes_r','mass_center'}(1,j)-CoM_r(j))^2); % Steiner Anteile #2 axis
            end
            %% Left Foot
            mass_l = seg{'talus_l','mass'}+seg{'calcn_l','mass'}+seg{'toes_l','mass'};
            % CoM calculations:
            CoM_l = [];
            for i = 1:3
                com_i = (seg{'talus_l','mass'}*seg{'talus_l','mass_center'}(1,i)+...
                    seg{'calcn_l','mass'}*seg{'calcn_l','mass_center'}(1,i)+...
                    seg{'toes_l','mass'}*seg{'toes_l','mass_center'}(1,i))/mass_l;
                CoM_l(i) = com_i;
            end
            % Inertia calculations:
            inertia_l = zeros(3,3);
            for i = 1:3
                j = mod(i+1,3);
                k = mod(i+2,3);
                j(j==0)=3;
                k(k==0)=3;
                inertia_l(i,i) = seg{'talus_l','inertia'}{1}(i,i) + ...
                    seg{'calcn_l','inertia'}{1}(i,i)+seg{'toes_l','inertia'}{1}(i,i) + ... %Sum of inertial masses
                    +seg{'talus_l','mass'}*((seg{'talus_l','mass_center'}(1,k)-CoM_l(k))^2) + ...
                    +seg{'calcn_l','mass'}*((seg{'calcn_l','mass_center'}(1,k)-CoM_l(k))^2) + ...
                    +seg{'toes_l','mass'}*((seg{'toes_l','mass_center'}(1,k)-CoM_l(k))^2)+ ... % Steiner Anteile
                    +seg{'talus_l','mass'}*((seg{'talus_l','mass_center'}(1,j)-CoM_l(j))^2) + ...
                    +seg{'calcn_l','mass'}*((seg{'calcn_l','mass_center'}(1,j)-CoM_l(j))^2) + ...
                    +seg{'toes_l','mass'}*((seg{'toes_l','mass_center'}(1,j)-CoM_l(j))^2); % Steiner Anteile #2 axis
            end
            %% Assigning Segments

            feet = table([mass_r;mass_l],[CoM_r;CoM_l],{inertia_r;inertia_l},...
                {'ankle_r';'ankle_l'},'VariableNames',{'mass','mass_center','inertia','parent_joint'},...
                'RowNames',{'foot_r','foot_l'});

            % Insert foot into the segments table
            seg([5,6,7,10,11,12],:) = [];
            seg = [seg;feet];
            new_order = [1,2,3,4,8,5,6,9,7];
            obj.segments = seg(new_order,:);
            
            %% Contact points, we use the location of the mtp joint for scale factors
            % Issue in OpenSim: Contact points are not affected by scaling
            generic_mtp = [0.1300   -0.0440]; 
            sf_l = (obj.joints{'subtalar_l','location'}(1:2) + obj.joints{'mtp_l','location'}(1:2))/generic_mtp;
            sf_r = (obj.joints{'subtalar_r','location'}(1:2) + obj.joints{'mtp_r','location'}(1:2))/generic_mtp;
          
            obj.CPs{'heel_l','position'} = obj.CPs{'heel_l','position'}*sf_l;
            obj.CPs{'heel_r','position'} = obj.CPs{'heel_r','position'}*sf_r;
            obj.CPs{'front_l','position'} = obj.CPs{'front_l','position'}*sf_l;
            obj.CPs{'front_r','position'} = obj.CPs{'front_r','position'}*sf_r;

            %% Assign correct segment names & positions to the Markers
            for i = 1:height(obj.markers)
                m_name = obj.markers.segment{i};
                switch m_name
                    case 'calcn_r'
                        obj.markers{i,'position'} = obj.markers{i,'position'}+obj.joints{'subtalar_r','location'};
                        obj.markers.segment{i} = 'foot_r';
                    case 'toes_r'
                        obj.markers{i,'position'} = obj.markers{i,'position'}+obj.joints{'subtalar_r','location'}+obj.joints{'mtp_r','location'};
                        obj.markers.segment{i} = 'foot_r';
                    case 'talus_r'
                        obj.markers.segment{i} = 'foot_r';                        
                    case 'calcn_l'
                        obj.markers{i,'position'} = obj.markers{i,'position'}+obj.joints{'subtalar_l','location'};
                        obj.markers.segment{i} = 'foot_l';
                    case 'toes_l'
                        obj.markers{i,'position'} = obj.markers{i,'position'}+obj.joints{'subtalar_l','location'}+obj.joints{'mtp_l','location'};
                        obj.markers.segment{i} = 'foot_l';
                    case 'talus_l'
                        obj.markers.segment{i} = 'foot_l';
                end
            end

            %% Delete unneccessary joints
            % ToDo: (here and above) search for the correct joint, not the
            % hardcoded index using obj.joints.Properties.RowNames
            obj.joints([5,6,10,11],:) = [];

            % Assign correct segment names
            obj.joints(4,16) = {'foot_r'};
            obj.joints(7,16) = {'foot_l'};

            %% ToDo: Feet segment indices are hardcoded
            obj.CPs.segmentindex = [5 8 5 8]';
        end

        % Declared function to compute moment arms using Opensim model
        [momentarms] = computeMomentArms(obj,range_muscleMoment,examine,name);
    end

    methods(Static, Access = protected)
        % Declared function to read Opensim file
        model = readOsim(osimfile)
        [opensimfile_tmp, hash] = scaleOsim(opensimfile, scale_factors, bodyweight, varargin)
    end

    methods(Static)
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

        % Declared function tomake MEX functions
        name_MEX = getMexFiles(osimName,clean, rebuild_contact, optimize)

        % Declared function to randomize model
        % obj = createRandomSubject(obj, symmusrand, massratio,lengthratio)

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
                Gait2d_osim.getMexFiles(obj.osim.name);
                obj.initMex;
                obj.add_listeners(); % When loading an object, listeners might get lost
            end
        end
        %> @endcond

    end
end



