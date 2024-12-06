% ======================================================================
%> @file @Gait2dc_Exo/Gait2dc_Exo.m
%> @brief Gait2dc class describing the Gait2dc_Exo model
%>
%> @author Anne Koelewijn
%> @date 2019
% ======================================================================

% ======================================================================
%> @brief The class describes the Gait2dc_Exo model
%>
%> @details
%> - Inherits from Gait2dc
%> - Adds an exoskeleton to one of the joints (both sides), which adds
%> a penalizing torque dependent on step frequency or center of mass
%>
%> @todo Write a test class inheriting from Gait2dcTest
% ======================================================================
classdef Gait2dc_Exo < Gait2dc
    
    properties
        %> 
        exoMoment
        type
        norm_value
        increase
    end
    
    methods

        %======================================================================
        %> @brief Default constructor setting default Gait2dc object
        %>
        %> @details
        %> Initializes the model and the mex function.
        %>
        %> The standard model can be called using:
        %> @code
        %> Gait2dc('gait2dc_par.xls')
        %> @endcode
        %> The excelfile must be in the matlab path.
        %>
        %> @param   varargin      Variable length input with:
        %>                          - excelfile     String: Excel file with path, name and extension
        %>                          - bodyheight    (optional) Double: Bodyheight in m
        %>                          - bodymass      (optional) Double: Bodymass in kg
        %> @retval  obj           Gait2dc class object
        %======================================================================
        function [obj] = Gait2dc_Exo(norm_value, type, increase, varargin)
            
            % call Gait2dc constructor
            obj = obj@Gait2dc(varargin{:});
            obj.norm_value = norm_value;
            obj.increase = increase;
            obj.type = type;
        end
        
        %======================================================================
        %> @brief Function to compute implicit differential equation for 2D musculoskeletal model
        %>
        %> @details
        %> This function calls the mex file of gait2dc.c:
        %> [f, dfdx, dfdxdot, dfdumus, dfdMextra]  = gait2dc('Dynamics',x,xdot,umus,Mextra);
        %>
        %> with the neural excitations umus and the extra (exoskeleton) torques Mextra.
        %>
        %> The dynamic residuals will be between fmin and fmax when inputs
        %> satisfy system dynamics: fmin <= f(x,dx/dt,umus,Mextra) <= fmax
        %>
        %> The last four outputs are optional and some computation time is saved if you do
        %> not request all of them.
        %>
        %> @todo Don't use fixed indices. Get them from the state vector table!
        %>
        %> @param   obj     Gait2dc class object
        %> @param   x       Double array: State of the model (Gait2dc.nStates x 1)
        %> @param   xdot    Double array: State derivatives (Gait2dc.nStates x 1)
        %> @param   u       Double array: Controls of the model (Gait2dc.nControls x 1)
		%> @param   input   Double: Value of the input to the exoskeleton controller (i.e., current duration)
        %>
        %> @retval  f       Double array: Dynamic residuals (Gait2dc.nConstraints x 1)
        %> @retval	dfdx	(optional) Double matrix: Transpose of Jacobian matrix df/dx 		(Gait2dc.nStates x Gait2dc.nConstraints)
        %> @retval	dfdxdot	(optional) Double matrix: Transpose of Jacobian matrix df/dxdot 	(Gait2dc.nStates x Gait2dc.nConstraints)
        %> @retval	dfdu	(optional) Double matrix: Transpose of Jacobian matrix df/du 		(Gait2dc.nControls x Gait3d.nConstraints)
        %======================================================================
        function [f, dfdx, dfdxdot, dfdu, dfdinput] = getDynamics(obj,x,xdot,u, input)
            
            noGlobDof = obj.nDofs-obj.nJoints;
            if strcmp(obj.type, 'dur')
                idx1 = find(strcmp(obj.states.name, 'knee_angle_r') ==1,2);
                idx1m = idx1(1);
                idx1 = idx1(2);
                idx2 = find(strcmp(obj.states.name, 'knee_angle_l') ==1,2);
                idx2m = idx2(1);
                idx2 = idx2(2);
                
				vel1 = x(idx1);
				vel2 = x(idx2);

				% Find exoskeleton moments and add to external moment vector
				[Mknee1, dM1ddur, dM1dvel1,dM1dacc1] = obj.findExoMoment(obj.norm_value, input, vel1, xdot(idx1), obj.increase);
				[Mknee2, dM2ddur, dM2dvel2,dM2dacc2] = obj.findExoMoment(obj.norm_value, input, vel2, xdot(idx2), obj.increase); 
				
				dMddur = [0;dM1ddur;0;0;dM2ddur;0];
                Mextra = [0;Mknee1;0;0;Mknee2;0]; % for joints but not for orientation and position
            end                        
            
            % Get dynamics   
            if nargout > 3
                [f, dfdx, dfdxdot, dfdu, dfdMextra] = gait2dc('Dynamics',x,xdot,u,Mextra);
                
                [dmdxdot, dmdx2] = deal(zeros(6,size(x,1)));
                dmdx2(idx1m-noGlobDof,idx1) = dmdx2(idx1m-noGlobDof,idx1)+dM1dvel1; %First index should not be hardcoded
                dmdx2(idx2m-noGlobDof,idx2) = dmdx2(idx2m-noGlobDof,idx2)+dM2dvel2;
                dmdxdot(idx1m-noGlobDof,idx1) = dmdxdot(idx1m-noGlobDof,idx1) + dM1dacc1; %Second index should not be hardcoded
                dmdxdot(idx2m-noGlobDof,idx2) = dmdxdot(idx2m-noGlobDof,idx2) + dM2dacc2;
              
                dfdx = dfdx + transpose(dfdMextra'*dmdx2);
                dfdxdot = dfdxdot + transpose(dfdMextra'*dmdxdot);
                if contains(obj.type, 'dur')
                    dfdinput = dfdMextra'*dMddur;
                elseif strcmp(obj.type, 'com')
                    dfdinput = dfdMextra'*dMdcom;
                end
            elseif nargout > 2
                [f, dfdx, dfdxdot] = gait2dc('Dynamics',x,xdot,u,Mextra);
            elseif nargout > 1 
                [f, dfdx] = gait2dc('Dynamics',x,xdot,u,Mextra);  
            else
                f = gait2dc('Dynamics',x,xdot,u,Mextra);
            end
        end
       
    end
    
    methods(Static)
        
        %======================================================================
        %> @brief Function to get the external moment from the exoskeleton
        %> @static
        %> @public
        %>
        %> @details
        %> We can use this static function to save coefficients we used in
        %> the past.
        %> 
        %>
        %> @param   sr_norm 	 Double: Reference step rate, normally the step rate
		%> 						 that is optimal for normal walking
        %> @param   dur     	 Double: Duration of current gait cycle
        %> @param   angvel  	 Double: Angular velocity of the knee
		%> @param   ddq			 Double: Angular acceleration of the knee
		%> @param   increase	 Boolean: 1 for penalize-high, 0 for penalize low condition
        %>
        %> @retval  Mexo    	 Double: Exoskeleton moment
        %> @retval	dMddur		 Double: gradient of exoskeleton moment wrt duration
        %> @retval	dMdangvel	 Double: gradient of exoskeleton moment wrt angular velocity
        %> @retval	dMexodxdot	 Double: gradient of exoskeleton moment wrt angular acceleration
        %======================================================================
        function [Mexo, dMddur, dMdangvel, dMexodxdot] = findExoMoment(sr_norm, dur, angvel, ddq, increase)

            if nargin < 5
                increase = 1; %penalizing torque increases with step rate
            end

            steprate = 30/dur;
            dsteprateddur = -30/dur^2;

            torque_c = 3.36/100; %torque constant from exoskeleton properties
            meanvel_c = 60; %constant to get around 6 Nm at the mean angular velocity

            %relation wrt step rate; 
            if increase
                sr_base = 0.85*sr_norm;
                slope = 1/(1.1*sr_norm - 0.85*sr_norm);
            else
                sr_base = 1.15*sr_norm;
                slope = 1/(0.9*sr_norm - 1.15*sr_norm);
            end
            
            %base equation for torque
            Mexo1 = torque_c*meanvel_c*angvel*((steprate-sr_base)*slope); 
            dMexo1dangvel = torque_c*meanvel_c*((steprate-sr_base)*slope);
            dMexo1dsr = torque_c*meanvel_c*angvel*slope;
            
            %Use only absolute value
            Mexo1_abs = sqrt(Mexo1^2+0.01^2);
            dMexo1dabs = Mexo1/Mexo1_abs;
            
            dMexo1adangvel = dMexo1dangvel*dMexo1dabs;
            dMexo1adsr = dMexo1dsr*dMexo1dabs;
            
            %Soft max at 12 Nm (previously 20)
            Mexo2 = -log(exp(-Mexo1_abs) + exp(-12) );
            dMexo2dMexo1 = exp(-Mexo1_abs)/(exp(-Mexo1_abs) + exp(-12) );

            %Opposite direction of torque
            Mexo = -Mexo2*atan(10*ddq)/pi*2; 
            dMexodMexo2 = -atan(10*ddq)/pi*2;
            dMexodxdot = -10*Mexo2/(1+(10*ddq)^2)/pi*2;
            dMexodMexo1 = dMexodMexo2*dMexo2dMexo1;

            dMddur = dMexodMexo1*dMexo1adsr*dsteprateddur;
            dMdangvel = dMexodMexo1*dMexo1adangvel;
        end
    end
end


