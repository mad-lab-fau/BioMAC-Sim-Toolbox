% ======================================================================
%> @file @Model/Model.m
%> @brief Matlab class (abstract) describing a musculoskeletal model
%>
%> @author Eva Dorschky, Marlies
%> @date November, 2017
% ======================================================================

% ======================================================================
%> @brief The abstract class describes the basics for a musculoskeletal model
%>
%> @todo We should set init = 0 as soon as we are loading an other model.
%> This should be also done for all models which are saved within an other
%> object (e.g. problem.model).
%> Furthermore, we should check init == 1 each time before calling the mex.
% ======================================================================
classdef (Abstract) Model < handle
    
    properties(SetAccess = protected)
        %> Double: Showing if the mex model is initialized 
        init
    end
    
    properties (Constant)
        %> Cell array with string: Gives the names for the GRF vector
        %> returned by getGRF()
        GRFNAMES = {'rightFx', 'rightFy', 'rightFz', 'rightMx', 'rightMy', 'rightMz', 'leftFx', 'leftFy', 'leftFz', 'leftMx', 'leftMy', 'leftMz' };  
    end
    
    properties (SetObservable, AbortSet, SetAccess = protected)
        %> Table: Segments (Segment properties should not be changed)
        segments
    end
    
    properties (SetObservable, AbortSet)
        %> Double array: Gravity vector in m/(s^2)
        gravity
        %> Double: Air drag coefficient in N/(m/s)^2  (default: 0.2128)
        drag_coefficient = 0.2128
        %> Double: Wind speed in direction of +X in m/s (default: 0)
        wind_speed = 0
        %> Double: Scale factor for extra torques in Nm (default: 100 since 
        %> we used this in previous simulations. However, it's probably overwritten
        %> in the construction (see Gait3d.m)!)
        mExtraScaleFactor = 100
        %> Table: Degrees of freedom (Range of Dofs should not be
        %> changed here. Set the bounds of the Problem instead.)
        dofs = table();
        %> Table: Joints
        joints = table();
        %> Table: Muscles
        muscles = table();
        %> Table: Torque actuators
        torques = table();
        %> Table: Contact points 
        CPs = table();
    end
    
    properties (Dependent, SetAccess = protected)
        %> Double: Number of muscles (height of Model.muscles)
        nMus
        %> Double: Number of torque actuators (height of Model.torques)
        nTor
        %> Double: Number of contact points (height of Model.CPs)
        nCPs
        %> Double: Number of degree of freedoms (height of Model.dofs)
        nDofs
        %> Double: Number of joints (height of Model.joints)
        nJoints
        %> Double: Number of segments (height of Model.segments)
        nSegments
        %> Table: Information on states of the model
    	states    
        %> Table: Information on controls of the model
        controls
        %> Table: Information on constraints implemented in the mex and here
        constraints
        %> Double: Number of states (height of Model.states)
        nStates
        %> Double: Number of controls (height of Model.controls)
        nControls
        %> Double: Number of constraints (height of Model.constraints)
        nConstraints
        %> Struct: Double arrays with indices for symmetry to map right to left
        idxSymmetry
        %> Double array: Indices for forward translation
        idxForward
        %> Double array: Indices for upward translation
        idxUpward
        %> Double array: Indices for sideward translation
        idxSideward
        %> Double array: Indices for forward translation, speed, position of contact points and force at contact points
        idxForwardAll
        %> Double array: Indices for sideward translation, speed, position of contact points and force at contact points
        idxSidewardAll
        %> Double array: Indices of dofs of torque actuators
        idxTorqueDof
    end
    
    %> @cond DO_NOT_DOCUMENT
    properties (Hidden, Access = protected) % use hidden states to store dependent variables to save computing time
        %> Hidden double: Number of muscles (height of Model.muscles)
        %> It has to be zero at the beginning before anything is added.
        hnMus = 0
        %> Hidden double: Number of torque actuators (height of Model.torques)
        %> It has to be zero at the beginning before anything is added.
        hnTor = 0
        %> Double: Number of contact points (height of Model.CPs)
        %> It has to be zero at the beginning before anything is added.
        hnCPs = 0
        %> Hidden double: Number of joints (height of Model.joints)
        %> It has to be zero at the beginning before anything is added.
        hnDofs = 0
        %> Hidden double: Number of joints (height of Model.joints)
        hnJoints = 0
        %> Hidden double: Number of segments (height of Model.segments)
        hnSegments = 0
        %> Hidden table: Information on states of the model
        hstates
        %> Hidden struct: Row indices for all types in states
        hidxStates
        %> Hidden double: Number of states (height of Model.states)
        hnStates = 0
        %> Hidden table: Information on controls of the model
        hcontrols
        %> Hidden struct: Row indices for all types in controls
        hidxControls
        %> Hidden double: Number of controls (height of Model.controls)
        hnControls
        %> Hidden table: Information on constraints implemented in the mex and here
        hconstraints
        %> Hidden double: Number of constraints (height of Model.constraints)
        hnConstraints = 0
        %> Hidden struct: Double arrays with indices for symmetry
        hidxSymmetry
        %> Hidden double array: Indices for forward translation
        hidxForward
        %> Hidden double array: Indices for upward translation
        hidxUpward
        %> Hidden double array: Indices for sideward translation
        hidxSideward
        %> Hidden double array: Indices for forward translation, speed, position of contact points and force at contact points
        hidxForwardAll
        %> Hidden double array: Indices for sideward translation, speed, position of contact points and force at contact points
        hidxSidewardAll
        %> Obsolete! Not used and only here for backwards compatability! (Hidden double array: Indices of dofs of arms)
        hidxArmdof
        %> Hidden double array: Indices of dofs of torque actuators
        hidxTorqueDof
    end
    %> @endcond

    methods (Abstract)
        %======================================================================
        %> @brief Abstract function to compute implicit differential equation for the model
        %======================================================================
        [f, dfdx, dfdxdot, dfdu] = getDynamics(obj,x,xdot,u)
        
        %======================================================================
        %> @brief Abstract function returning the ground reaction forces for the system in state x
        %======================================================================
        [grf, dgrfdx] = getGRF(obj, x)
           
        %======================================================================
        %> @brief Abstract function returning muscle forces for the system in state x
        %======================================================================
        [muscleForces] = getMuscleforces(obj, x)
        
        %======================================================================
        %> @brief Abstract function returning muscle forces of CE for the system in state x
        %======================================================================
        [muscleCEForces] = getMuscleCEforces(obj, x, xdot)
       
        %======================================================================
        %> @brief Abstract function returns power generated by muscle contractile elements, for the system in state x
        %======================================================================
        [powers] = getMuscleCEpower(obj, x, xdot)
      
        %======================================================================
        %> @brief Abstract function returns joint moments, for the system in state x
        %======================================================================
        [M, dMdx, dMdu] = getJointmoments(obj, x, u)
        
        %======================================================================
        %> @brief Abstract function to simulate acceleration and gyroscope signals
        %======================================================================
        [s, ds_dq, ds_dqd, ds_dqdd] = simuAccGyro(obj, data, q, qd, qdd)
        
        %======================================================================
        %> @brief Abstract function to show model as stick figure
        %======================================================================
        showStick(obj,x)
        
        
    end
    
    methods(Abstract, Access = protected)
        %======================================================================
        %> @brief Abstract function to initialize with default parameters
        %======================================================================
        initModel(obj, vargin)
        
        %======================================================================
        %> @brief Abstract function to initialize model mex file
        %======================================================================
        initMex(obj)
        
        %======================================================================
        %> @brief Abstract function performed to update the parameter of the mex and the tables
        %======================================================================
        update_mexParameter(obj,src,evnt)
        
        %======================================================================
        %> @brief Abstract function defining the table Model.states
        %======================================================================
        update_states(obj)
        
        %======================================================================
        %> @brief Abstract function defining the table Model.controls
        %======================================================================
        update_controls(obj)
        
        %======================================================================
        %> @brief Abstract function defining the table Model.constraints
        %======================================================================
        update_constraints(obj)
        
        %======================================================================
        %> @brief Abstract function defining Model.idxSymmetry
        %======================================================================
        update_idxSymmetry(obj)
        
    end
    
    
    methods (Access = protected)
        
        %======================================================================
        %> @brief Function defining Model.idxTorqueDof
        %>
        %> @details
        %> Searches for indices of dof which have torque actuators in
        %> Model.torques.
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxTorqueDof(obj)
            
            % Find indices of dofs with torques
            if ~isempty(obj.dofs) && ~isempty(obj.torques)
                iarms = find(ismember(obj.dofs.Properties.RowNames, obj.torques.dof));
            else
                iarms = [];
            end
            
            % Set indices
            obj.hidxTorqueDof = iarms;
           
        end
        
        %======================================================================
        %> @brief Function defining Model.hidxStates
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxStates(obj)

            % Get all types
            types = unique(obj.states.type, 'stable');

            % Go over all types and get indices
            idxStates = struct();
            for iType = 1 : length(types)
                idxStates.(types{iType}) = find(strcmp(obj.states.type, types{iType}));
            end

            % Set indices
            obj.hidxStates = idxStates;

        end

        %======================================================================
        %> @brief Function defining Model.hidxControls
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxControls(obj)

            % Get all types
            types = unique(obj.controls.type, 'stable');

            % Go over all types and get indices
            idxControls = struct();
            for iType = 1 : length(types)
                idxControls.(types{iType}) = find(strcmp(obj.controls.type, types{iType}));
            end

            % Set indices
            obj.hidxControls = idxControls;

        end

        %======================================================================
        %> @brief Function defining Model.idxForward
        %>
        %> @details
        %> Provides the indices of state variables that have forward translation.
        %> States have to be updated before!
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxForward(obj)
            
            % forward translation in dofs
            idxDofs = obj.extractState('q', 'pelvis_tx');
            
            % forward translation in global x of all contact points
            idxCPs = obj.extractState('xc')';
            
            obj.hidxForward = [idxDofs, idxCPs];
        end 
        
        %======================================================================
        %> @brief Function defining Model.idxUpward
        %>
        %> @details
        %> Provides the indices of state variables that have upward translation.
        %> States have to be updated before!
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxUpward(obj)
            
            % upward translation in dofs
            idxDofs = obj.extractState('q', 'pelvis_ty');
            
            % upward translation in global y of all contact points
            idxCPs = obj.extractState('yc')';
            
            obj.hidxUpward = [idxDofs, idxCPs];
        end 
        
        %======================================================================
        %> @brief Function defining Model.idxSideward
        %>
        %> @details
        %> Provides the indices of state variables that have sideward translation.
        %> States have to be updated before!
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxSideward(obj)
            
            % sideward translation in dofs
            idxDofs = obj.extractState('q', 'pelvis_tz');
            
            % sideward translation in global z of all contact points
            idxCPs = obj.extractState('zc')';
            
            obj.hidxSideward= [idxDofs, idxCPs];
        end
        
        %======================================================================
        %> @brief Function defining Model.idxForwardAll
        %>
        %> @details
        %> Provides the indices of state variables that have forward translation, 
        %> speed, position of contact points and force at contact points.
        %> States and idxForward have to be updated before!
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxForwardAll(obj)
            
            % forward translation speed in dofs
            idxDofsDot = obj.extractState('qdot', 'pelvis_tx');
            
            % forward force in global x of all contact points
            idxCPFs = obj.extractState('Fx')';
            
            obj.hidxForwardAll = [obj.idxForward, idxDofsDot, idxCPFs]; % We do not care about the order
        end 
        
        %======================================================================
        %> @brief Function defining Model.idxSidewardAll
        %>
        %> @details
        %> Provides the indices of state variables that have sideward translation, 
        %> speed, position of contact points and force at contact points.
        %> States and idxSideward have to be updated before!
        %>
        %> @param   obj     Model class object
        %======================================================================
        function update_idxSidewardAll(obj)
            
            % sideward translation speed in dofs
            idxDofsDot = obj.extractState('qdot', 'pelvis_tz');
            
            % sideward force in global x of all contact points
            idxCPFs = obj.extractState('Fz')';
            
            obj.hidxSidewardAll = [obj.idxSideward, idxDofsDot, idxCPFs]; % We do not care about the order
        end
        
        
    end

    methods
        
        %======================================================================
        %> @brief Function to obtain index of state with a specific type (and name)
        %>
        %> @param   obj         Model class object 
        %> @param   type        String: Type of the state
        %> @param   statename   (optional) String or cell array of strings: Name of the state
        %> @retval  iState      Double array: Indices in Model.states matching type (and name)
        %====================================================================== 
        function iState = extractState(obj,type,statename)
            if isempty(obj.hidxStates) % Only needed when old results are loaded
                obj.update_idxStates;
            end
            if ~isfield(obj.hidxStates, type)
                iState = [];
            elseif nargin == 2
                iState = obj.hidxStates.(type);
            elseif nargin == 3
                idxType = obj.hidxStates.(type);
                namesForType = obj.states.name(idxType);
                if ischar(statename)
                    iState = idxType(strcmp(namesForType, statename));
                elseif iscell(statename)
                    iState = zeros(size(statename));
                    for iName1 = 1 : size(statename, 1)
                        for iName2 = 1 : size(statename, 2)
                            iState(iName1, iName2) = idxType(strcmp(namesForType,statename{iName1, iName2}));
                        end
                    end
                end
            end
        end
        
        %======================================================================
        %> @brief Function to obtain index of contol with a specific type (and name)
        %>
        %> @param   obj         Model class object 
        %> @param   type        String: Type of the control
        %> @param   controlname (optional) String or cell array of strings: Name of the control
        %> @retval  iControl    Double array: Indices in Model.controls matching type (and name)
        %====================================================================== 
        function iControl = extractControl(obj,type,controlname)
            if isempty(obj.hidxControls) % Only needed when old results are loaded
                obj.update_idxControls;
            end
            if ~isfield(obj.hidxControls, type)
                iControl = [];
            elseif nargin == 2
                iControl = obj.hidxControls.(type);
            elseif nargin == 3
                idxType = obj.hidxControls.(type);
                namesForType = obj.controls.name(idxType);
                if ischar(controlname)
                    iControl = idxType(strcmp(namesForType, controlname));
                elseif iscell(controlname)
                    iControl = zeros(size(controlname));
                    for iName1 = 1 : size(controlname, 1)
                        for iName2 = 1 : size(controlname, 2)
                            iControl(iName1, iName2) = idxType(strcmp(namesForType, controlname{iName1, iName2}));
                        end
                    end
                end
            end
        end

        
        %======================================================================
        %> @brief Function to set segment mass
        %>
        %> @param   obj             Model class object 
        %> @param   segmentName     String: Name of the segment
        %> @param   segmentMass     Double: Mass of the segment in kg
        %====================================================================== 
        function setSegmentMass(obj, segmentName, segmentMass)
            obj.segments.mass(segmentName) = segmentMass;
        end
        
        
        %======================================================================
        %> @brief Function to set segment Table
        %>
        %> @param   obj             Model class object 
        %> @param   segmentTable     Table: segmentTable
       
        %=====================================================================
        function setSegments(obj,segmentTable)
            obj.segments=segmentTable;
        end
        
        %=====================================================================
        %> @brief Function to set body mass, does not influence the
        %         simulation
        %>
        %> @param   obj             Model class object 
        %> @param   segmentTable     Double: bodymass
       
        %=====================================================================
        function setbodymass(obj,bodymass)
            obj.bodymass=bodymass;
        end

        %=====================================================================
        %> @brief Function to set the muscles table
        %>
        %> @param   obj             Model class object 
        %> @param   segmentTable     Table: muscles
       
        %=====================================================================
        function setMuscles(obj,muscles)
            obj.muscles=muscles;
        end


        %> @cond DO_NOT_DOCUMENT  
        %======================================================================
        %> @brief Function to save object of Model 
        %>
        %> @details 
        %> Sets the flag Model.init to false
        %>
        %> @param   obj    Model class object which sould be saved
        %> @retval  sobj   Model class object which is saved
        %====================================================================== 
        function sobj = saveobj(obj)
            sobj = obj;
            sobj.init = 0;
        end
        
        % Getter and setter function:
        %======================================================================
        %> @brief Function returning Model.states
        %>
        %> @param   obj     Model class object
        %> @retval  states  Table: Information on states of the model
        %======================================================================
        function x = get.states(obj)
            x = obj.hstates;
        end
        
        %======================================================================
        %> @brief Function returning Model.controls
        %>
        %> @param   obj         Model class object
        %> @retval  controls    Table: Information on the controls of the model
        %======================================================================
        function controls = get.controls(obj)
            controls = obj.hcontrols;
        end
        
        %======================================================================
        %> @brief Function returning Model.constraints
        %>
        %> @param   obj         Model class object
        %> @retval  constraints Table: Information on the constraints
        %======================================================================
        function constraints = get.constraints(obj)
            constraints = obj.hconstraints;
        end
        
        %======================================================================
        %> @brief Function returning Model.joints
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of joints
        %======================================================================
        function n = get.nJoints(obj)
            n = obj.hnJoints;
        end
        
        %======================================================================
        %> @brief Function returning Model.nMus
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of muscles
        %======================================================================
        function n = get.nMus(obj)
            n = obj.hnMus;
        end
        
        %======================================================================
        %> @brief Function returning Model.nTor
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of torque actuators
        %======================================================================
        function n = get.nTor(obj)
            n = obj.hnTor;
        end

        %======================================================================
        %> @brief Function returning Model.nCPs
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of contact points
        %======================================================================
        function n = get.nCPs(obj)
            n = obj.hnCPs;
        end
        
        %======================================================================
        %> @brief Function returning Model.nDofs
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of dofs
        %======================================================================
        function n = get.nDofs(obj)
            n = obj.hnDofs;
        end
        
        %======================================================================
        %> @brief Function returning Model.nSegements
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of segments
        %======================================================================
        function n = get.nSegments(obj)
            n = obj.hnSegments;
        end
        
        %======================================================================
        %> @brief Function returning Model.nStates
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of states
        %======================================================================
        function n = get.nStates(obj)
            n = obj.hnStates;
        end
        
        %======================================================================
        %> @brief Function returning Model.nControls
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of control variables
        %======================================================================
        function n = get.nControls(obj)
            n = obj.hnControls;
        end
        
        %======================================================================
        %> @brief Function returning Model.nConstraints
        %>
        %> @param   obj     Model class object
        %> @retval  n       Double: Number of constraints
        %======================================================================
        function n = get.nConstraints(obj)
            n = obj.hnConstraints;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxSymmetry
        %>
        %> @param   obj     Model class object
        %> @retval  sym     Double vector: Symmetry indices to map right to left
        %======================================================================
        function sym = get.idxSymmetry(obj)
            sym = obj.hidxSymmetry;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxForward
        %>
        %> @param   obj     Model class object
        %> @retval  ifor    Double vector: Indices of global dofs with translation in +x direction
        %======================================================================
        function ifor = get.idxForward(obj)
            ifor = obj.hidxForward;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxUpward
        %>
        %> @param   obj     Model class object
        %> @retval  isid    Double vector: Indices of global dofs with translation in +y direction
        %======================================================================
        function isid = get.idxUpward(obj)
            isid = obj.hidxUpward;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxSideward
        %>
        %> @param   obj     Model class object
        %> @retval  isid    Double vector: Indices of global dofs with translation in +z direction
        %======================================================================
        function isid = get.idxSideward(obj)
            isid = obj.hidxSideward;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxForwardAll
        %>
        %> @param   obj     Model class object
        %> @retval  ifor    Double vector: Indices of all dofs in +x direction
        %======================================================================
        function ifor = get.idxForwardAll(obj)
            ifor = obj.hidxForwardAll;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxSideward
        %>
        %> @param   obj     Model class object
        %> @retval  isid    Double vector: Indices of all dofs in +z direction
        %======================================================================
        function isid = get.idxSidewardAll(obj)
            isid = obj.hidxSidewardAll;
        end
        
        %======================================================================
        %> @brief Function returning Model.idxTorqueDof
        %>
        %> @param   obj     Model class object
        %> @retval  itor    Double vector: Indices with dofs having torque actuators
        %======================================================================
        function itor = get.idxTorqueDof(obj)
            itor = obj.hidxTorqueDof;
        end

        %======================================================================
        %> @brief Function setting Model.init
        %>
        %> @param   obj     Model class object
        %> @param  init   Boolean
        %======================================================================
        function set.init(obj, init)
            obj.init = init;
        end
        %> @endcond


        
        
    end
        
end


