%======================================================================
%> @file Gait3dTest.m
%> @brief Maltab m-file containing test class for class Gait3d
%>
%> @author Marlies Nitschke
%> @date December, 2017
%======================================================================

%======================================================================
%> @brief Test class to test the class Gait3d
%>
%> @details
%> Run single test selected by name: 
%> @code 
%> testCase = Gait3dTest;
%> testCase.model = Gait3d('gait3d.osim');
%> res = run(testCase, 'derivativetest_Dynamics');
%> @endcode
%>
%> Overview:
%> | Testname/Tag                    | Derivative | Freefall | Isometric | Memory | Speed | Neutral | getDynamics | getFkin | getGRF | getJointmoments | getMuscleforces | simuAccGyro | showStick |
%> |---------------------------------|------------|----------|-----------|--------|-------|---------|-------------|---------|--------|-----------------|-----------------|-------------|-----------|
%> | derivativetest_Dynamics         | x          |          |           |        |       |         | x           |         |        |                 |                 |             |           |
%> | derivativetest_SimuAccGyro_Acc  | x          |          |           |        |       |         |             |         |        |                 |                 | x           |           |
%> | derivativetest_SimuAccGyro_Gyro | x          |          |           |        |       |         |             |         |        |                 |                 | x           |           |
%> | derivativetest_Moments          | x          |          |           |        |       |         |             |         |        | x               |                 |             |           |
%> | derivativetest_GRF              | x          |          |           |        |       |         |             |         | x      |                 |                 |             |           |
%> | derivativetest_Fkin             | x          |          |           |        |       |         |             | x       |        |                 |                 |             |           |
%> | test_simulateFreefall           |            | x        |           |        |       |         | x           |         |        |                 |                 |             |           |
%> | test_showStickNeutral           |            |          |           |        |       | x       |             |         |        |                 |                 |             | x         |
%> | test_dynamics                   |            |          |           |        |       |         | x           |         | x      | x               |                 |             |           |
%> | test_memory                     |            |          |           | x      |       |         | x           |         |        |                 |                 |             |           |
%> | test_speedOfSimuAccGyro         |            |          |           |        | x     |         |             |         |        |                 |                 | x           |           |
%> | test_speedOfMex                 |            |          |           |        | x     |         | x           |         |        |                 |                 |             |           |
%> | test_s_hip_flexion_gyro         |            |          |           |        |       |         |             |         |        |                 |                 | x           |           |
%> | test_s_pelvis_rotation_gyro     |            |          |           |        |       |         |             |         |        |                 |                 | x           |           |
%> | test_s_hip_flexion_acc          |            |          |           |        |       |         |             |         | x      |                 |                 | x           |           |
%> | test_s_pelvis_rotation_acc      |            |          |           |        |       |         |             |         | x      |                 |                 | x           |           |
%> | test_s_static_upright_acc       |            |          |           |        |       | x       |             |         |        |                 |                 | x           |           |
%> | test_s_static_upright_gyro      |            |          |           |        |       | x       |             |         |        |                 |                 | x           |           |
%> | test_isometric                  |            |          | x         |        |       |         | x           |         |        | x               | x               |             |           |
%> | test_Fkin                       |            |          |           |        |       |         |             | x       |        |                 |                 |             | x         |
%> | test_grf                        |            |          |           |        |       |         |             |         | x      |                 |                 |             |           |
%>
%======================================================================
classdef Gait3dTest < ModelTest
    
    properties
        %> Struct: Options for Ipopt for test_grf()
        grfoptions
        %> Struct: Temporary storage for test_grf()
        grftmp
    end
    
    
    %> @name Tests: Gait3dTest
    %> @{
    methods (Test, TestTags = {'getGRF'})
        %=================================================================
        %> @brief Function for testing the ground reaction force model
        %>
        %> @details
        %> This test evaluates the contact model. To test the vertical
        %> force, the model is put in a series of positions with different
        %> hip height. The vertical ground reaction force at the heel is
        %> plotted as a function of hip height. The result (left side of
        %> Figure) shows the total force-deformation relationship of the
        %> contact model, i.e. deformation of the contact point relative
        %> to the foot, and the interaction between the foot and ground.
        %> To test the horizontal force, the model is put in three of these
        %> hip positions, and the horizontal velocity is varied. The
        %> resulting horizontal force is simulated and shown on the right
        %> side of the Figure.
        %>
        %> @image html Gait3dTest_test_grf.png width=800px
        %>
        %> @todo
        %> Add error checking
        %=================================================================
        function test_grf(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % vertical force, as function of vertical position at various vertical speeds
            y = 0.9501:0.0005:0.97;
            ny = size(y,2);
            vy = -1:0.5:1;
            nvy = size(vy,2);
            Fy = zeros(size(y,2),size(vy,2));
            for i = 1:nvy
                for j = 1:ny
                    [~,Fy(j,i)] = testCase.solvegrf(y(j),vy(i),0);
                end
            end
            figure('Name', 'test_grf');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            subplot(1,2,1);
            plot(y, Fy);
            legend([repmat('speed ',nvy,1) num2str(vy')]);
            xlabel('hip height (m)');
            ylabel('Fy (BW)');
            title('Vertical deformation test');
            
            % horizontal force, as function of horizontal speed at various vertical positions
            y = 0.95:0.005:0.97;
            ny = size(y,2);
            vx = -.2:0.01:.2;
            nvx = size(vx,2);
            Fx = zeros(nvx,ny);
            for i = 1:ny
                for j = 1:nvx
                    [Fx(j,i),~] = testCase.solvegrf(y(i),0,vx(j));
                end
            end
            subplot(1,2,2);
            plot(vx,Fx);
            legend([repmat('hip height ',ny,1) num2str(y')]);
            xlabel('horizontal speed (m/s)');
            ylabel('right Fx (BW)');
            title('Horizontal slip test');
            
        end
       
    end
    
    methods (Test, TestTags = {'Derivative','getFkin'})
        %-----------------------------------------------------------------
        %> @brief Function to test the derivatives of Gait3d.getFkin()
        %>
        %> @details
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Fkin
        %> dFKdq: Max. difference is:   0.0000000514 (  0.1216789450 vs.   0.1216788936)  at  'femur_r: R21'(idx:19) and 'pelvis_tilt'(idx:1) 
        %> dFKdotdq: Max. difference is:   0.0000001331 (  0.9619361218 vs.   0.9619359886)  at  'radius_r: R12'(idx:173) and 'elbow_flex_r'(idx:27) 
        %> @endcode 
        %>
        %> @todo 
        %> Add error checking
        %-----------------------------------------------------------------
        function derivativetest_Fkin(testCase)
            iCom = 0;
            % setup the Model object and data
            testCase = testCase.setup_Test_Random(iCom);
            
            % test the Jacobian of the forward kinematic model
            q = testCase.generateRandType('q',1);
            qd = testCase.generateRandType('qdot',1);
            [FK, dFKdq, dFKdotdq] = testCase.model.getFkin(q,qd);
            figure('Name', 'derivativetest_Fkin');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            subplot(1,2,1);
            spy(dFKdq,'.',5);xlim([1 size(dFKdq, 2)]); ylim([1 size(dFKdq, 1)]);
            title('sparsity pattern of dFK/dq');
            subplot(1,2,2);
            spy(dFKdotdq,'.',5);xlim([1 size(dFKdotdq, 2)]); ylim([1 size(dFKdotdq, 1)]);
            title('sparsity pattern of dFKdot/dq');
            
            hh = 1e-7;
            FKdot = dFKdq * qd;     % compute dfk/dt = dfk/qd * dq/dt
            dFKdq_num = zeros(size(dFKdq));
            dFKdotdq_num = zeros(size(dFKdotdq));
            for i = 1:testCase.model.nDofs
                qsave = q(i);
                q(i) = q(i) + hh;
                [FKnew, dFKnewdq] = testCase.model.getFkin(q);
                dFKdq_num(:,i) = (FKnew - FK)/hh;
                FKdotnew = dFKnewdq * qd;
                dFKdotdq_num(:,i) = (FKdotnew - FKdot)/hh;
                q(i) = qsave;
            end
            
            % compare all the Jacobian to their numerical estimtes
            fprintf('\n');
            disp('Derivatives for Fkin');
            qNames = testCase.model.states.name(strcmp(testCase.model.states.type,'q')); 
            delimiter(1:length(FK)) = {': '};
            p_R_Names = {'p1', 'p2', 'p3', 'R11', 'R12', 'R13', 'R21', 'R22', 'R23', 'R31', 'R32', 'R33'}; 
            p_R_Names = repmat(p_R_Names, 1, length(FK)/length(p_R_Names))';
            segmentNames = testCase.model.segments.Properties.RowNames(~strcmp(testCase.model.segments.Properties.RowNames, 'ground'));
            segmentNames = reshape(repmat(segmentNames, 1, length(FK)/length(segmentNames))', length(FK), 1);
            FKNames = strcat(segmentNames, delimiter', p_R_Names);
            testCase.matcompare(dFKdq   , dFKdq_num   , 'dFKdq'   , FKNames, qNames);
            testCase.matcompare(dFKdotdq, dFKdotdq_num, 'dFKdotdq', FKNames, qNames);
         
        end
        
    end
    
    methods (Test, TestTags = {'getFkin', 'showStick'})
        %=================================================================
        %> @brief Function for testing the forward kinematic model Gait3d.getFkin().
        %>
        %> @details
        %> This test sequentially animates the kinematic degrees of freedom. 
        %> First, the model is placed in a neutral position. Then each 
        %> degree of freedom is animated. After that the Jacobian dFk/dq 
        %> is computed i.e., the deritive of every segment position and 
        %> orientation with respect to each kinematic degree of freedom. 
        %> The result is compared to finite differences and the sparsity 
        %> of the jacobian is plotted.
        %> 
        %> @todo 
        %> Add error checking
        %>
        %> @todo Would be nice to do this interactively with slider GUI like in Opensim
        %=================================================================
        function test_Fkin(testCase)

            iCom = 0;
            % setup the Model object and data
            testCase = testCase.setup_Test_Random(iCom);
            
            % sequentially animates the kinematic degrees of freedom
            figure('Name', 'test_Fkin');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            for i = 1:testCase.model.nDofs
                x = testCase.model.states.xneutral;
                for q = 0:0.1:1
                    x(i) = 0.5*sin(2*pi*q);
                    clf;
                    testCase.model.showStick(x);
                    title(strrep(testCase.model.dofs.Properties.RowNames{i}, '_', '\_'));
                    axis([-.6 .6 -.6 .6 0 1.5])
                    view(3);    % use the default 3D viewpoint
                    drawnow;
                end
            end
            
        end
            
    end
    
    methods (Test, TestTags = {'Isometric', 'getDynamics', 'getJointmoments', 'getMuscleforces'})
        %=================================================================
        %> @brief Function for testing the isometric joint moments
        %>
        %> @image html Gait3dTest_test_isometricMuscles.png width=800px
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_isometric(testCase)
            % setup the Model object 
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
 
            dofs = testCase.model.dofs.Properties.RowNames;
            dofs = dofs(~strcmp(testCase.model.dofs.joint, 'ground_pelvis')); % Remove the ones to the global cooridnate system
            dofs = dofs(~endsWith(dofs, '_l')); % Remove left side;
            idxDofs = find(ismember(testCase.model.dofs.Properties.RowNames, dofs));
            
            figure('Name', 'test_isometric');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            plotrows = 3;
            plotcols = 5;

            for i = 1:numel(idxDofs)
                idof = idxDofs(i);
                [angles, pasmom, actmax, actmin] = testCase.moments(idof);
                subplot(plotrows,plotcols,i);
                angles = angles*180/pi;
                plot(angles,[pasmom actmax actmin]);
                xlabel('angle (deg)');
                ylabel('moment (Nm)');
                title(strrep(testCase.model.dofs.Properties.RowNames{idof}, '_','\_'));
                set(gca, 'XLim', [min(angles) max(angles)])
                drawnow
            end
            legend('passive','maximum','minimum'); % legend only on last subplot
        end
        
    end
    
    methods (Test, TestTags = {'Neutral', 'simuAccGyro'})
        %-----------------------------------------------------------------
        %> @brief Function to test the simulated acc signal for static standing
        %>
        %> @details
        %> Tests that the acceleration at calcn_r is equal to gravity
        %> during static standing.
        %-----------------------------------------------------------------
        function test_s_static_upright_acc(testCase)
            % setup the Model object 
            iCom = 4; % 3 acc  (one in each direction)
            testCase = testCase.setup_Test_Random(iCom);
            
            % adapt data structure
            %Not needed anymore after fixing setup_Test_Random():
            %{
            idxRightFoot = find(strcmp(testCase.model.segments.Properties.RowNames, 'calcn_r'));
            testCase.data.variables(1).info.segment = idxRightFoot;     
            testCase.data.variables(2).info.segment = idxRightFoot;     
            testCase.data.variables(3).info.segment = idxRightFoot;   
            %}
            
            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            qd = zeros(testCase.model.nDofs,1);
            qdd = zeros(testCase.model.nDofs,1);
            % simulate acc signal
            [s, dsdq, dsdqd, dsdqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % verify
            testCase.verifyEqual(norm(s), norm(testCase.model.gravity),'AbsTol', testCase.TOL_S, 'The norm of s is not correct.');
            testCase.verifyEqual(s(1), -testCase.model.gravity(1),'AbsTol', testCase.TOL_S, 's in x direction is not correct.');
            testCase.verifyEqual(s(2), -testCase.model.gravity(2),'AbsTol', testCase.TOL_S, 's in y direction is not correct.');
            testCase.verifyEqual(s(3), -testCase.model.gravity(3),'AbsTol', testCase.TOL_S, 's in z direction is not correct.');
        end
        
        %-----------------------------------------------------------------
        %> @brief Function to test the simulated gyro signal for static standing
        %>
        %> @details
        %> Tests that the gyropscope signal at calcn_r is equal to 0
        %> during static standing.
        %-----------------------------------------------------------------
        function test_s_static_upright_gyro(testCase)
            % setup the Model object 
            iCom = 5; % 3 gyro  (one in each direction)
            testCase = testCase.setup_Test_Random(iCom);
            
            % adapt data structure
            %Not needed anymore after fixing setup_Test_Random():
            %{
            idxRightFoot = find(strcmp(testCase.model.segments.Properties.RowNames, 'calcn_r'));
            testCase.data.variables(1).info.segment = idxRightFoot;     
            testCase.data.variables(2).info.segment = idxRightFoot;     
            testCase.data.variables(3).info.segment = idxRightFoot; 
            %}
            
            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            qd = zeros(testCase.model.nDofs,1);
            qdd = zeros(testCase.model.nDofs,1);
            % simulate gyro signal
            [s, dsdq, dsdqd, dsdqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % verify
            s_expected = [0, 0, 0];
            testCase.verifyEqual(norm(s), norm(s_expected),'AbsTol', testCase.TOL_S, 'The norm of s is not correct.');
            testCase.verifyEqual(s(1), s_expected(1),'AbsTol', testCase.TOL_S, 's in x direction is not correct.');
            testCase.verifyEqual(s(2), s_expected(2),'AbsTol', testCase.TOL_S, 's in y direction is not correct.');
            testCase.verifyEqual(s(3), s_expected(3),'AbsTol', testCase.TOL_S, 's in z direction is not correct.');
        end
        
    end
    
    methods (Test, TestTags = {'simuAccGyro', 'getFkin'})
        %-----------------------------------------------------------------
        %> @brief Function to test the simulated acc signal for hip flexion
        %> 
        %> @details
        %> Simulates the acceleration signal at calcn_r for hip flexion
        %> with 100 r/s. The expected acceleration is computed based on 
        %> the radius (computed from positions of segments) and known 
        %> rotational velocity.
        %-----------------------------------------------------------------
        function test_s_hip_flexion_acc(testCase)
            % setup the Model object 
            iCom = 4; % 3 acc  (one in each direction)
            testCase = testCase.setup_Test_Random(iCom);
            
            % adapt data structure
            idxRightFoot = find(strcmp(testCase.model.segments.Properties.RowNames, 'calcn_r'));
            %Not needed anymore after fixing setup_Test_Random():
            %{
            testCase.data.variables(1).info.segment = idxRightFoot;     
            testCase.data.variables(2).info.segment = idxRightFoot;     
            testCase.data.variables(3).info.segment = idxRightFoot;
            %}

            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            qd = zeros(testCase.model.nDofs,1);
            qd(strcmp(testCase.model.dofs.Properties.RowNames, 'hip_flexion_r')) = 100;       % 100 rad/s, we should get about 10^4 m/s2 upward acceleration
            qdd = zeros(testCase.model.nDofs,1);
            % simulate acc signal
            [s, dsdq, dsdqd, dsdqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % get true value
            FK = testCase.model.getFkin(q);
            dxRightFemur = find(strcmp(testCase.model.segments.Properties.RowNames, 'femur_r'));
            r = FK((dxRightFemur-2)*12 + (1:3)) - FK((idxRightFoot-2)*12 + (1:3)); 
            s_expected = -testCase.model.gravity + r'*(100)^2;% gravity + centripetal
            s_expected(3) = 0; % 0 because rotation is around z axis
            % verify
            testCase.verifyEqual(norm(s), norm(s_expected),'AbsTol', testCase.TOL_S, 'The norm of s is not correct.');
            testCase.verifyEqual(s(1), s_expected(1),'AbsTol', testCase.TOL_S, 's in x direction is not correct.');
            testCase.verifyEqual(s(2), s_expected(2),'AbsTol', testCase.TOL_S, 's in y direction is not correct.');
            testCase.verifyEqual(s(3), s_expected(3),'AbsTol', testCase.TOL_S, 's in z direction is not correct.');
        end
        
        %-----------------------------------------------------------------
        %> @brief Function to test the simulated acc signal for pelvis rotation    
        %> 
        %> @details
        %> Simulates the acceleration signal at calcn_r for pelvis rotation
        %> around 45 Deg and a hip flexion rate of 100 rad/s. The expected 
        %> acceleration is computed based on the radius (computed from positions
        %> of segments), the known hip rotation and known rotational velocity.
        %-----------------------------------------------------------------
        function test_s_pelvis_rotation_acc(testCase)
            % setup the Model object 
            iCom = 4; % 3 acc  (one in each direction)
            testCase = testCase.setup_Test_Random(iCom);
            
            % adapt data structure
            idxRightFoot = find(strcmp(testCase.model.segments.Properties.RowNames, 'calcn_r'));
            %Not needed anymore after fixing setup_Test_Random():
            %{
            testCase.data.variables(1).info.segment = idxRightFoot;     
            testCase.data.variables(2).info.segment = idxRightFoot;     
            testCase.data.variables(3).info.segment = idxRightFoot; 
            %}

            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            rot_y = 45;
            q(strcmp(testCase.model.dofs.Properties.RowNames, 'hip_rotation_r')) = rot_y/180*pi; % pelvis rotated around 45deg
            qd = zeros(testCase.model.nDofs,1);
            angVel_z = 100;
            qd(strcmp(testCase.model.dofs.Properties.RowNames, 'hip_flexion_r')) = angVel_z;       % 100 rad/s, (hip flexion is DOF 7)
            qdd = zeros(testCase.model.nDofs,1);
            % simulate acc signal
            [s, dsdq, dsdqd, dsdqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % get true value
            R_y = [cosd(rot_y),     0      ,sind(rot_y);...
                       0      ,     1      ,    0      ;...
                   -sind(rot_y),    0      ,cosd(rot_y)]';  
            FK = testCase.model.getFkin(q);
            dxRightFemur = find(strcmp(testCase.model.segments.Properties.RowNames, 'femur_r'));
            r = FK((dxRightFemur-2)*12 + (1:3)) - FK((idxRightFoot-2)*12 + (1:3)); 
            s_expected = -testCase.model.gravity + r'*(100)^2;% gravity + centripetal
            s_expected(3) = 0; % 0 because rotation is around z axis
            s_expected = R_y * s_expected';
            
            % verify
            testCase.verifyEqual(norm(s), norm(s_expected),'AbsTol', testCase.TOL_S, 'The norm of s is not correct.');
            testCase.verifyEqual(s(1), s_expected(1),'AbsTol', testCase.TOL_S, 's in x direction is not correct.');
            testCase.verifyEqual(s(2), s_expected(2),'AbsTol', testCase.TOL_S, 's in y direction is not correct.');
            testCase.verifyEqual(s(3), s_expected(3),'AbsTol', testCase.TOL_S, 's in z direction is not correct.');
        end
        
    end
    
    methods (Test, TestTags = {'simuAccGyro'})
        %-----------------------------------------------------------------
        %> @brief Function to test the simulated gyro signal for hip flexion
        %> 
        %> @details
        %> Simulates the gyropscope signal at calcn_r for hip flexion
        %> with 100 r/s. 
        %-----------------------------------------------------------------
        function test_s_hip_flexion_gyro(testCase)
            % setup the Model object 
            iCom = 5; % 3 gyro  (one in each direction)
            testCase = testCase.setup_Test_Random(iCom);
            
            % adapt data structure
            %Not needed anymore after fixing setup_Test_Random():
            %{
            idxRightFoot = find(strcmp(testCase.model.segments.Properties.RowNames, 'calcn_r'));
            testCase.data.variables(1).info.segment = idxRightFoot;     
            testCase.data.variables(2).info.segment = idxRightFoot;     
            testCase.data.variables(3).info.segment = idxRightFoot; 
            %}

            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            qd = zeros(testCase.model.nDofs,1);
            qd(strcmp(testCase.model.dofs.Properties.RowNames, 'hip_flexion_r')) = 100;       % 100 rad/s
            qdd = zeros(testCase.model.nDofs,1);
            % simulate gyro signal
            [s, dsdq, dsdqd, dsdqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % get true value
            s_expected = [0, 0, 100];
            % verify
            testCase.verifyEqual(norm(s), norm(s_expected),'AbsTol', testCase.TOL_S, 'The norm of s is not correct.');
            testCase.verifyEqual(s(1), s_expected(1),'AbsTol', testCase.TOL_S, 's in x direction is not correct.');
            testCase.verifyEqual(s(2), s_expected(2),'AbsTol', testCase.TOL_S, 's in y direction is not correct.');
            testCase.verifyEqual(s(3), s_expected(3),'AbsTol', testCase.TOL_S, 's in z direction is not correct.');
        end
        
        %-----------------------------------------------------------------
        %> @brief Function to test the simulated acc signal for pelvis rotation    
        %> 
        %> @details
        %> Simulates the acceleration signal at calcn_r for pelvis rotation
        %> around 45 Deg and a hip flexion rate of 100 rad/s. The expected 
        %> gyroscope signal is computed based on the known hip rotation and 
        %> known rotational velocity.
        %-----------------------------------------------------------------
        function test_s_pelvis_rotation_gyro(testCase)
             % setup the Model object 
            iCom = 5; % 3 gyro  (one in each direction)
            testCase = testCase.setup_Test_Random(iCom);
            
            % adapt data structure
            %Not needed anymore after fixing setup_Test_Random():
            %{
            idxRightFoot = find(strcmp(testCase.model.segments.Properties.RowNames, 'calcn_r'));
            testCase.data.variables(1).info.segment = idxRightFoot;
            testCase.data.variables(2).info.segment = idxRightFoot;
            testCase.data.variables(3).info.segment = idxRightFoot;
            %}

            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            rot_y = 45;
            q(strcmp(testCase.model.dofs.Properties.RowNames, 'hip_rotation_r')) = rot_y/180*pi; % pelvis rotated around 45deg
            qd = zeros(testCase.model.nDofs,1);
            angVel_z = 100;
            qd(strcmp(testCase.model.dofs.Properties.RowNames, 'hip_flexion_r')) = angVel_z;       % 100 rad/s, (hip flexion is DOF 7)
            qdd = zeros(testCase.model.nDofs,1);
            % get true value
            R_y = [cosd(rot_y),     0      ,sind(rot_y);...
                       0      ,     1      ,    0      ;...
                   -sind(rot_y),    0      ,cosd(rot_y)]';  
            s_expected = R_y * [0, 0, angVel_z]';
            % simulate gyro signal
            [s, dsdq, dsdqd, dsdqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % verify
            testCase.verifyEqual(norm(s), norm(s_expected),'AbsTol', testCase.TOL_S, 'The norm of s is not correct.');
            testCase.verifyEqual(s(1), s_expected(1),'AbsTol', testCase.TOL_S, 's in x direction is not correct.');
            testCase.verifyEqual(s(2), s_expected(2),'AbsTol', testCase.TOL_S, 's in y direction is not correct.');
            testCase.verifyEqual(s(3), s_expected(3),'AbsTol', testCase.TOL_S, 's in z direction is not correct.');
           
        end
        
    end
    %> @}
    
    methods (Access = protected)
        
        %-----------------------------------------------------------------
        %> @brief Function to setup the data
        %>
        %> @param  testCase
        %> @param  iCom
        %> @retval testCase
        %-----------------------------------------------------------------
        function testCase = setup_Test_Random(testCase, iCom)
            
            % generate random sensor data
            testCase.data = struct();
            if iCom == 1 % only acc
                % a random acc sensor for each segment
                sensorSegments = {'torso', 'femur_l', 'tibia_l', 'calcn_l', 'femur_r', 'tibia_r', 'calcn_r'}; % we only used this
                sizeArr=length(sensorSegments);
                nameArr=cell(sizeArr, 1);
                typeArr=cell(sizeArr, 1);
                unitArr=cell(sizeArr, 1);
                segment=cell(sizeArr, 1);
                direction=zeros(sizeArr,3);
                position=zeros(sizeArr,3);
                sim=cell(sizeArr, 1);
                for iSen = 1 : length(sensorSegments)
                    nameArr{iSen,1}=sensorSegments{iSen};
                    typeArr{iSen,1} = 'acc';
                    unitArr{iSen,1} = 'm/s^2';
                    segment{iSen,1} = sensorSegments{iSen};
                    axis = rand(3, 1);
                    direction(iSen, :) = axis/norm(axis);
                    position(iSen, :) = rand(3, 1);
                    sim{iSen,1} = [];
                end
                variableTable = table(nameArr, typeArr, unitArr, segment, position, direction, sim, 'VariableNames', ["name", "type", "unit", "segment", "position", "direction", "sim"]);
                testCase.data.variables = variableTable;

                testCase.data.nVars.acc = length(sensorSegments);
                testCase.data.idxVar.acc = 1 : length(sensorSegments);
                testCase.data.nVars.gyro = 0;
                testCase.data.idxVar.gyro = [];
                testCase.data.nVars.all = length(sensorSegments);
            end
            if iCom == 2 % only gyro
                % a random gyro sensor for each segment
                sensorSegments = {'torso', 'femur_l', 'tibia_l', 'calcn_l', 'femur_r', 'tibia_r', 'calcn_r'}; % we only used this
                sizeArr=length(sensorSegments);
                nameArr=cell(sizeArr, 1);
                typeArr=cell(sizeArr, 1);
                unitArr=cell(sizeArr, 1);
                segment=cell(sizeArr, 1);
                direction=zeros(sizeArr,3);
                position=zeros(sizeArr,3);
                sim=cell(sizeArr, 1);
                for iSen = 1 : length(sensorSegments)
                    nameArr{iSen,1}=sensorSegments{iSen};
                    typeArr{iSen,1} = 'gyro';
                    unitArr{iSen,1} = 'deg/s';
                    segment{iSen,1} = sensorSegments{iSen};
                    axis = rand(3, 1);
                    direction(iSen, :) = axis/norm(axis);
                    position(iSen, :) = rand(3, 1);
                    sim{iSen,1} = [];
                end
                variableTable = table(nameArr, typeArr, unitArr, segment, position, direction, sim, 'VariableNames', ["name", "type", "unit", "segment", "position", "direction", "sim"]);
                testCase.data.variables = variableTable;

                testCase.data.nVars.acc = 0;
                testCase.data.idxVar.acc = [];
                testCase.data.nVars.gyro = length(sensorSegments);
                testCase.data.idxVar.gyro = 1:length(sensorSegments);
                testCase.data.nVars.all = length(sensorSegments);
            end
            if iCom == 3
                % 21 acc and 21 gyro sensors
                sensorSegments = {'torso', 'femur_l', 'tibia_l', 'calcn_l', 'femur_r', 'tibia_r', 'calcn_r'};
                if iCom == 3 % acc in x and y direction and gyro in z direction
                    type = {'acc', 'acc', 'acc', 'gyro', 'gyro', 'gyro'};
                    direction = {[1; 0; 0], [0; 1; 0], [0; 0; 1], [1; 0; 0], [0; 1; 0], [0; 0; 1]};
                    
                    sizeArr=length(sensorSegments)*length(type);
                    nameArr=cell(sizeArr, 1);
                    typeArr=cell(sizeArr, 1);
                    unitArr=cell(sizeArr, 1);
                    segment=cell(sizeArr, 1);
                    directionArr=zeros(sizeArr,3);
                    positionArr=zeros(sizeArr,3);
                    sim=cell(sizeArr, 1);
                    index=1;

                    for iSig = 1 : length(type)
                        for iSen = 1 : length(sensorSegments)
                            nameArr{iSen + (iSig-1)*  length(sensorSegments), 1}=sensorSegments{iSen};
                            typeArr{iSen + (iSig-1)*  length(sensorSegments), 1}=type{iSig};
                            if index<=21
                                unitArr{iSen + (iSig-1)*  length(sensorSegments), 1}='m/s^2';
                            elseif index>21
                                unitArr{iSen + (iSig-1)*  length(sensorSegments), 1}='deg/s';
                            end
                            index=index+1;
                            segment{iSen + (iSig-1)*  length(sensorSegments), 1}=sensorSegments{iSen};
                            directionArr(iSen + (iSig-1)*  length(sensorSegments), :)=direction{iSig};
                            positionArr(iSen + (iSig-1)*  length(sensorSegments), :)=rand(3,1);
                            sim{iSen + (iSig-1)*  length(sensorSegments), 1} = [];
                        end
                    end
                    variableTable = table(nameArr, typeArr, unitArr, segment, positionArr, directionArr, sim, 'VariableNames', ["name", "type", "unit", "segment", "position", "direction", "sim"]);
                    testCase.data.variables = variableTable;

                    testCase.data.nVars.acc = 3*length(sensorSegments);
                    testCase.data.idxVar.acc = 1:3*length(sensorSegments);
                    testCase.data.nVars.gyro = 3*length(sensorSegments);
                    testCase.data.idxVar.gyro = 3*length(sensorSegments)+ (1 : 3*length(sensorSegments));
                    testCase.data.nVars.all = 6* length(sensorSegments);
                end                
            end
            if iCom == 4  || iCom == 5 % 3 acc or 3 gyro (one in each direction)
                if iCom == 4
                    type = 'acc';
                    unit = 'm/s^2';
                    noType = 'gyro';
                else
                    type = 'gyro';
                    unit = 'deg/s';
                    noType = 'acc';
                end
                sensorSegment='calcn_r';
                nameArr={};
                typeArr={};
                unitArr={};
                segment={};
                directionArr=[];
                positionArr=[];
                sim={};

                for iSen = 1 : 3
                    nameArr{iSen, 1}=sensorSegment;
                    typeArr{iSen, 1}=type;
                    unitArr{iSen, 1}=unit;
                    segment{iSen, 1}=sensorSegment;
                    if iSen == 1
                        directionArr(iSen,:)=[1,0,0];
                    elseif iSen == 2
                        directionArr(iSen,:)=[0,1,0];
                    else
                        directionArr(iSen,:)=[0,0,1];
                    end
                    positionArr(iSen, :)=zeros(3,1);
                    sim{iSen, 1}=[];
                end
                variableTable = table(nameArr, typeArr, unitArr, segment, positionArr, directionArr, sim, 'VariableNames', ["name", "type", "unit", "segment", "position", "direction", "sim"]);
                testCase.data.variables = variableTable;

                testCase.data.nVars.(type) = 3;
                testCase.data.idxVar.(type) = 1:3;
                testCase.data.nVars.(noType) = 0;
                testCase.data.idxVar.(noType) = [];
                testCase.data.nVars.all = 3;
            end
            
            
        end
        
        %======================================================================
        %> @brief Subfunction for GRF testing
        %>
        %> @details
        %> Function to solves grf by search for static equilibrium in contact variables
        %>
        %> @param testCase
        %> @param y
        %> @param vy
        %> @param vx
        %>
        %> @retval Fx
        %> @retval Fy
        %> @retval x
        %======================================================================
        function [Fx,Fy,x] = solvegrf(testCase, y, vy, vx)
            % solves the GRF when using a static contact model
            
            if isempty(testCase.grfoptions)
                testCase.grftmp.ixg = find(ismember(testCase.model.states.type, {'Fx', 'Fy', 'Fz', 'xc', 'yc', 'zc', 'xf', 'yf', 'zf'}));		% determine indices of grf state variables within x
                testCase.grftmp.ifg = find(ismember(testCase.model.constraints.type, 'CP'));		% determine indices of grf equations within f
                % determine the Jacobian structure for random inputs
                testCase.grftmp.y = rand;						% y of pelvis
                testCase.grftmp.vx = rand;                    % x velocity of pelvis
                testCase.grftmp.vy = rand;                    % y velocity of pelvis
                testCase.grftmp.xg = rand(size(testCase.grftmp.ixg));         % random values for the contact variables
                J = testCase.grfjacobian(testCase.grftmp.xg);
                testCase.grftmp.grfjacstruct = double(J~=0);
                
                % IPOPT settings
                testCase.grfoptions.lb = testCase.model.states.xmin(testCase.grftmp.ixg); % Lower bound on the variables.
                testCase.grfoptions.ub = testCase.model.states.xmax(testCase.grftmp.ixg); % Upper bound on the variables.
                testCase.grfoptions.cl = testCase.model.constraints.fmin(testCase.grftmp.ifg); % Lower bounds on the constraints.
                testCase.grfoptions.cu = testCase.model.constraints.fmax(testCase.grftmp.ifg); % Upper bounds on the constraints.
                testCase.grfoptions.ipopt.tol         = 1e-5;
                testCase.grftmp.grffuncs.objective         = @(x) 0;
                testCase.grftmp.grffuncs.gradient          = @(x) zeros(size(x));
                testCase.grftmp.grffuncs.constraints       = @testCase.grfconfun;
                testCase.grftmp.grffuncs.jacobian          = @testCase.grfjacobian;
                testCase.grftmp.grffuncs.jacobianstructure = @() testCase.grftmp.grfjacstruct;
                testCase.grfoptions.ipopt.hessian_approximation = 'limited-memory';
                testCase.grfoptions.ipopt.print_level = 0;
            end
            
            % find the equilibrium
            testCase.grftmp.y = y;
            testCase.grftmp.vy = vy;
            testCase.grftmp.vx = vx;
            [testCase.grftmp.xg, info] = ipopt(testCase.grftmp.xg,testCase.grftmp.grffuncs,testCase.grfoptions);
            if (info.status ~= 0)
                info.status
                %                 keyboard
            end
            x = testCase.model.states.xneutral;
            x(5) = testCase.grftmp.y;
            x(testCase.grftmp.ixg) = testCase.grftmp.xg;
            grf = testCase.model.getGRF(x);
            Fx = grf(1);						% Fx on right foot
            Fy = grf(2);						% Fy on right foot
            % [y vy vx Fx Fy]
            
        end
        
        % ----------------------------------------------------------------
        %> @brief Subfunction for GRF testing
        %>
        %> @param  testCase
        %> @param  X
        %> @retval c
        % ----------------------------------------------------------------
        function [c] = grfconfun(testCase, X)
            x = testCase.model.states.xneutral;
            xd = zeros(testCase.model.nStates,1);
            u  = zeros(testCase.model.nControls,1);
            x(5) = testCase.grftmp.y;
            xd(testCase.grftmp.ixg(4:9:end)) = testCase.grftmp.vx;			% put x velocity in all contact points
            xd(testCase.grftmp.ixg(5:9:end)) = testCase.grftmp.vy;			% put y velocity in all contact points
            x(testCase.grftmp.ixg) = X;
            c = testCase.model.getDynamics(x,xd,u);
            c = c(testCase.grftmp.ifg);							% only use those that are associated with grf equations
        end
        
        % ----------------------------------------------------------------
        %> @brief Subfunction for GRF testing
        %>
        %> @param  testCase
        %> @param  X
        %> @retval J
        % ----------------------------------------------------------------
        function [J] = grfjacobian(testCase,X)
            x = testCase.model.states.xneutral;
            xd = zeros(testCase.model.nStates,1);
            u  = zeros(testCase.model.nControls,1);
            x(5) = testCase.grftmp.y;
            xd(testCase.grftmp.ixg(4:9:end)) = testCase.grftmp.vx;			% put x velocity in all contact points
            xd(testCase.grftmp.ixg(5:9:end)) = testCase.grftmp.vy;			% put y velocity in all contact points
            x(testCase.grftmp.ixg) = X;
            [~, dcdx] = testCase.model.getDynamics(x,xd,u);
            J = dcdx(testCase.grftmp.ixg, testCase.grftmp.ifg)';
        end
        
        %-----------------------------------------------------------------
        %> @brief Function to generate active and passive joint moment curves
        %>
        %> @param  testCase
        %> @param  idof
        %> @retval angles
        %> @retval pasmom
        %> @retval actmax
        %> @retval actmin
        %-----------------------------------------------------------------
        function [angles, pasmom, actmax, actmin] = moments(testCase,idof)
            % generates active and passive joint moment curves

            Npts = 100;
            qmin = testCase.model.dofs{idof,'range_passiveMoment'}(1)-10/180*pi;
            qmax = testCase.model.dofs{idof,'range_passiveMoment'}(2)+10/180*pi;
            angles = linspace(qmin, qmax, Npts)';
            pasmom = zeros(size(angles));
            actmin = zeros(size(angles));
            actmax = zeros(size(angles));

            % passive muscles
            x = testCase.model.states.xneutral;
            x(testCase.model.extractState('s')) = 0.5;   % better initial guess for contraction state
            for i = 1:numel(angles)
                x(idof) = angles(i);
                [pasmom(i),x] = testCase.maxmom(x,idof,0);
            end

            % active muscles
            % we can use common initial guess for max and min, these are different
            % muscles
            x(testCase.model.extractState('s')) = 0.5;   % better initial guess for contraction state
            for i = 1:numel(angles)
                x(idof) = angles(i);
                [actmax(i),x] = testCase.maxmom(x,idof,1);
            end
            x(testCase.model.extractState('s')) = 0.5;   % better initial guess for contraction state
            for i = 1:numel(angles)
                x(idof) = angles(i);
                [actmin(i),x] = testCase.maxmom(x,idof,-1);
            end

        end
        
        %-----------------------------------------------------------------
        %> @brief Function to find the maximal moment for a joint 
        %>
        %> @details
        %> Function to find the maximal moment for a joint with the model 
        %> in posture given by x while muscles are required to be in
        %> equilibrium contraction state.
        %>
        %> @param   testCase
        %> @param	x			posture of the model
        %> @param	idof    	which joint moment we want to maximize
        %> @param   direction	+1 for max moment, -1 for min moment
        %>
        %> @retval	M			the maximum moment
        %> @retval	x			the system state that produces the max moment
        %>				        (can be used as initial guess for next function call)
        %-----------------------------------------------------------------
        function [M,x] = maxmom(testCase, x, idof, direction)
            
            % set activation states of all muscles to zero
            idxActivation = testCase.model.extractState('a');
            x(idxActivation) = 0;

            % if active moment is required, find moment arms for muscles that will
            % increase moment(idof) in the specified direction
            % set activation state to 1.0 for those muscles
            if direction ~= 0 && testCase.model.nMus > 0
                [~,~,MA] = testCase.model.getMuscleforces(x);
                x(idxActivation(direction*MA(:,idof) > 0)) = 1.0;
            end

            % find the equilibrium contraction state for all muscles
            % use Newton's method to solve x from f(x,0,u) = 0 (u does not matter) 
            u = zeros(testCase.model.nControls,1);
            xd = zeros(size(x));	
              
            idxLE = testCase.model.extractState('s'); % indices for contraction states in x 
            idxContDym = find(strcmp(testCase.model.constraints.type, 'muscles') & ...
                strcmp(testCase.model.constraints.equation, 'contraction dynamics: f = Fsee - Fce - Fpee'));  % and equations in f
            maxiterations = 30;
            tolerance = 1e-5;
            for i = 1:maxiterations
                [f, dfdx] = testCase.model.getDynamics(x,xd,u);
                F = f(idxLE);
                J = dfdx(idxLE,idxContDym);
                x(idxLE) = x(idxLE) - J\F;
                if norm(F) < tolerance
                    break
                end
            end
            if i == maxiterations
                warning(['maxmom: could not find equilibrium contraction states: ', num2str(x(idof)/pi*180) ' ; ' num2str(idof), ' ; ' num2str(direction) '\n']);
            end

            % find the joint moment with the model in this state
            allmoments = testCase.model.getJointmoments(x, u);
            
            M = allmoments(idof);
        end	
        
    end
    
    
end