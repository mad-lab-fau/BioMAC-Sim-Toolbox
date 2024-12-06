%======================================================================
%> @file Gait2dcTest.m
%> @brief Maltab m-file containing test class for class Gait2dc
%>
%> @author Marlies Nitschke
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Test class to test the class Gait2dc
%>
%> @details
%> Run single test selected by name: 
%> @code 
%> testCase = Gait2dcTest;
%> testCase.model = Gait2dc('gait2dc_par.xls');
%> res = run(testCase, 'derivativetest_Dynamics');
%> @endcode
%>
%> Overview:
%> | Testname/Tag                    | Apparel | Derivative | Freefall | Isometric | Memory | Speed | Neutral | getDynamics | getGRF | getJointmoments | simuAccGyro | showStick |
%> |---------------------------------|---------|------------|----------|-----------|--------|-------|---------|-------------|--------|-----------------|-------------|-----------|
%> | derivativetest_Dynamics         |         | x          |          |           |        |       |         | x           |        |                 |             |           |
%> | derivativetest_SimuAccGyro_Acc  |         | x          |          |           |        |       |         |             |        |                 | x           |           |
%> | derivativetest_SimuAccGyro_Gyro |         | x          |          |           |        |       |         |             |        |                 | x           |           |
%> | derivativetest_Moments          |         | x          |          |           |        |       |         |             |        | x               |             |           |
%> | derivativetest_GRF              |         | x          |          |           |        |       |         |             | x      |                 |             |           |
%> | test_simulateFreefall           |         |            | x        |           |        |       |         | x           |        |                 |             |           |
%> | test_showStickNeutral           |         |            |          |           |        |       | x       |             |        |                 |             | x         |
%> | test_dynamics                   |         |            |          |           |        |       |         | x           | x      | x               |             |           |
%> | test_memory                     |         |            |          |           | x      |       |         | x           |        |                 |             |           |
%> | test_speedOfSimuAccGyro         |         |            |          |           |        | x     |         |             |        |                 | x           |           |
%> | test_speedOfMex                 |         |            |          |           |        | x     |         | x           |        |                 |             |           |
%> | test_s_hip_flexion_gyro         |         |            |          |           |        |       |         |             |        |                 | x           |           |
%> | test_s_static_upright_acc       |         |            |          |           |        |       | x       |             |        |                 | x           |           |
%> | test_s_static_upright_gyro      |         |            |          |           |        |       | x       |             |        |                 | x           |           |
%> | test_dynamicGRF                 |         |            |          |           |        |       |         |             | x      |                 |             |           |
%> | test_isokineticMuscles          |         |            |          |           |        |       | x       | x           |        | x               |             |           |
%> | test_isometricMuscles           |         |            |          | x         |        |       |         | x           |        | x               |             |           |
%> | test_apparel                    | x       |            |          |           |        |       |         | x           |        | x               |             |           |
%> | test_grf                        |         |            |          |           |        |       |         |             | x      |                 |             |           |
%>
%======================================================================
classdef Gait2dcTest < ModelTest
       
    %> @name Tests: Gait2dcTest
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
        %> @image html Gait2dcTest_test_grf.png width=800px
        %> 
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_grf(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            xneutral = testCase.model.states.xneutral;
            
            % vertical force, as function of vertical position at various vertical speeds
            y = -0.02:0.0001:0.02;
            ny = size(y,2);
            vy = 0;
            nvy = size(vy,2);
            Fy = zeros(size(y,2),size(vy,2));
            for i = 1:nvy
                for j = 1:ny
                    GRF = testCase.solvegrf(xneutral(2) - xneutral(52) + y(j),vy(i),0,0,0);
                    Fy(j,i) = GRF(2);		% vertical GRF of right foot (heel)
                end
            end
            figure('Name', 'test_grf');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            subplot(1,2,1);
            plot(xneutral(2) - xneutral(52) + y,Fy);
            legend([repmat('speed ',nvy,1) num2str(vy')]);
            xlabel('hip height (m)');
            ylabel('right Fy (BW)');
            title('Vertical deformation test');
            
            % horizontal force, as function of horizontal speed at various vertical positions
            vx = -.2:0.001:.2;
            nvx = size(vx,2);
            y = -0.02:0.01:0;
            ny = size(y,2);
            Fx = zeros(nvx,ny);
            for i = 1:ny
                for j = 1:nvx
                    GRF = testCase.solvegrf(xneutral(2) - xneutral(52) + y(i),0,0,vx(j),vx(j));
                    Fx(j,i) = GRF(1);		% horizontal GRF of right foot (heel)
                end
            end
            subplot(1,2,2);
            plot(vx,Fx);
            legend([repmat('hip height ',ny,1) num2str((xneutral(2) - xneutral(52) + y)')]);
            xlabel('horizontal speed (m/s)');
            ylabel('right Fx (BW)');
            title('Horizontal slip test');
            
        end
        
    end
    
    methods (Test, TestTags = {'Apparel','getDynamics', 'getJointmoments'})
        %=================================================================
        %> @brief Function for testing the model using apparel
        %>
        %> @details
        %> Tests additional joint moments to model apparel.
        %> (See "2D Model with Contact (gait2dc) Reference Manual for Version 1.0")
        %>
        %> @image html Gait2dcTest_test_apparel.png width=800px
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_apparel(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % add apparel model to the parameters matrix
            k = 10000;					% stiffness of structure (N/m)
            L0 = -0.10;					% elongation of the structure at full extension
            dH = -0.15;					% moment arm at hip (m), negative because structure is posterior to hip
            dK = 0.10;					% moment arm at knee (m), positive because structure is anterior to knee
            polynomial = [ ...
                -k * dH * L0		1	0	0	0	0	0; ...
                -k * dK * L0		0	1	0	0	0	0; ...
                0.5 * k * dH^2		2	0	0	0	0	0; ...
                0.5 * k * dK^2		0	2	0	0	0	0; ...
                k * dH*dK			1	1	0	0	0	0; ...
                -k * dH * L0		0	0	0	1	0	0; ...
                -k * dK * L0		0	0	0	0	1	0; ...
                0.5 * k * dH^2		0	0	0	2	0	0; ...
                0.5 * k * dK^2		0	0	0	0	2	0; ...
                k * dH*dK			0	0	0	1	1	0; ...
                ];
            polynomial = [polynomial;polynomial]; % Copy for other foot
            evalc('testCase.model.strainEnergyTerms = polynomial;');
            
            figure('Name', 'test_apparel');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            % hip moment-angle curves, at 30-deg intervals in knee angle
            testCase.isometric_curves(1,-30:5:90,2,-120:30:0);
            % knee moment-angle curves, at 30-deg intervals in hip angle
            testCase.isometric_curves(2,-120:5:0,1,-30:30:90);
            
        end
        
    end
    
    methods (Test, TestTags = {'Isometric', 'getDynamics', 'getJointmoments'})
        %=================================================================
        %> @brief Function for testing the muscles strength
        %>
        %> @details
        %> This puts the model into various joint angle combinations and activates 
        %> the muscle set three times for each joint. First, activating the muscles 
        %> that contribute to positive moment. Second, activating the muscles that 
        %> contribute to a negative moment. Third, no activation, to obtain the 
        %> passive moment.
        %> (See "2D Model with Contact (gait2dc) Reference Manual for Version 1.0")
        %>
        %> @image html Gait2dcTest_test_isometricMuscles.png width=800px
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_isometricMuscles(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % moment-angular velocity curves      
            figure('Name', 'test_isometricMuscles');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            % hip moment-angle curves, at 30-deg intervals in knee angle
            testCase.isometric_curves(1,-30:5:90,2,-120:30:0);
            % knee moment-angle curves, at 30-deg intervals in hip angle
            testCase.isometric_curves(2,-120:5:0,1,-30:30:90);
            % ankle moment-angle curves, at 30-deg intervals in knee angle
            testCase.isometric_curves(3,-60:5:60,2,-120:30:0);
        end
        
    end
        
    methods (Test, TestTags = {'Neutral', 'getDynamics', 'getJointmoments'})
        %=================================================================
        %> @brief Function for testing the muscles strength
        %>
        %> @details
        %> For the isokinetic test, the model is placed in the neutral position 
        %> (hips and ankles at zero degrees, knees at 5 degrees flexion). At each 
        %> joint, angular velocity is then varied from -1000 to +1000 deg/s. At 
        %> each angular velocity, the state of the system is found that produces 
        %> steady state muscle forces (dF/dt = 0). The muscle forces are then 
        %> converted to joint moments and plotted. This is done separately for the 
        %> flexors and extensors.
        %> (See "2D Model with Contact (gait2dc) Reference Manual for Version 1.0")
        %>
        %> @image html Gait2dcTest_test_isokineticMuscles.png width=800px
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_isokineticMuscles(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % moment-angular velocity curves      
            figure('Name', 'test_isokineticMuscles');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            % moment-angular velocity curves
            testCase.isokinetic_curves(1,-1000:50:1000);
            testCase.isokinetic_curves(2,-1000:50:1000);
            testCase.isokinetic_curves(3,-1000:50:1000);
        end
        
    end
    
    methods (Test, TestTags = {'getGRF'})
        %=================================================================
        %> @brief Function for testing the GRF
        %>
        %> @details
        %> The foot will be moved along a prescribed trajectory, which is circular 
        %> with diminishing radius. The speed is one revolution per second.
        %> (See "2D Model with Contact (gait2dc) Reference Manual for Version 1.0")
        %>
        %> @code{.unparsed}
        %> Test dynamic GRFs
        %> Number of time steps:           1100
        %> Number of function evaluations: 23955
        %> @endcode
        %> @image html Gait2dcTest_test_dynamicGRF.png width=800px
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_dynamicGRF(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % set solver parameters for IMstepgrf
            tend = 10.0;		% duration of simulation
            h = 0.01;
            times = 0:h:tend;
            options.tol = 1e-8;
            options.maxiterations = 100;
            
            Xamplitude = 0.012;
            Yamplitude = 0.012;
            period = 1.0;	% 1 s second period for the sine
            
            % Get indices
            idxPelvis_ty  = find(strcmp(testCase.model.states.type, 'q') & ...
                                 strcmp(testCase.model.states.name, 'pelvis_ty'));
            idxHeelrCP_y  = find(strcmp(testCase.model.states.type, 'yc') & ...
                                 strcmp(testCase.model.states.name, 'heel_r'));
            idxAllCPs     = find(strcmp(testCase.model.states.type, 'xc') | ...
                                 strcmp(testCase.model.states.type, 'yc') | ...
                                 strcmp(testCase.model.states.type, 'Fx') | ...
                                 strcmp(testCase.model.states.type, 'Fy'));
            idxAllCPs_y   = find(strcmp(testCase.model.states.type, 'yc'));
            idxAllCPs_y_inCPs = find(ismember(idxAllCPs, idxAllCPs_y));
            
            xneutral = testCase.model.states.xneutral;
            setpointx2 = xneutral(idxPelvis_ty) - xneutral(idxHeelrCP_y) - 0.005;	% setpoint is where contact points just touch the ground
            xinit_CPs = xneutral(idxAllCPs);		% initial conditions for contact variables
            xinit_CPs(idxAllCPs_y_inCPs) = xinit_CPs(idxAllCPs_y_inCPs) - xneutral(idxHeelrCP_y);		% move initial cond for contact points down also, so foot is not deformed at t=0
            
            
            x = repmat(xneutral,1,numel(times)); % make skeleton trajectory
            decayfactor = exp(-times/5);
            x(1,:) = Xamplitude*decayfactor.*sin(2*pi*times/period);
            x(2,:) = Yamplitude*decayfactor.*cos(2*pi*times/period) + setpointx2;		% use cos to start at the top!
            x = [repmat(x(:,1),1,100) x];		% add a static period at the beginning
            times = [(-100:-1)*h times];
            xinit_CPs(idxAllCPs_y_inCPs) = xinit_CPs(idxAllCPs_y_inCPs) + Yamplitude;
            [tout, xout, info, neval] = testCase.stepgrf(x,xinit_CPs,times,options);
            fprintf('\n');
            disp('Test dynamic GRFs');
            fprintf('Number of time steps:           %d\n', numel(tout)-1);
            fprintf('Number of function evaluations: %d\n', neval);
            x(idxAllCPs,:) = xout;					% put the contact variables in their right place in x
            
            % plot
            figure('Name', 'test_dynamicGRF:Animation');
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            
            for iTime=1:10:size(x,2)-1
                testCase.model.showStick(x(:,iTime));
                hold on
                plot(x(idxAllCPs(idxAllCPs_y_inCPs)-1,iTime+1),x(idxAllCPs(idxAllCPs_y_inCPs),iTime+1),'o','Markersize',10);
                hold off
                axis('equal');
                axis([-0.1 0.2 -0.05 0.1])
                pause(0.04);
            end
            Fxy = zeros(numel(times)-1,2);
            pxy = zeros(numel(times)-1,2);
            for iTime=1:size(x,2)-1
                grf = testCase.model.getGRF(x(:,iTime+1))';
                Fxy(iTime,:) = grf(1:2)';
                pxy(iTime,:) = x(1:2,iTime+1)';
            end
            
            % plot curves
            figure('Name', 'test_dynamicGRF');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            subplot(2,2,1)
            plot(pxy(:,1),pxy(:,2));
            axis('equal');
            xlabel('x_{hip} (m)');
            ylabel('y_{hip} (m)');
            title('prescribed input motion');
            
            subplot(2,2,2)
            plot(Fxy(:,1),Fxy(:,2));
            axis('equal');
            xlabel('Fx (BW)');
            ylabel('Fy (BW)');
            title('simulated forces');
            
            subplot(2,2,3)
            plot(pxy(:,1),Fxy(:,1));
            xlabel('x_{hip} (m)');
            ylabel('Fx (BW)');
            title('horizontal deformation & force');
            
            subplot(2,2,4)
            plot(pxy(:,2),Fxy(:,2));
            xlabel('y_{hip} (m)');
            ylabel('Fy (BW)');
            title('vertical deformation & force');
            
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
            testCase.data.variables.name{1} = 'foot_r';     
            testCase.data.variables.name{2} = 'foot_r';     
            testCase.data.variables.name{3} = 'foot_r';   
            
            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            qd = zeros(testCase.model.nDofs,1);
            qdd = zeros(testCase.model.nDofs,1);
            % simulate acc signal
            s = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % verify
            testCase.verifyEqual(norm(s), norm(testCase.model.gravity),'AbsTol', testCase.TOL_S, 'The norm of s is not correct.');
            testCase.verifyEqual(s(1), 0                     ,'AbsTol', testCase.TOL_S, 's in x direction is not correct.');
            testCase.verifyEqual(s(2), testCase.model.gravity,'AbsTol', testCase.TOL_S, 's in y direction is not correct.');
            testCase.verifyEqual(s(3), 0                     ,'AbsTol', testCase.TOL_S, 's in z direction is not correct.');
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
            testCase.data.variables.name{1} = 'foot_r';     
            testCase.data.variables.name{2} = 'foot_r';     
            testCase.data.variables.name{3} = 'foot_r'; 
            
            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            qd = zeros(testCase.model.nDofs,1);
            qdd = zeros(testCase.model.nDofs,1);
            % simulate acc signal
            s = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % verify
            s_expected = [0, 0, 0];
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
            testCase.data.variables.name{1} = 'foot_r';     
            testCase.data.variables.name{2} = 'foot_r';     
            testCase.data.variables.name{3} = 'foot_r';
            testCase.data.variables.segment{1} = 'foot_r';     
            testCase.data.variables.segment{2} = 'foot_r';     
            testCase.data.variables.segment{3} = 'foot_r';
            
            % generate state
            q = testCase.model.states.xneutral(1:testCase.model.nDofs);
            qd = zeros(testCase.model.nDofs,1);
            qd(strcmp(testCase.model.dofs.Properties.RowNames, 'hip_flexion_r')) = 100;       % 100 rad/s
            qdd = zeros(testCase.model.nDofs,1);
            % simulate acc signal
            s = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % get true value
            s_expected = [0, 0, 100];
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
            nType = 1;
            if iCom == 3
                nVars = 3*testCase.model.nSegments;
            elseif  iCom == 4  || iCom == 5
                nVars = 3;
            else
                nVars = testCase.model.nSegments;
            end
            type = cell(nVars,1);
            name = cell(nVars,1);
            segment = cell(nVars,1);
            position = zeros(nVars,3);
            direction = zeros(nVars,3);
            testCase.data.variables = table(type,name,segment,position,direction);
            %testCase.data.variables.Properties.VariableNames = {'type','name','position','direction'};
            if iCom == 1 % only acc
                % a random acc sensor for each segment
                for iSen = 1 : testCase.model.nSegments
                    testCase.data.variables.type{iSen} = 'acc';
                    testCase.data.variables.name{iSen} = testCase.model.segments.Properties.RowNames{iSen};
                    testCase.data.variables.segment{iSen} = testCase.model.segments.Properties.RowNames{iSen};
                    axis = rand(3, 1);
                    testCase.data.variables.direction(iSen,:) = axis/norm(axis);
                    testCase.data.variables.position(iSen,:) = rand(3, 1);
                end
%                 testCase.data.nVars.acc = testCase.model.nSegments;
%                 testCase.data.idxVar.acc = 1 : testCase.model.nSegments;
%                 testCase.data.nVars.gyro = 0;
%                 testCase.data.idxVar.gyro = [];
%                 testCase.data.nVars.all = testCase.model.nSegments;
            end
            if iCom == 2  % only gyro
                % a random gyro sensor for each segment
                for iSen = 1 : testCase.model.nSegments
                    testCase.data.variables.type{iSen} = 'gyro';
                    testCase.data.variables.name{iSen} = testCase.model.segments.Properties.RowNames{iSen};
                    testCase.data.variables.segment{iSen} = testCase.model.segments.Properties.RowNames{iSen};
                    axis = rand(3, 1);
                    testCase.data.variables.direction(iSen,:) = axis/norm(axis);
                    testCase.data.variables.position(iSen,:) = rand(3, 1);
                end
%                 testCase.data.nVars.acc = 0;
%                 testCase.data.idxVar.acc = [];
%                 testCase.data.nVars.gyro = testCase.model.nSegments;
%                 testCase.data.idxVar.gyro = 1 : testCase.model.nSegments;
%                 testCase.data.nVars.all = testCase.model.nSegments;
            end
            if iCom == 3 % acc in x and y direction and gyro in z direction
                type = {'acc', 'acc', 'gyro'};
                direction = {[1; 0; 0], [0; 1; 0], [0; 0; 1]};
                for iSig = 1 : length(type)
                    for iSen = 1 : testCase.model.nSegments
                       testCase.data.variables.type{(iSen + (iSig-1)*  testCase.model.nSegments)} = type{iSig};
                       testCase.data.variables.name{(iSen + (iSig-1)*  testCase.model.nSegments)} = testCase.model.segments.Properties.RowNames{iSen};
                       testCase.data.variables.segment{(iSen + (iSig-1)*  testCase.model.nSegments)} = testCase.model.segments.Properties.RowNames{iSen};
                       testCase.data.variables.direction((iSen + (iSig-1)*  testCase.model.nSegments),:) = direction{iSig};
                       testCase.data.variables.position((iSen + (iSig-1)*  testCase.model.nSegments),:) = rand(3, 1);
                    end
                end
%                 testCase.data.nVars.acc = 2*testCase.model.nSegments;
%                 testCase.data.idxVar.acc = 1:2*testCase.model.nSegments;
%                 testCase.data.nVars.gyro = testCase.model.nSegments;
%                 testCase.data.idxVar.gyro = 2*testCase.model.nSegments + (1 : testCase.model.nSegments);
%                 testCase.data.nVars.all = 3* testCase.model.nSegments;             
            end
            if iCom == 4  || iCom == 5 % 3 acc or 3 gyro (one in each direction)
                if iCom == 4
                    type = 'acc';
                    noType = 'gyro';
                else
                    type = 'gyro';
                    noType = 'acc';
                end
                for iSen = 1 : 3
                    testCase.data.variables.type{iSen} = type;
                    testCase.data.variables.name{iSen} = 'pelvis';
                    testCase.data.variables.segment{iSen} = 'pelvis';
                    testCase.data.variables.direction(iSen,:) = zeros(3, 1);
                    testCase.data.variables.direction(iSen,iSen) = 1;
                    testCase.data.variables.position(iSen,:) = zeros(3, 1);
                end
%                 testCase.data.nVars.(type) = 3;
%                 testCase.data.idxVar.(type) = 1:3;
%                 testCase.data.nVars.(noType) = 0;
%                 testCase.data.idxVar.(noType) = [];
%                 testCase.data.nVars.all = 3;
            end
            
        end
        
        %======================================================================
        %> @brief Function to solves grf by search for static equilibrium in contact variables
        %>
        %> @param testCase
        %> @param y
        %> @param vy
        %> @param xx
        %> @param vx
        %> @param vxc
        %>
        %> @retval grf Ground reaction force
        %======================================================================
        function [grf] = solvegrf(testCase,y,vy,xx,vx,vxc)
            
            x = testCase.model.states.xneutral;% put model in upright position
            x(1) = xx;
            x(2) = y;						% hip at height y
            x(10) = vx;						% hip horizontal velocity
            x(11) = vy;						% hip vertical velocity
            
            xdot = zeros(testCase.model.nStates,1);
            ix = 2*testCase.model.nDofs + 2*testCase.model.nMus + 4*(1:testCase.model.nCPs)-3;		% indices of X coordinates of contact point within state vector x
            xdot(ix) = vxc;									% give contact points the horizontal speed vxc
            xdot(ix+1) = vy;								% and vertical speed vy
            
            u = zeros(testCase.model.nControls,1);	            % no stim, we don't care about activation dynamics here
            icontact = 2*testCase.model.nDofs+2*testCase.model.nMus+(1:4*testCase.model.nCPs);	% index of all contact variables within x
            
            options.ftol = 1e-5;						% tolerance of iterative Newton solver
            options.maxiterations = 100;
            if isempty(testCase.grf_xc)
                testCase.grf_xc = x(icontact);			% initial guess
            end
            [testCase.grf_xc, f, info] = testCase.fsolve1(@grf_equilibrium, testCase.grf_xc, options);
            if (info ~= 0)
                grf = NaN(6,1);						% no solution found
            else
                grf = testCase.model.getGRF(x);		% determine the grf at this state of the system
            end
            
            return
            
            %======================================================================
            %> @brief Helperfunction to test grf
            %>
            %> @param xcontact
            %>
            %> @retval F
            %> @retval J
            %=====================================================================
            function [F,J] = grf_equilibrium(xcontact)
                x(icontact) = xcontact;
                [f,dfdx] = testCase.model.getDynamics(x,xdot,u);
                F = f(icontact);
                J = dfdx(icontact,icontact)';
            end
        end
        
        
        %======================================================================
        %> @brief Function to plot isometric curves 
        %>
        %> @details
        %> Produces isometric moment-angle curves for a joint, one for each value 
        %> in a range of angles in second joint.
        %> Third joint angle is kept zero.
        %>
        %> @param testCase
        %> @param joint1
        %> @param range1
        %> @param joint2
        %> @param range2
        %======================================================================
        function isometric_curves(testCase, joint1, range1, joint2, range2)
            
            joints = testCase.model.joints.Properties.RowNames;
            angles = zeros(testCase.model.nJoints,1)+1e-6;
            angvel = zeros(testCase.model.nJoints,1);
            pascurves = [];
            poscurves = [];
            negcurves = [];
            legends = {};
            for angle2 = range2+1e-6
                angles(joint2) = angle2;
                pasmoments = [];
                posmoments = [];
                negmoments = [];
                for angle1 = range1+1e-6
                    angles(joint1) = angle1;
                    pasmoments = [pasmoments testCase.maxmoment(joint1, angles, angvel, 0)];
                    posmoments = [posmoments testCase.maxmoment(joint1, angles, angvel, 1)];
                    negmoments = [negmoments testCase.maxmoment(joint1, angles, angvel, -1)];
                end
                pascurves = [pascurves  pasmoments'];
                poscurves = [poscurves  posmoments'];
                negcurves = [negcurves  negmoments'];
                legends = [legends ; [char(joints(joint2)) ' ' num2str(angle2)] ];
            end
            
            % plot total moments on left side of figure
            subplot(3,3,3*joint1-2);
            plot(range1, poscurves);hold on;
            plot(range1, negcurves);
            labels;
            title([char(joints(joint1)) ': total moment']);
            
            % plot passive moments in middle column of figure
            subplot(3,3,3*joint1-1);
            plot(range1, pascurves);hold on;
            labels;
            title([char(joints(joint1)) ': passive moment']);
            
            % subtract passive moments and plot in rightmost column of figure
            subplot(3,3,3*joint1);
            plot(range1, poscurves-pascurves);hold on;
            plot(range1, negcurves-pascurves);
            labels;
            title([char(joints(joint1)) ': active = total -- passive']);
            
            legend(legends);
            

            %======================================================================
            %> @brief Adds labels to plot.
            %======================================================================
            function labels
                a = get(gca);
                axis([min(range1) max(range1) a.YLim]);
                plot([0 0],a.YLim,'k:');
                plot(a.XLim,[0 0],'k:');
                hold off;
                xlabel('angle (deg)');
                ylabel('moment (Nm)');
            end

        end
        
        
        %======================================================================
        %> @brief Function to plot isokinetic moment-angular velocity curves
        %>
        %> @param testCase
        %> @param joint
        %> @param range
        %======================================================================
        function isokinetic_curves(testCase, joint, range)
            
            joints = testCase.model.joints.Properties.RowNames;
            angles = zeros(testCase.model.nJoints,1)+1e-6;
            angles([2 5]) = -5;			% knees 5 deg flexed to avoid high passive moment
            angvel = zeros(testCase.model.nJoints,1);
            pasmoments = [];
            posmoments = [];
            negmoments = [];
            for vel = range
                angvel(joint) = vel;
                pasmoments = [pasmoments testCase.maxmoment(joint, angles, angvel, 0)];
                posmoments = [posmoments testCase.maxmoment(joint, angles, angvel, 1)];
                negmoments = [negmoments testCase.maxmoment(joint, angles, angvel, -1)];
            end
            
            % plot total moments on left side of figure
            subplot(3,3,3*joint-2);
            plot(range, posmoments);hold on;
            plot(range, negmoments);
            labels;
            title([char(joints(joint)) ': total moment']);
            
            % plot passive moments in middle column of figure
            subplot(3,3,3*joint-1);
            plot(range, pasmoments);hold on;
            labels;
            title([char(joints(joint)) ': passive moment']);
            
            % subtract passive moments and plot in rightmost column of figure
            subplot(3,3,3*joint);
            plot(range, posmoments-pasmoments);hold on;
            plot(range, negmoments-pasmoments);
            labels;
            title([char(joints(joint)) ': active = total -- passive']);
            
            %======================================================================
            %> @brief Adds labels to plot.
            %======================================================================
            function labels
                a = get(gca);
                axis([min(range) max(range) a.YLim]);
                plot([0 0],a.YLim,'r:');
                plot(a.XLim,[0 0],'r:');
                hold off;
                xlabel('angular velocity (deg/s)');
                ylabel('moment (Nm)');
            end
            
        end
        
        
        %======================================================================
        %> @brief Function to simulate maximum moment for test_apparel()
        %>
        %> @details
        %> Simulates maximum moment at one joint, as function of all joint angles and angular velocities.
        %>
        %> @param testCase
        %> @param joint		For which joint we will calculate moment.
        %> @param angles	The six joint angles (deg).
        %> @param angvel	The six angular velocities (deg/s).
        %> @param sign		0: passive, 1: max positive moment, -1: max negative moment.
        %> @retval mom
        %======================================================================
        function [mom] = maxmoment(testCase, joint, angles, angvel, sign)
            
            angles = angles*pi/180;		% convert to radians
            angvel = angvel*pi/180;
            
            % determine moment arms so we know which muscles to activate
            % here (2D model) we have constant moment arms.
            % we should in a later version ask the MEX function what the moment arms are at this posture
            momentarms_oneside = [	0.05	0		0		; ...
                -0.062	0		0		; ...
                -0.072	-0.034	0		; ...
                0.034	0.05	0		; ...
                0		0.042	0		; ...
                0		-0.02	-0.053	; ...
                0		0		-0.053	; ...
                0		0		0.037	];
            % make a matrix that contains moment arms for both sides (in case needed...)
            momentarms = [	momentarms_oneside 		zeros(testCase.model.nMus/2,testCase.model.nJoints/2) ; ...
                zeros(testCase.model.nMus/2,testCase.model.nJoints/2) 	momentarms_oneside];
            
            Act =  sign*momentarms(:,joint) > 0;	% vector that has a 1 for all muscles we want to activate
            
            % determine lengthening velocities of the muscle-tendon complexes, normalize to Lceopt
            Vmuscle = -(momentarms * angvel) ./ testCase.model.muscles.lceopt;
            
            % determine the Lce's at which there is contraction equilibrium (dF/dt = 0, or Lcedot = Vmuscle)
            x = [zeros(3,1); angles; zeros(3,1); angvel; zeros(testCase.model.nMus,1); Act; zeros(4*testCase.model.nCPs,1)];
            xdot = [zeros(2*testCase.model.nDofs,1) ; Vmuscle ; zeros(testCase.model.nMus,1) ; zeros(4*testCase.model.nCPs,1)];	% we want these state derivatives
            u = zeros(testCase.model.nControls,1);				% no stim, we don't care about activation dynamics here
            
            % solve the equilibrium equation with Matlab fzero function, one muscle at a time
            for imus=1:testCase.model.nMus
                Lce = 1.0;		% initial guess for this muscle's Lce
                [Lce, Fval, Flag] = fzero(@contraction_equilibrium, Lce);
            end
            
            if (flag < 0)
                fprintf('maxmoment: muscle contraction equilibrium not found within max number of iterations.\n');
                keyboard
            end
            
            % now determine the joint moments at this state of the system
            moments = testCase.model.getJointmoments(x,u);
            moments = moments(~strcmp(testCase.model.dofs.joint, 'ground_pelvis')); % compare Gait2dc_acc.getJointmoments() => Dofs of joints
            mom = moments(joint);
              
            
            %======================================================================
            %> @brief Function to get contraction equilibrium for test_apparel()
            %>
            %> @param Lce
            %> @retval F
            %======================================================================
            function [F] = contraction_equilibrium(Lce)
                x(2*testCase.model.nDofs+imus) = Lce;
                f =  testCase.model.getDynamics(x,xdot,u);
                F = f(2*testCase.model.nDofs+imus);
            end
            
            
        end
        
        
        %======================================================================
        %> @brief Function to solve implicit differential equation for simulaton for contact points
        %>
        %> @details
        %> Solves the implicit differential equation f(x,dx/dt) = 0 for the contact variables.
        %>
        %> @param testCase
        %> @param xskeleton
        %> @param xinit
        %> @param times
        %> @param options
        %>
        %> @retval tout
        %> @retval xout
        %> @retval info
        %> @retval neval
        %======================================================================
        function [tout, xout, info, neval] = stepgrf(testCase, xskeleton,xinit,times,options)
            
            nsteps = numel(times);
            
            % get indices in the state vector
            idxAllCPs     = find(strcmp(testCase.model.states.type, 'xc') | ...
                strcmp(testCase.model.states.type, 'yc') | ...
                strcmp(testCase.model.states.type, 'Fx') | ...
                strcmp(testCase.model.states.type, 'Fy'));			% indices of the contact equations in f
            idxq          = find(strcmp(testCase.model.states.type, 'q'));
            idxqdot       = find(strcmp(testCase.model.states.type, 'qdot'));
            
            x = xinit;
            xdot = zeros(size(x));
            xsk = xskeleton(:,1);
            
            xout = x;
            tout = times(1);
            
            info = 0;
            neval = 0;
            for iStep = 2 : nsteps
                dt = times(iStep)-times(iStep-1);
                u = zeros(testCase.model.nControls,1);
                xnew = x + xdot*dt + rand(size(xdot));
                xsknew = xskeleton(:,iStep);
                xsknew(idxAllCPs) = xnew;
                for jIter = 1:options.maxiterations
                    xsknew(idxqdot) = (xsknew(idxq)-xsk(idxq))/dt;		% make sure xsknew has the skeleton velocities
                    [f,fx,fxdot] = testCase.model.getDynamics(xsknew, (xsknew-xsk)/dt, u);
                    neval = neval + 1;
                    Jc = fx(idxAllCPs,idxAllCPs)' + fxdot(idxAllCPs,idxAllCPs)'/dt;			% transpose because of how gait2dc was written
                    fc = f(idxAllCPs);
                    neval = neval + 1;
                    normfc = norm(fc);
                    dx = -Jc\fc;
                    xnew = xnew + dx;
                    xsknew(idxAllCPs) = xnew;
                    if (max(abs(fc)) < options.tol )
                        break
                    end
                    % backtrack in factors of 2 if no improvement from previous iteration
                    % 14 times will nearly go down to machine precision, we never need more
                    for k=1:14
                        % BE
                        f = gait2dc('Dynamics',xsknew, (xsknew-xsk)/dt, u);
                        fc = f(idxAllCPs);
                        neval = neval + 1;
                        if (norm(fc) < normfc)
                            break
                        end
                        dx = dx/2;
                        xnew = xnew - dx;
                        xsknew(idxAllCPs) = xnew;
                    end
                end
                if (jIter==options.maxiterations)
                    info = -1;
                    return
                end
                xdot = (xnew-x)/dt;
                x = xnew;
                xsk = xsknew;
                xout = [xout  x];
                tout = [tout; times(iStep)];
            end
        end    
    end
       
    
    
end