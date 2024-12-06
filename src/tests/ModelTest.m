%======================================================================
%> @file ModelTest.m
%> @brief Maltab m-file containing abstract test class for class Model
%>
%> @author Marlies Nitschke
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Abstract test class to test the class Model
%>
%> Tests:
%> - derivativetest_Dynamics: tests the derivatives of getDynamics()
%> - derivativetest_GRF: tests the derivatives of getGRF()
%> - derivativetest_Moments: tests the derivatives of getJointmoments()
%> - derivativetest_SimuAccGyro_Acc: tests the acc derivatives of simuAccGyro()
%> - derivativetest_SimuAccGyro_Gyro tests the gyro derivatives of simuAccGyro()
%> - test_speedOfMex: tests the execution time of the multibody dynamics
%> - test_memory: tests the memory usage of the multibody dynamics
%> - test_dynamics: tests the correctness of the dynamic model
%> - test_grf: tests the ground reaction force model
%> - test_showStickNeutral: tests the function showStick() for neutral position
%> - test_simulateFreefall: simulates freefall to test muscle dynamics
%>
%>
%> @todo Consider using Parameterized Tests later (introduced 2014a). http://de.mathworks.com/help/matlab/matlab_prog/create-basic-parameterized-test.html
%> @todo Consider to tag the tests later (introduced 2015a). https://de.mathworks.com/help/matlab/matlab_prog/tag-unit-tests.html#
%>       E.G.: (Test, TestTags = {'Acc_2D', 'objfun'}) => use a tag for the conditions (e.g Acc_2D) and one for the test (e.g. objfun)
%>       Then you can test single tags separalty. For example only tests related to tracking of accelerometer or only tests for objfun().
%======================================================================
classdef (Abstract)ModelTest < matlab.unittest.TestCase
    
    properties
        %> Model object
        model % There could be a smarter solution using parameterized tests!!!
    end
    properties(Access = protected) 
        %> Data for tracking tests
        data
        %> Used in test_grf() to keep the previous result as initial guess
        grf_xc
    end
    properties (Constant,Access = protected)
        %> Tolerance to pass derivative tests: maxerror/value < tol
        TOL_DERIVATIVETEST =  10^(-3);
        %> Tolerance to pass tests for Model.simuAccGyro()
        TOL_S = 10^(-3); 
        %> Names of different test conditions
        COMMANDS =  {'Acc', 'Gyro'};
    end
    
    %> @name Tests: ModelTest
    %> @{
    methods (Test, TestTags = {'Derivative','getDynamics'})
        
        %-----------------------------------------------------------------
        %> @brief Function to test the derivatives of Model.getDynamics()
        %>
        %> @details
        %> Gait2dc
        %> 
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Dynamics
        %> dfdx: Max. difference is:   0.0001211260 (3450.3251510427 vs. 3450.3250299167)  at  'dofs: hip_flexion_l: equations of motion from Autolev'(idx:16) and 's: Glutei_l'(idx:28) 
        %> dfdxdot: Max. difference is:   0.0000023285 ( -0.0212161829 vs.  -0.0212185114)  at  'dofs: hip_flexion_l: equations of motion from Autolev'(idx:16) and 'qdot: pelvis_ty'(idx:11) 
        %> dfdu: Max. difference is:   0.0000068123 (-152.4926607122 vs. -152.4926675245)  at  'muscles: Rectus_l: activation dynamics: da/dt - (u-a)(c1*u + c2) = 0'(idx:46) and 'Rectus_l'(idx:12) 
        %> @endcode 
        %>
        %> @image html Gait2dcTest_derivativetest_Dynamics.png width=800px
        %> 
        %>
        %> Gait3d
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Dynamics
        %> dfdx: Max. difference is:   0.0027930663 (-46136.0346080336 vs. -46136.0318149673)  at  'muscles: add_brev_l: contraction dynamics: f = Fsee - Fce - Fpee'(idx:122) and 's: add_brev_l'(idx:122) 
        %> dfdxdot: Max. difference is:   0.0010051358 (-494.0920577725 vs. -494.0930629083)  at  'muscles: per_long_r: contraction dynamics: f = Fsee - Fce - Fpee'(idx:106) and 's: per_long_r'(idx:106) 
        %> dfdu: Max. difference is:   0.0000086362 (350.9248257347 vs. 350.9248170985)  at  'muscles: glut_max3_l: activation dynamics: da/dt - rate * (u-a) = 0'(idx:223) and 'glut_max3_l'(idx:65) 
        %> @endcode 
        %>
        %> @image html Gait3dTest_derivativetest_Dynamics.png width=800px
        %-----------------------------------------------------------------
        function derivativetest_Dynamics(testCase)
            iCom = 0;
            % setup the Model object and data
            testCase = testCase.setup_Test_Random(iCom);
            
            % Do the test: Evaluate the Jacobians at random inputs and show sparsity
            % pattern
            x = testCase.generateRandType('states');
            xd = testCase.generateRandType('states'); %> @todo How do we set this correctly? We do not have bounds for the derivative!
            u = testCase.generateRandType('controls');
            [f, dfdx, dfdxdot, dfdu] = testCase.model.getDynamics(x,xd,u);
            dfdx = dfdx';
            dfdxdot = dfdxdot';
            dfdu = dfdu';
            figure('Name', 'derivativetest_Dynamics');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            subplot(1, 3, 1); spy(dfdx);title('dfdx');xlim([1 size(dfdx, 2)]); ylim([1 size(dfdx, 1)]);
            subplot(1, 3, 2); spy(dfdxdot);title('dfdxdot');xlim([1 size(dfdxdot, 2)]); ylim([1 size(dfdxdot, 1)]);
            subplot(1, 3, 3); spy(dfdu);title('dfdu');xlim([1 size(dfdu, 2)]); ylim([1 size(dfdu, 1)]);
            
            % Make sure that there are no hard zeros in the "nonzero
            % elements" of the Jacobians
            if nnz(dfdx) ~= nnz(dfdx*1)
                disp('WARNING: dfdx contains hard zeros');
                disp('Hit ENTER to continue');
                pause;
            end
            if nnz(dfdxdot) ~= nnz(dfdxdot*1)
                disp('WARNING: dfdxdot contains hard zeros');
                disp('Hit ENTER to continue');
                pause;
            end
            if nnz(dfdu) ~= nnz(dfdu*1)
                disp('WARNING: dfdu contains hard zeros');
                disp('Hit ENTER to continue');
                pause;
            end
            
            %Check the Jacobians by comparing to finite difference
            %calculations
            hh = 1e-7;
            dfdx_num = zeros(size(dfdx));
            dfdxdot_num = zeros(size(dfdxdot));
            dfdu_num = zeros(size(dfdu));
            for i = 1:testCase.model.nStates
                tmp = x(i);
                x(i) = x(i) + hh;
                fnew = testCase.model.getDynamics(x,xd,u);
                dfdx_num(:, i) = (fnew-f)/hh;
                x(i) = tmp;
                
                tmp = xd(i);
                xd(i) = xd(i) + hh;
                fnew = testCase.model.getDynamics(x,xd,u);
                dfdxdot_num(:, i) = (fnew -f)/hh;
                xd(i) = tmp;
            end
            
            for i = 1:testCase.model.nControls
                tmp = u(i);
                u(i) = u(i) +hh;
                fnew = testCase.model.getDynamics(x,xd,u);
                dfdu_num(:, i) = (fnew-f)/hh;
                u(i) = tmp;
            end
            
            % compare all the Jacobians to their numerical estimates
            fprintf('\n');
            disp('Derivatives for Dynamics');
            delimiter(1:testCase.model.nConstraints) = {': '};
            fNames = strcat(testCase.model.constraints.type, delimiter', ...
                            testCase.model.constraints.name, delimiter', ...
                            testCase.model.constraints.equation); clear delimiter;
            delimiter(1:testCase.model.nStates) = {': '};
            xNames = strcat(testCase.model.states.type, delimiter', testCase.model.states.name);
            uNames = testCase.model.controls.name; 
            testCase.matcompare(dfdx   , dfdx_num   , 'dfdx'   ,fNames,xNames);
            testCase.matcompare(dfdxdot, dfdxdot_num, 'dfdxdot',fNames,xNames);
            testCase.matcompare(dfdu   , dfdu_num   , 'dfdu'   ,fNames,uNames);
        end
        
    end
    
    methods (Test, TestTags = {'Derivative','getGRF'})
        %-----------------------------------------------------------------
        %> @brief Function to test the derivatives of Model.getGRF()
        %>
        %> @details
        %> Gait2dc
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for GRF
        %> dGRFdx: Max. difference is:   0.0000000055 ( -0.0688392154 vs.  -0.0688392099)  at  'right Mz'(idx:6) and 'Fy: toe_r'(idx:58) 
        %> @endcode 
        %>
        %> @image html Gait2dcTest_derivativetest_GRF.png width=800px
        %>
        %> Gait3d
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for GRF
        %> dGRFdx: Max. difference is:   0.0000000187 ( -0.9337272476 vs.  -0.9337272289)  at  'left My'(idx:11) and 'zc: CPL_l'(idx:337) 
        %> @endcode 
        %>
        %> @image html Gait3dTest_derivativetest_GRF.png width=800px
        %-----------------------------------------------------------------
        function derivativetest_GRF(testCase)
            iCom = 0;
            % setup the Model object and data
            testCase = testCase.setup_Test_Random(iCom);
            
            % Do the test: Test the correctness of the GRF Jacobian at a random state
            x = testCase.generateRandType('states');
            [GRF, dGRFdx] = testCase.model.getGRF(x);
            dGRFdx = dGRFdx';
            figure('Name', 'derivativetest_GRF');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            spy(dGRFdx); title('dGRFdx');xlim([1 size(dGRFdx, 2)]); ylim([1 size(dGRFdx, 1)]);
            
            hh = 1e-7;
            dGRFdx_num = zeros(size(dGRFdx));
            for i = 1:testCase.model.nStates
                tmp = x(i);
                x(i) = x(i) + hh;
                GRFnew = testCase.model.getGRF(x);
                dGRFdx_num(:, i) = (GRFnew - GRF)/hh;
                x(i) = tmp;
            end
            
            % make sure there are no hard zeros (this crashes IPOPT
            % beacause the sparsity structure is not consistent with the
            % sparse matrix contents)
            if nnz(dGRFdx) ~= nnz(dGRFdx*1)
                disp('WARNING: dGRFdx contains hard zeros');
                disp('Hit ENTER to continue');
                pause;
            end
            
            % compare all the Jacobian to their numerical estimates
            fprintf('\n');
            disp('Derivatives for GRF');
            grfNames = testCase.model.GRFNAMES;
            delimiter(1:testCase.model.nStates) = {': '};
            xNames = strcat(testCase.model.states.type, delimiter', testCase.model.states.name);
            testCase.matcompare(dGRFdx   , dGRFdx_num   , 'dGRFdx', grfNames, xNames);
        end
        
    end
    
    methods (Test, TestTags = {'Derivative','getJointmoments'})
        %-----------------------------------------------------------------
        %> @brief Function to test the derivatives of Model.getJointmoments()
        %>
        %> @details
        %> Gait2dc
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Moments
        %> dMdx: Max. difference is:   0.0120983138 (41480.6475985067 vs. 41480.6355001929)  at  'hip_flexion_l'(idx:7) and 's: Glutei_l'(idx:28) 
        %> dMdu: Max. difference is:   0.0000000000 (  0.0000000000 vs.   0.0000000000)  at  'pelvis_tx'(idx:1) and 'Iliopsoas_r'(idx:1) 
        %> @endcode 
        %>
        %> @image html Gait2dcTest_derivativetest_Moments.png width=800px
        %>
        %> Gait3d
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Moments
        %> dMdx: Max. difference is:   0.1738820117 (-350944.6392260600 vs. -350944.8131080717)  at  'hip_flexion_r'(idx:7) and 'q: hip_flexion_r'(idx:7) 
        %> dMdu: Max. difference is:   0.0000000000 (  0.0000000000 vs.   0.0000000000)  at  'pelvis_tx'(idx:1) and 'Iliopsoas_r'(idx:1) 
        %> @endcode 
        %>
        %> @image html Gait3dTest_derivativetest_Moments.png width=800px
        %-----------------------------------------------------------------
        function derivativetest_Moments(testCase)
            iCom = 0;
            % setup the Model object and data
            testCase = testCase.setup_Test_Random(iCom);
            
            % Do the test: Test the correctness of the joint moment Javobian at a random state
            x = testCase.generateRandType('states');
            u = testCase.generateRandType('controls');
            [Mout, dMdx, dMdu] = testCase.model.getJointmoments(x, u);
            dMdx = dMdx';
            dMdu = dMdu';
            figure('Name', 'derivativetest_Moments');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            subplot(1, 2, 1); spy(dMdx); title('dMdx');xlim([1 size(dMdx, 2)]); ylim([1 size(dMdx, 1)]);
            subplot(1, 2, 2); spy(dMdu); title('dMdu');xlim([1 size(dMdu, 2)]); ylim([1 size(dMdu, 1)]);
            
            hh = 1e-7;
            dMdx_num = zeros(size(dMdx));
            dMdu_num = zeros(size(dMdu));
            for i = 1:testCase.model.nStates
                tmp = x(i);
                x(i) = x(i) + hh;
                Mnew = testCase.model.getJointmoments(x,u);
                dMdx_num(:, i) = (Mnew -Mout)/hh;
                x(i) = tmp;
            end
            
            for i = 1:testCase.model.nControls
                tmp = u(i);
                u(i) = u(i) +hh;
                Mnew = testCase.model.getJointmoments(x,u);
                dMdu_num(:, i) = (Mnew -Mout)/hh;
                u(i) = tmp;
            end
            
            % make sure there are no hard zeros (this crashes IPOPT
            % because the sparsity structure is not consistent with the
            % sparse matrix contents)
            if nnz(dMdx) ~= nnz(dMdx*1)
                disp('WARNING: dMdx contains hard zeros');
                disp('Hit ENTER to continue');
                pause;
            end
            if nnz(dMdu) ~= nnz(dMdu*1)
                disp('WARNING: dMdu contains hard zeros');
                disp('Hit ENTER to continue');
                pause;
            end
            
            % compare all the Jacobian to their numerical estimates
            fprintf('\n');
            disp('Derivatives for Moments');
            MNames = testCase.model.dofs.Properties.RowNames;
            delimiter(1:testCase.model.nStates) = {': '};
            xNames = strcat(testCase.model.states.type, delimiter', testCase.model.states.name);
            uNames = testCase.model.controls.name;
            testCase.matcompare(dMdx   , dMdx_num   , 'dMdx', MNames, xNames);
            testCase.matcompare(dMdu   , dMdu_num   , 'dMdu', MNames, uNames);
        end
        
    end
    
    methods (Test, TestTags = {'Derivative','simuAccGyro'})
        %-----------------------------------------------------------------
        %> @brief Function to test the acc derivatives of Model.simuAccGyro()
        %>
        %> @details
        %> Gait2dc
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Acc
        %> dsdq: Max. difference is:   0.0000000227 (  1.7175791416 vs.   1.7175791189)  at  '4'(idx:4) and 'hip_flexion_r'(idx:4) 
        %> dsdqd: Max. difference is:   0.0000000229 ( -0.8737816552 vs.  -0.8737816781)  at  '4'(idx:4) and 'pelvis_tilt'(idx:3) 
        %> dsdqdd: Max. difference is:   0.0000000271 (  0.2436197648 vs.   0.2436197377)  at  '4'(idx:4) and 'pelvis_tx'(idx:1) 
        %> @endcode 
        %>
        %> @image html Gait2dcTest_derivativetest_SimuAccGyro_Acc.png width=800px
        %>
        %> Gait3d
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Acc
        %> dsdq: Max. difference is:   0.0000016256 (  1.7379826469 vs.   1.7379810213)  at  '10'(idx:7) and 'pelvis_tilt'(idx:1) 
        %> dsdqd: Max. difference is:   0.0000000929 ( -1.5895990794 vs.  -1.5895991723)  at  '10'(idx:7) and 'hip_rotation_l'(idx:16) 
        %> dsdqdd: Max. difference is:   0.0000000209 (  0.0951314259 vs.   0.0951314050)  at  '5'(idx:4) and 'pelvis_list'(idx:2)  
        %> @endcode 
        %>
        %> @image html Gait3dTest_derivativetest_SimuAccGyro_Acc.png width=800px
        %-----------------------------------------------------------------
        function derivativetest_SimuAccGyro_Acc(testCase)
            iCom = 1;
            % setup the Model object and data
            testCase = testCase.setup_Test_Random(iCom);
            % do test
            testCase.doTest_derivativetest_SimuAccGyro(iCom);
        end
        
        %-----------------------------------------------------------------
        %> @brief Function to test the gyro derivatives of odel.simuAccGyro()
        %>
        %> @details
        %> Gait2dc
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Gyro
        %> dsdq: Max. difference is:   0.0000000000 (  0.0000000000 vs.   0.0000000000)  at  '1'(idx:1) and 'pelvis_tx'(idx:1) 
        %> dsdqd: Max. difference is:   0.0000000061 (  1.0000000000 vs.   0.9999999939)  at  '4'(idx:4) and 'pelvis_tilt'(idx:3) 
        %> dsdqdd: Max. difference is:   0.0000000000 (  0.0000000000 vs.   0.0000000000)  at  '1'(idx:1) and 'pelvis_tx'(idx:1) 
        %> @endcode 
        %>
        %> @image html Gait2dcTest_derivativetest_SimuAccGyro_Gyro.png width=800px
        %>
        %> Gait3d
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Derivatives for Gyro
        %> dsdq: Max. difference is:   0.0000000843 (  0.6910369209 vs.   0.6910368366)  at  '3'(idx:2) and 'hip_adduction_r'(idx:8) 
        %> dsdqd: Max. difference is:   0.0000000081 (  0.8115557434 vs.   0.8115557515)  at  '5'(idx:4) and 'hip_rotation_r'(idx:9) 
        %> dsdqdd: Max. difference is:   0.0000000000 (  0.0000000000 vs.   0.0000000000)  at  '2'(idx:1) and 'pelvis_tilt'(idx:1) 
        %> @endcode 
        %>
        %> @image html Gait3dTest_derivativetest_SimuAccGyro_Gyro.png width=800px
        %-----------------------------------------------------------------
        function derivativetest_SimuAccGyro_Gyro(testCase)
            iCom = 2;
            % setup the Model object and data
            testCase = testCase.setup_Test_Random(iCom);
            % do test
            testCase.doTest_derivativetest_SimuAccGyro(iCom);
        end
  
    end
    
    methods (Test, TestTags = {'Speed','getDynamics'})
        %=================================================================
        %> @brief Function for testing the speed of the MEX function.
        %>
        %> @details
        %> This tests the execution time of the multibody dynamics,
        %> without and with derivatives.  The MEX function is executed
        %> 1000 times with random inputs, and the resulting time is
        %> divided by 1000.
        %> 
        %> Gait2dc
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Speed test
        %> Evaluation of dynamics function without derivatives:    1.100 ms
        %> Evaluation of dynamics function with derivatives:       0.831 ms
        %> @endcode 
        %>
        %> Gait3d
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Speed test
        %> Evaluation of dynamics function without derivatives:    3.340 ms
        %> Evaluation of dynamics function with derivatives:       4.998 ms
        %> @endcode 
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_speedOfMex(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % create random data
            N = 1000;
            x = testCase.generateRandType('states');
            xd = testCase.generateRandType('states'); %> @todo How do we set this correctly? We do not have bounds for the derivative!
            u = testCase.generateRandType('controls');
            
            % run tests
            t1 = zeros(N,1);
            for i=1:N
                tic;
                [f] = testCase.model.getDynamics(x,xd,u);
                t1(i) = toc;
            end
            evaltime1 = 1000*mean(t1);
            
            t2 = zeros(N,1);
            for i=1:N
                tic;
                [f, dfdx, dfdxdot, dfdu] = testCase.model.getDynamics(x,xd,u);
                t2(i) = toc;
            end
            evaltime2 = 1000*mean(t2);
            
            % Prints 95th percentiles of the worst evaluation times
            perc1 = sort(t1);
            perc1 = perc1(int64(0.95*N))*1000;
            perc2 = sort(t2);
            perc2 = perc2(int64(0.95*N))*1000;
            
            fprintf('\n');
            fprintf('Speed test\t\t\t\t\t\tavg\t\t95th perc. \n');
            fprintf('Evaluation of dynamics function without derivatives: %8.3f ms     %8.3f ms\n', evaltime1,perc1);
            fprintf('Evaluation of dynamics function with derivatives:    %8.3f ms     %8.3f ms\n', evaltime2,perc2);
        end
        
    end
     
    methods (Test, TestTags = {'Speed','simuAccGyro'})
        %=================================================================
        %> @brief Function for testing the speed of Model.simuAccGyro().
        %>
        %> @details
        %> This tests the execution time of Model.simuAccGyro(), 
        %> without and with derivatives.  The simuAccGyro() is executed 
        %> 1000 times with random inputs, and the resulting time is 
        %> divided by 1000.  
        %>
        %> Gait2dc
        %>
        %> Example output for 3 (acc x, acc y, acc z) times 7 (trunk, legs, feet) IMU signals:
        %> @code{.unparsed}
        %> Speed test of SimuAccGyro()
        %> Evaluation of dynamics function without Jacobians:    1.714 ms
        %> Evaluation of dynamics function with Jacobians:       3.047 ms
        %> @endcode 
        %>
        %> Gait3d
        %>
        %> Example output for 6 (all acc and gyro) times 7 (trunk, legs, feet) IMU signals::
        %> @code{.unparsed}
        %> Speed test of simuAccGyro()
        %> Evaluation of dynamics function without Jacobians:    4.066 ms
        %> Evaluation of dynamics function with Jacobians:      26.996 ms
        %> @endcode
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_speedOfSimuAccGyro(testCase)
            % setup the Model object 
            iCom = 3; % all sensors which are used in reality 
            testCase = testCase.setup_Test_Random(iCom);
            
            % create random data
            N = 1000;
            q = testCase.generateRandType('q', 1);
            qd = testCase.generateRandType('qdot', 1);
            qdd = testCase.generateRandType('qdot', 1); %> @todo How do we set this correctly? We do not have bounds for the derivative qdd!
            
            % run tests
            tic;
            for i=1:N
                [s] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            end
            evaltime1 = 1000*toc/N;

            tic;
            for i=1:N
                [s, ds_dq, ds_dqd, ds_dqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            end
            evaltime2 = 1000*toc/N;
            
            fprintf('\n');
            fprintf('Speed test of simuAccGyro()\n');
            fprintf('Evaluation of dynamics function without Jacobians: %8.3f ms\n', evaltime1);
            fprintf('Evaluation of dynamics function with Jacobians:    %8.3f ms\n', evaltime2);
        end
        
    end
        
    methods (Test, TestTags = {'Memory','getDynamics'})
        %=================================================================
        %> @brief Function for testing memory use of the MEX function
        %>
        %> @details
        %> This tests the memory usage of the multibody dynamics.
        %> The MEX function is executed 500000 times with random inputs,
        %> and the resulting memory usage is plotted.
        %>
        %> This function is only working at Windows!
        %>
        %> Gait2dc
        %> @image html Gait2dcTest_test_memory.png width=800px
        %>
        %> Gait3d
        %> @image html Gait3dTest_test_memory.png width=800px
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_memory(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            fprintf('\n');
            fprintf('Memory test\n');
            
            % we do this only for the dynamics command
            N = 500000;
            userview = memory;
            mem_baseline = userview.MemUsedMATLAB;
            mem = [];
            tic;
            for i = 1:N
                x = testCase.generateRandType('states');
                xd = testCase.generateRandType('states'); %> @todo How do we set this correctly? We do not have bounds for the derivative!
                u = testCase.generateRandType('controls');
                [f, dfdx, dfdxdot, dfdu] = testCase.model.getDynamics(x,xd,u);
                if toc > 5.0
                    userview = memory;
                    mem = [mem ; i userview.MemUsedMATLAB];
                    fprintf('Done %d function calls out of %d\n', i, N);
                    tic;
                end
            end
            figure('Name', 'test_memory')
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            plot(mem(:,1),(mem(:,2)-mem_baseline)/1024/1024,'-o');
            xlabel('MEX function evaluations');
            ylabel('Memory used (MB)');
        end
        
    end
    
    methods (Test, TestTags = {'getDynamics', 'getJointmoments', 'getGRF'})
        %=================================================================
        %> @brief Function for testing the correctness of the dynamic model.
        %>
        %> @details
        %> Free fall dynamics test:
        %> This test shows the correctness of the dynamic model
        %> f = f(x,xdot,u,M) = 0.
        %> Therefore, the model is put in a freefall acceleration state,
        %> so residuals should be zero when acceleration is -g for DOF 5
        %> and zero elsewhere. Moreover, the joint moments of the model
        %> in its neutral position are computed and should be zero.
        %>
        %>
        %> Gait2dc
        %>
        %> Dynamics violations:
        %> @image html Gait2dcTest_test_dynamics_dynamics.png width=800px
        %> Joint moments:
        %> @image html Gait2dcTest_test_dynamics_moments.png width=800px
        %> GRFs:
        %> @image html Gait2dcTest_test_dynamics_GRFs.png width=800px
        %>
        %> Gait3d
        %>
        %> Dynamics violations:
        %> @image html Gait3dTest_test_dynamics_dynamics.png width=800px
        %> Joint moments:
        %> @image html Gait3dTest_test_dynamics_moments.png width=800px
        %> GRFs:
        %> @image html Gait3dTest_test_dynamics_GRFs.png width=800px
        %=================================================================
        function test_dynamics(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % put the model in a freefall acceleration state, so residuals should be zero when acceleration is -g for DOF 5 and zero elsewhere.
            x = testCase.model.states.xneutral;
            xd = zeros(size(x));
            xd(testCase.model.nDofs + 2) = -9.81;
            u = zeros(testCase.model.nControls,1);
            f = testCase.model.getDynamics(x, xd, u);
            
            % compute how far f is outside of the dynamics constraints
            f = min(0, f - testCase.model.constraints.fmin) + max(0, f - testCase.model.constraints.fmax);
            
            fprintf('\n');
            disp('Test dynamics');
            
            figure('Name', 'test_dynamics:dynamics violation')
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            plot(f);
            xlim([1 length(f)]);
            xlabel('equation number');
            ylabel('dynamics violation');
            title('dynamics violation');
            disp('In Figure ''dynamics violation'', check that all dynamics violations are close enough to zero.');
            answer = input('Y for ok, N for failed\n','s');
            testCase.verifyMatches(lower(answer), 'y', 'You entered that dynamics violations are close enough to zero');
            
            % also check the Jointmoments function with the model in its neutral position
            M = testCase.model.getJointmoments(x, u);
            figure('Name', 'test_dynamics:joint moment')
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            plot(M);
            xlim([1 length(M)]);
            xlabel('DOF number');
            ylabel('joint moment');
            title('joint moment');
            disp('In Figure ''joint moment'', check that the joint moments are close enough to zero.');
            answer = input('Y for ok, N for failed\n','s');
            testCase.verifyMatches(lower(answer), 'y', 'You entered that, joint Moments are close enough to zero');
            
            % also check the GRFs function with the model in its neutral position
            GRF = testCase.model.getGRF(x);
            figure('Name', 'test_dynamics:GRF')
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            plot(GRF);
            xlim([1 length(GRF)]);
            ylabel('GRF');
            title('GRF');
            disp('In Figure ''GRF'', check that the GRFs are close enough to zero.');
            answer = input('Y for ok, N for failed\n','s');
            testCase.verifyMatches(lower(answer), 'y', 'You entered that, joint Moments are not close enough to zero');
        end
        
    end
    
    methods (Test, TestTags = {'getGRF'})
        %=================================================================
        % Declared function for testing the ground reaction force model
        %=================================================================
        test_grf(testCase)
        
    end
    
    methods (Test, TestTags = {'Neutral', 'showStick'})
        %=================================================================
        %> @brief Function for testing showStick()
        %>
        %> @details
        %> This tests plots a stick figure for the neutral position of the
        %> model.
        %>
        %> Gait2dc
        %> @image html Gait2dcTest_test_showStickNeutral.png width=800px
        %>
        %> Gait3d
        %> @image html Gait3dTest_test_showStickNeutral.png width=800px
        %>
        %> @todo 
        %> Add error checking
        %=================================================================
        function test_showStickNeutral(testCase)
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % get neutral position
            xneutral = testCase.model.states.xneutral;
            
            figure('Name','test_showStickNeutral');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            % plot
            testCase.model.showStick(xneutral);
            
        end
    end
        
    methods (Test, TestTags = {'Freefall', 'getDynamics'})
        %=================================================================
        %> @brief Function to simulate freefall
        %>
        %> @details
        %> The model is placed in a neutral position. We drop the model on
        %> the ground. This is done twice, once with passive muscles and
        %> once with muscles that are 50% activated. It uses the implicit
        %> midpoint Euler method (van den Bogert et al., 2011) with a step
        %> size of 10 ms.
        %>
        %> The second simulation (with activated muscles) is then repeated
        %> with a step size of 1 ms. The results are compared with the 10 ms
        %> stepsize.
        %>
        %> Gait2dc
        %>
        %> Example outpu:
        %> @code{.unparsed}
        %> Running test_simulateFreefall:
        %> 1. Running passive freefall simulation...
        %> Number of time steps:           400
        %> Number of function evaluations: 2572
        %> 2. Running active freefall simulation...
        %> Number of time steps:           400
        %> Number of function evaluations: 2285
        %> 3. Running active freefall simulation with smaller steps...
        %> Number of time steps:           4000
        %> Number of function evaluations: 16259
        %> @endcode
        %> @image html Gait2dcTest_test_simulateFreefall.png width=800px
        %>
        %> Gait3d
        %>
        %> Example output:
        %> @code{.unparsed}
        %> Running test_simulateFreefall:
        %> 1. Running passive freefall simulation...
        %> Number of time steps:           400
        %> Number of function evaluations: 3120
        %> 2. Running active freefall simulation...
        %> Number of time steps:           1
        %> Number of function evaluations: 1129
        %> 3. Running active freefall simulation with smaller steps...
        %> Number of time steps:           4000
        %> Number of function evaluations: 17836
        %> @endcode
        %>
        %> @todo The second simulation is not solving for 3D. And switching from
        %> Windows to Linux caused that also the first simulation is not
        %> solving. A difference in the compilers maybe caused this...
        %=================================================================
        function test_simulateFreefall(testCase)
            fprintf('\n');
            disp('Running test_simulateFreefall:');
            
            % duration of simulation
            tend = 4.0;
            
            % set solver parameters for IMstep
            h = 0.01;
            nsteps = round(tend/h);
            options.tol = 1e-6;
            options.maxiterations = 100;
            
            % setup the Model object
            iCom = 0; % no data
            testCase = testCase.setup_Test_Random(iCom);
            
            % get neutral position
            xinit = testCase.model.states.xneutral;
            
            % run passive simulation
            disp('1. Running passive freefall simulation...');
            [tout1, xout1, info, neval] = testCase.step(@testCase.nostim, xinit,[0 tend],nsteps,options);
            fprintf('Number of time steps:           %d\n', numel(tout1)-1);
            fprintf('Number of function evaluations: %d\n', neval);
            figure('Name', 'test_simulateFreefall:passive simulation');clf;
            for i=1:numel(tout1)
                testCase.model.showStick(xout1(:,i));
                pause(0.02);
            end
            testCase.verifyGreaterThanOrEqual(info, 0, 'Running active freefall simulation: A time step did not converge. info was below 0')
            
            
            % run simulation with 50% muscle stim
            disp('2. Running active freefall simulation...');
            [tout2, xout2, info, neval] = testCase.step(@testCase.halfstim, xinit,[0 tend],nsteps,options);
            fprintf('Number of time steps:           %d\n', numel(tout2)-1);
            fprintf('Number of function evaluations: %d\n', neval);
            figure('Name', 'test_simulateFreefall:50% muscle stim');clf;
            for i=1:numel(tout2)
                testCase.model.showStick(xout2(:,i));
                pause(0.02);
            end
            testCase.verifyGreaterThanOrEqual(info, 0, 'Running active freefall simulation with smaller steps: A time step did not converge. info was below 0')
            
            
            % run same simulation with 1/10th step size
            disp('3. Running active freefall simulation with smaller steps...');
            [tout3, xout3, info, neval] = testCase.step(@testCase.halfstim, xinit,[0 tend],10*nsteps,options);
            fprintf('Number of time steps:           %d\n', numel(tout3)-1);
            fprintf('Number of function evaluations: %d\n', neval);
            testCase.verifyGreaterThanOrEqual(info, 0, 'Running active freefall simulation with smaller steps: A time step did not converge. info was below 0')
            
            
            % compare the solutions at two different step sizes
            figure('Name','test_simulateFreefall:Step size comparison');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            nRows = ceil(sqrt(testCase.model.nDofs));
            nCol = ceil(testCase.model.nDofs / nRows);
            for iDof = 1 : testCase.model.nDofs
                subplot(nRows, nCol,iDof)
                plot(tout2,xout2(iDof,:)',tout3,xout3(iDof,:)');
                set(gca,'XLim',[0 tend])
                title(testCase.model.dofs.Properties.RowNames{iDof}, 'Interpreter', 'None');
            end
            legend(['h = ' num2str(tend/nsteps)],['h = ' num2str(tend/nsteps/10)]);
            
            disp('In Figure ''test_simulateFreefall:Step size comparison'',check if the differences for the step sizes is small.');
            answer = input('Y for ok, N for failed\n','s');
            testCase.verifyMatches(lower(answer), 'y', 'You entered that, the differences for the step sizes is too large');
        end
        
    end
    %> @}
    
    methods (Access = protected)
        
        %-----------------------------------------------------------------
        % Declared function to setup the Model object and data
        %-----------------------------------------------------------------
        testCase = setup_Test_Random(testCase, iCom)
        
        %-----------------------------------------------------------------
        %> @brief Function with test content to test Model.simuAccGyro()
        %-----------------------------------------------------------------
        function testCase = doTest_derivativetest_SimuAccGyro(testCase, iCom)
            % test the Jacobians of the sensor models
            q = testCase.generateRandType('q', 1);
            qd = testCase.generateRandType('qdot', 1);
            qdd = testCase.generateRandType('qdot', 1); %> @todo How do we set this correctly? We do not have bounds for the derivative qdd!
            [s, dsdq, dsdqd, dsdqdd] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
            % show the sparsity patterns
            figure('Name',['derivativtest_SimuAccGyro' testCase.COMMANDS{iCom}])
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);pause(0.1);
            subplot(1, 3, 1);
            spy(dsdq, '.', 5);xlim([1 size(dsdq, 2)]); ylim([1 size(dsdq, 1)]);
            title('sparsity pattern of ds/dq');
            subplot(1, 3, 2);
            spy(dsdqd, '.', 5);xlim([1 size(dsdqd, 2)]); ylim([1 size(dsdqd, 1)]);
            title('sparsity pattern of ds/dqd');
            subplot(1, 3, 3);
            spy(dsdqdd, '.', 5);xlim([1 size(dsdqdd, 2)]); ylim([1 size(dsdqdd, 1)]);
            title('sparsity pattern of ds/dqdd');
            
            % numerical estimate of the sensor Jacobians
            hh = 1e-7;
            dsdq_num = zeros(size(dsdq));
            dsdqdot_num = zeros(size(dsdqd));
            dsdqdd_num = zeros(size(dsdqdd));
            for i = 1:testCase.model.nDofs
                tmp = q(i);
                q(i) = q(i) + hh;
                [snew] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
                dsdq_num(:, i) = (snew-s)/hh;
                q(i) = tmp;
                
                tmp = qd(i);
                qd(i) = qd(i) + hh;
                [snew] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
                dsdqdot_num(:, i) = (snew-s)/hh;
                qd(i) = tmp;
                
                tmp = qdd(i);
                qdd(i) = qdd(i) + hh;
                [snew] = testCase.model.simuAccGyro(testCase.data.variables, q, qd, qdd);
                dsdqdd_num(:, i) = (snew-s)/hh;
                qdd(i) = tmp;
            end
            
            % compare all the Jacobians to their numerical estimates
            fprintf('\n');
            disp(['Derivatives for ' testCase.COMMANDS{iCom}]);
            sNames = testCase.data.variables.name;
            sNames = cellfun(@num2str,sNames,'UniformOutput',false);
            qNames = testCase.model.states.name(strcmp(testCase.model.states.type, 'q'));        
            testCase.matcompare(dsdq   , dsdq_num     , 'dsdq'  , sNames, qNames);
            testCase.matcompare(dsdqd  , dsdqdot_num  , 'dsdqd' , sNames, qNames);
            testCase.matcompare(dsdqdd , dsdqdd_num   , 'dsdqdd', sNames, qNames);
        end
        
        %======================================================================
        % Declared function to solve grf by search for static equilibrium in contact variables.
        %======================================================================
        [grf] = solvegrf(testCase,y,vy,xx,vx,vxc)
        
        %======================================================================
        %> @brief Function to solve implicit differential equation for simulaton
        %>
        %> @details
        %> Solves the implicit differential equation fmin <= f(x,dx/dt,u) <= fmax using the 
        %> midpoint Euler method, with constant stepsize, and Newton iteration.
        %>
        %> @param testCase
        %> @param stimfun               A function u = fun(t) to stimulate the muscles.
        %> @param x0                    Initial system state.
        %> @param trange                [starttime enddtime]
        %> @param nsteps                How many time steps to take.
        %> @param options               Struct with fields:
        %>                              - tol: Tolerance
        %>                              - maxiterations: Maximum number of Newton iterations
        %>
        %> @retval tout                 Vector of time points (ntimes x 1).
        %> @retval xout                 Solution vectors (nstates x ntimes).
        %> @retval info                 0: success, -1: maximum number of iterations exceeded.
        %> @retval neval                Number of evaluations.
        %======================================================================
        function [tout, xout, info, neval] = step(testCase, stimfun, x0, trange, nsteps, options)
            
            t = trange(1);
            dt = diff(trange)/nsteps;
            
            x = x0;
            xdot = zeros(size(x));
            
            tout = t;
            xout = x;
            
            fmin = testCase.model.constraints.fmin;
            fmax = testCase.model.constraints.fmax;
            
            info = 0;
            neval = 0;
            for i=1:nsteps
                u = stimfun(t+dt/2);			% controls at midpoint of time interval
                xnew = x + xdot*dt;				% initial guess of state at t+dt
                for j = 1:options.maxiterations
                    % midpoint formula
                    % [f,fx,fxdot] = gait2dc('Dynamics',(xnew+x)/2, (xnew-x)/dt, u);
                    % J = fx'/2 + fxdot'/dt;		% transpose because of how gait2dc was written
                    % backward Euler formula
                    [f,fx,fxdot] = testCase.model.getDynamics(xnew, (xnew-x)/dt, u);
                    J = fx' + fxdot'/dt;			% transpose because of how gait2dc was written
                    neval = neval + 1;
                    sumViol = norm((fmin-f < 0) .*(fmin-f) + (fmax-f > 0) .*(fmax-f)); % norm of f if f is out of bounds
                    dx = -J\f;
                    xnew = xnew + dx;
                    % check if it is fmin-tol < f < fmax+tol
                    if fmin - options.tol < f &  f < fmax+options.tol
                        break
                    end
                    % backtrack in factors of 2 if no improvement from previous iteration
                    % 14 times will nearly go down to machine precision, we never need more
                    for k=1:14
                        % midpoint
                        % f = gait2dc('Dynamics',(xnew+x)/2, (xnew-x)/dt, u);
                        % BE
                        f = testCase.model.getDynamics(xnew, (xnew-x)/dt, u);
                        sumViol_new = norm((fmin-f < 0) .*(fmin-f) + (fmax-f > 0) .*(fmax-f)); % norm of f if f is out of bounds
                        neval = neval + 1;
                        if (sumViol_new < sumViol)
                            break
                        end
                        dx = dx/2;
                        xnew = xnew - dx;
                    end
                end
                if (j==options.maxiterations)
                    info = -1;
                    return
                end
                xdot = (xnew-x)/dt;
                x = xnew;
                t =  t + dt;
                tout = [tout ; t];
                xout = [xout  x];
            end           
            
        end

        %======================================================================
        %> @brief Helperfunction for muscle stimulation function for passive simulation
        %>
        %> @param testCase
        %> @param t
        %>
        %> @retval u
        %======================================================================
        function [u] = nostim(testCase, t)
            u = zeros(testCase.model.nControls,1);
        end
        
        %======================================================================
        %> @brief Helperfunction for muscle stimulation function for 50% muscle stimulation
        %>
        %> @todo Does it make sense to use 0.5 for all? We should
        %> differentiate between muscle exctitation and torques
        %>
        %> @param testCase
        %> @param t
        %>
        %> @retval u
        %======================================================================
        function [u] = halfstim(testCase, t)
            u = 0.5+zeros(testCase.model.nControls,1);
        end
        
        %-----------------------------------------------------------------
        %> @brief Function to verify the derivative outcomes
        %> 
        %> @param testCase
        %> @param value         Array or Matrix: Derivatives of Model class
        %> @param valueRef      Array or Matrix: Reference derivatives
        %> @param kind          String: Name of derivative (e.g. dfdx)
        %> @param namesDim1     (optional) Cell array with strings: Names defining the first dimension of the value
        %> @param namesDim2     (optional) Cell array with strings: Names defining the second dimension of the value
        %-----------------------------------------------------------------
        function matcompare(testCase,value,valueRef, kind, namesDim1, namesDim2)
            % error checking
            if ~isequal(size(value), size(valueRef))
               error('ModelTest:matcompare', 'a and b must have equal size to be compared'); 
            end
            
            if ~isvector(value) && ~isvector(valueRef) % matrices
                % get maximal error
                [maxerr,irow] = max(abs(value-valueRef)); 
                [maxerr,icol] = max(maxerr);
                irow = irow(icol);
                
                % get message str
                str{1} = sprintf('Max. difference is too high: %14.10f (%14.10f vs. %14.10f)  at ', ...
                    full(maxerr), full(value(irow,icol)),full(valueRef(irow,icol)));
                if nargin > 4 && iscell(namesDim1) && length(namesDim1) == size(value, 1)
                    str{2} = sprintf('''%s''(idx:%d) and', namesDim1{irow}, irow);
                else
                    str{2} = sprintf('%d and', irow);
                end
                if nargin > 5 && iscell(namesDim2) && length(namesDim2) == size(value, 2)
                    str{3} = sprintf('''%s''(idx:%d) \n', namesDim2{icol}, icol);
                else
                    str{3} = sprintf('%d \n', icol);
                end
                str = [kind ': ' strjoin(str)];
                
                % get ratio
                if maxerr == 0
                    ratio = 0;
                else
                    ratio = abs(full(maxerr)/full(valueRef(irow,icol)));
                end
                
            elseif isvector(value) && isvector(valueRef)
                % get maximal error
                [maxerr,irow] = max(abs(value-valueRef)); 
                
                % get message str
                str{1} = sprintf('Max. difference is too high: %14.10f (%14.10f vs. %14.10f)  at ', ...
                    full(maxerr), full(value(irow)),full(valueRef(irow)));
                if nargin > 4 && iscell(namesDim1) && length(namesDim1) == size(value, 1)
                    str{2} = sprintf('''%s''(idx:%d) \n', namesDim1{irow}, irow);
                else
                    str{2} = sprintf('%d \n', irow);
                end
                str = [kind ': ' strjoin(str)];
                
                % get ratio
                if maxerr == 0
                    ratio = 0;
                else
                    ratio = abs(full(maxerr)/full(valueRef(irow)));
                end
                
            else
                error('ModelTest:matcompare: a and b must be a matrix or a vector');
            end
            
            testCase.verifyLessThan( ratio, testCase.TOL_DERIVATIVETEST, str);
            fprintf(strrep(str, ' too high', ''));
        end
        
        
        %-----------------------------------------------------------------
        %> @brief Function to generate a vector with random numbers of a
        %> specific type (e.g. q or u)
        %> 
        %> @details
        %> Usage:
        %> @code
        %> x = testCase.generateRandType('States');
        %> u = testCase.generateRandType('Controls');
        %> q = testCase.generateRandType('q', 1);
        %> M = testCase.generateRandType('torque', 0);
        %> @endcode
        %>
        %> @param  testCase
        %> @param  type          String: Either 'States' or 'Controls' or a type which has to be identical to the type specified in model.states or in model.controls
        %> @param  isState       (optional) Boolean: Needed if type is not 'States' or 'Controls'. If true, the type is in the model.states. Otherwise we search in the model.controls.
        %> @retval randType      Double vector: Vector containing uniformly sampled random numbers of that type 
        %>                       in the interval [xmin, xmax] which is specified in model.states or model.controls
        %-----------------------------------------------------------------
        function randType = generateRandType(testCase,type,isState)
            
            if strcmp(type, 'States') || strcmp(type, 'states')
                % Get lower and upper bound for states
                LB = testCase.model.states.xmin;
                UB = testCase.model.states.xmax;
            elseif strcmp(type, 'Controls') || strcmp(type, 'controls')
                % Get lower and upper bound for controls
                LB = testCase.model.controls.xmin;
                UB = testCase.model.controls.xmax;
            else
                if isState % Search in state table
                    % Get indices in state table
                    idx = testCase.model.extractState(type);
                    % Get lower and upper bound for this type
                    LB = testCase.model.states.xmin(idx);
                    UB = testCase.model.states.xmax(idx);
                else % Search in control table
                    % Get indices in state table
                    idx = testCase.model.extractState(type);
                    % Get lower and upper bound for this type
                    LB = testCase.model.states.xmin(idx);
                    UB = testCase.model.states.xmax(idx);
                end
            end
            
            % Generate uniformly sampled random numbers in the interval [LB, UB]
            randType = LB + (UB - LB) .* rand(size(UB));
           
        end
        
    end
    
    methods(Static, Access = protected)
        
        %======================================================================
        %> @brief Function to solves f(x)=0 using Newton's method.
        %>
        %> @param fun
        %> @param x0
        %> @param options
        %>
        %> @retval x
        %> @retval f
        %> @retval info
        %> @retval iterations
        %======================================================================
        function [x,f,info, iterations] = fsolve1(fun, x0, options)
            % solves f(x)=0 using Newton's method
            
            % initialize
            info = 0;
            iterations = 0;
            x = x0;
            [f, J] = fun(x);
            bestfnorm = norm(f);
            
            % Start loop
            while (1)
                iterations = iterations + 1;
                % check if max. number of iterations exceeded
                if (iterations > options.maxiterations)
                    info = 1;
                    return
                end
                
                % compute the (Gauss-)Newton step dx
                % dx = -(J'*J)\(J'*f);
                dx = -J\f;
                
                % do the step dx
                x = x + dx;
                
                % Make sure we have improvement in the merit function norm(f)
                [f,J] = fun(x);
                fnorm = norm(f);
                while (fnorm > bestfnorm)
                    iterations = iterations + 1;
                    % check if max. number of iterations exceeded
                    if (iterations > options.maxiterations)
                        info = 1;
                        return
                    end
                    dx = dx/2;
                    x = x-dx;
                    [f,J] = fun(x);
                    fnorm = norm(f);
                end
                bestfnorm = fnorm;
                
                % are we done?
                if fnorm < options.ftol
                    info = 0;
                    return;
                end
                
            end
            
        end
        
        
    end
    
    
    
end
