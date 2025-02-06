% ======================================================================
%> @file Collocation_NN.m
%> @brief Class for python import 
%
%> @author Johannes Derrer
%> @date November,2021
% ======================================================================

% ======================================================================
%> @brief The class
% ======================================================================
classdef Collocation_NN < Collocation
    
    properties
        %%%%%%% for python import for NN
        net;
        np;
    end
    
    methods
        %======================================================================
        %> @brief Constructor setting default Collocation options.
        %>
        %> @param   model     Model class object
        %> @param   nNodes    (optional) Double: Number of collocation nodes 
        %> @param   Euler     (optional) String: Discretization method: backward euler 'BE' or midpoint euler 'ME' 
        %> @param   logfile   (optional) String: Logfile including the path. The input can be skipped with an empty string.
        %> @param   plotLog   (optional) Bool: If true, plot log in call of Collocation.objfun() (stick figure, objective, constrains) (default: 0)
        %>
        %> @retval    obj   Collocation class object
        %======================================================================
        function [obj] = Collocation_NN(model, nNodes, Euler, logfile, plotLog)
            if nargin < 5
                plotLog = 0;
            end
            obj = obj@Collocation(model, nNodes, Euler, logfile, plotLog);

            %%%%%%% for python import for NN

            obj.np = py.importlib.import_module("numpy");
            modelFolder = "C:\Users\annek\Documents\MATLAB\Bio-Sim-Toolbox\rnn_test_model\";
            module = py.importlib.import_module("keras.models");
            obj.net = module.load_model(modelFolder);
        end
    end
end
