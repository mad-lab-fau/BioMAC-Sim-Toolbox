% ======================================================================
%> @file @Collocation_Exo/Collocation_Exo.m
%> @brief Matlab class describing collocation using the exoskeleton
%>
%> @author Anne
%> @date October, 2019
% ======================================================================

% ======================================================================
%> @brief The class for defining an optimal control of musculoskeletal
%>
%> @details
%> - Inherits from Collocation
%> - Adds some dependency on the duration
% ======================================================================
classdef Collocation_Exo < Collocation
    
    properties
        
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
        function [obj] = Collocation_Exo(varargin)
            
            % call Gait2dc constructor
            obj = obj@Collocation(varargin{:});
            
        end
    end
end
      
