classdef Gait3d_smoothsphere < Gait3d
    
    
    methods
        function obj = Gait3d_smoothsphere(varargin)
            % Call superclass constructor to initialize standard properties
            obj = obj@Gait3d(varargin{:});
            
            % Build custom MEX for smoothsphere
            nameMEX = Gait3d_smoothsphere.getMexFiles(obj.osim.name);
            obj.hdlMEX = str2func(nameMEX);
            
            % Re-initialize MEX state with custom MEX
            obj.initMex();
        end
    end
    
    methods (Static)
        name_MEX = getMexFiles(osimName, rebuild_contact, debugMode)
    end
end
