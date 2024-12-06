%======================================================================
%> @file @Gait2d_osim/scaleOsim.m
%> @brief Gait2d_osim function to scale OpenSim model in Matlab
%> @details
%> Details: Gait2d_osim::scaleOsim()
%>
%> @author Markus Gambietz
%> @date October, 2022
%======================================================================

%======================================================================
%> @brief Function to scale OpenSim model
%> @protected
%> @static
%>
%> @details
%> Saves the model struct using strrep(osimfile, '.osim', '.mat').
%>
%> @param  opensimfile     String: Filename of opensim file including path and extension
%> @param   scale_factors   String, List or Float: .trc marker
%scalefile, list of scalefactors: [nSegments = 14] or [nSegments = 14 * 3],
%(3d case); or float: bodyheight
%> @param   bodyweight  Float: bodyweight in kg
%> @retval opensimfile_tmp     String: Path to the tmp_osim model
%> @retval hash String: Hash to check whether the model was created before
%> @todo: .trc/xml scaling; load scalefactors from table?
%======================================================================
function [opensimfile_tmp, hash] = scaleOsim(opensimfile, scale_factors, bodymass, varargin)

%% Create hash
import java.security.*;
import java.math.*;
import org.opensim.modeling.*

md = MessageDigest.getInstance('MD5');
% Use *10^4 so that close numbers don't create the same hashes
if isscalar(scale_factors) || ~isa(scale_factors,'double')
    dig = md.digest(double(10000*bodymass)) + md.digest(10000*double(scale_factors)); %
else
    % Issue: The hash generation for r does not properly work
    warning('scaleOsim: hash generation for matrices is broken')
    r = strjoin(string(scale_factors));
    r = strrep(strrep(r,' ',''),'.','');
    dig = md.digest(double(10000*bodymass)) + md.digest(0.001*double(r));
end
bi = BigInteger(1, dig);
hash = char(bi.toString(16));

%% Various in args
if nargin > 3
    % Overwrite the hash with a modelname (be careful when using this xD)
    hash = varargin{1};
end



%% Check if the hash already exists, if so, return the tmp model already
% tmp models are saved to a temporary folder in case many are needed
opensimfile_tmp = strrep(which(opensimfile),'.osim',strcat('_tmp/',hash,'.osim'));
if isfile(opensimfile_tmp)
    return;
end
if ~isfolder(strrep(which(opensimfile),'.osim','_tmp/'))
    mkdir(strrep(which(opensimfile),'.osim','_tmp/'));
end


%% Depending on if it is a file, do manual or automatic scaling
osim_model = Gait2d_osim.readOsim(opensimfile);
switch class(scale_factors)
    % Manual scaling in case doubles are used
    case 'double'
        % If boyheight is given, convert scalefactor to sf.
        if isscalar(scale_factors)
            scale_factors = ones(height(osim_model.segments),3)*scale_factors / 1.8;
            % build a matrix
            % If the scalefactors matrix just has one dimension
        elseif size(scale_factors,2) == 1
            scale_factors = repmat(scale_factors,1,3);
        end

        %% Initialize
        tmp_mod = Model(which(opensimfile));
        tmp_mod.initSystem;

        %% Setup scale Tool
        scaleTool = ScaleTool(which('setup_manualscale.xml'));
        scaleTool.setSubjectMass(bodymass);
        scaleTool.getModelScaler().setOutputModelFileName(opensimfile_tmp);
        scaleTool.getModelScaler().setScaleSetFile(which('scaleset.xml'));
        %% add the scales
        for i = 1:size(osim_model.segments.Properties.RowNames) % Loop over segments
            skale = Scale();
            skale.setSegmentName(osim_model.segments.Properties.RowNames{i});
            sf = Vec3(scale_factors(i,1),scale_factors(i,2),scale_factors(i,3));
            skale.setScaleFactors(sf);
            scaleTool.getModelScaler().addScale(skale)
        end

        %% run the scale tool & call the model
        scaleTool.getModelScaler().processModel(tmp_mod,which(opensimfile),bodymass);

    case {'char', 'string'}
        scale_factors = char(scale_factors);
        %% Scale tool scaling: mostly copied from my MT
        if strcmp(scale_factors(end-3:end),'.trc')
            inputOSIM = which(opensimfile);
            myModel= Model(inputOSIM);
        	myModel.initSystem;

        	%set settings and files;
        	scaleTool = ScaleTool(which('Setup_Scaling_Shoulder_2.xml'));
        	timeArray = ArrayDouble();
        	%scale the first 0.3 seconds;
        	timeArray.set(0,0);
        	timeArray.set(1,0.15);
        	scaleTool.setSubjectMass(bodymass);
            scaleTool.getModelScaler().setApply(1);
        	scaleTool.getModelScaler().setMarkerFileName(which(scale_factors));
        	scaleTool.getModelScaler().setTimeRange(timeArray);
        	scaleTool.getModelScaler().setOutputModelFileName(opensimfile_tmp);
        	scaleTool.getModelScaler().setOutputScaleFileName(strcat(opensimfile_tmp,'.xml'));
        	%Run scaling;
        	scaleTool.getModelScaler().processModel(myModel,'',bodymass)
        	%Marker placer on the scaled model;
        	myModel= Model(opensimfile_tmp);
        	myModel.initSystem();
        	%Repeat the settings from the scale tool;
        	timeArray.set(0, 0);
        	timeArray.set(1, 0.15);
        	scaleTool.getMarkerPlacer().setApply(1);
        	scaleTool.getMarkerPlacer().setTimeRange(timeArray);
        	scaleTool.getMarkerPlacer().setOutputModelFileName(opensimfile_tmp);
        	scaleTool.getMarkerPlacer().setStaticPoseFileName(which(scale_factors));
        	scaleTool.getMarkerPlacer().setOutputMotionFileName(strcat(opensimfile_tmp,'.mot'));
        	%Run marker placer
        	scaleTool.getMarkerPlacer().processModel(myModel,'')

        elseif strcmp(scale_factors(end-3:end),'.xml')
            %% Initialize
            tmp_mod = Model(which(opensimfile));
            tmp_mod.initSystem;

            %% Setup scale Tool
            scaleTool = ScaleTool(which('setup_manualscale.xml'));
            scaleTool.setSubjectMass(bodymass);
            scaleTool.getModelScaler().setOutputModelFileName(opensimfile_tmp);
            scaleTool.getModelScaler().setScaleSetFile(which(scale_factors));
            %% run the scale tool & call the model
            scaleTool.getModelScaler().processModel(tmp_mod,which(opensimfile),bodymass)
        else
            error('Gait2d_osim.scaleOsim: Unknown Filetype')
        end
end