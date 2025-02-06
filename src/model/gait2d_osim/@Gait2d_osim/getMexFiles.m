%======================================================================
%> @file @Gait2d_osim/getMexFiles.m
%> @brief Gait2d_osim function to make MEX functions
%> @details
%> Details: Gait2d_osim::getMexFiles()
%> The code is copied from gait3d entirely, with adjustments for gait2d_osim
%>
%> @author Eva Dorschky, Ton van den Bogert, Marlies Nitschke
%> @date May, 2022
%======================================================================

%======================================================================
%> @brief Function to make MEX functions
%> @static
%> @public
%>
%> @details
%> Generates Autolev source file for model and builds the MEX function
%>
%> @param  osimName             String: Name of Opensim model. Currently, the following is supported:
%>                                - 'gait10dof18musc' and 'gait2d_osim'
%> @param  rebuild_contact      (optional) Bool: If true, enforce rebuilding contact_raw.c (default: 0)
%> @param  debugMode            (optional) Bool: If true, code will be compiled in debug mode and not optimized. (default: 0)
%> @retval name_MEX             String: Name of the mex function which has to be called
%======================================================================
function name_MEX = getMexFiles(osimName, rebuild_contact, debugMode)

clear mex						% to make sure mex files are not loaded, so we can delete/overwrite

% define names for Autolev code and c-code based on the model name
switch osimName
    case {'gait10dof18musc','gait10dof18musc.osim','gait10dof18musc_wrap.osim'}
        name_MEX = 'gait10dof18musc';
        warning('This model is not validated.')
    case {'gait10dof18musc_ll','gait10dof18musc_ll.osim','gait10dof18musc_wrap_ll.osim','gait2d.osim'}
        name_MEX = 'gait2d_osim';
    otherwise
        error('Gait2d_osim:getMexFiles','C files do not exist for OpenSim Model %s',osimName);
end

if isdeployed % don't compile if function is running in deployed mode (e.g. on HPC cluster)
    return
end

if nargin < 3
    rebuild_contact = 0;
end

if nargin > 3 && debugMode
    optimize = 0;
    mexoptOpt = '-g';
else
    optimize = 1;
end

% change directory
currentFolder = pwd;
methodFolder = strrep(mfilename('fullpath'), [filesep , mfilename], ''); %.../Gait2d_osim
cfilesFolder = strrep(methodFolder, '@Gait2d_osim', 'c_files');          %.../c_files
cd(cfilesFolder);

% get ending of object code based on the system
if ispc % true for windows
    objext = 'obj';
elseif isunix % true for Linux and Mac
    objext = 'o';
    mexoptUNIX = 'GCC=''/usr/bin/gcc''';
else
    error('System is unknown');
end

% define names of Autolev-files, c-files, mex-files, etc.
name_model_path        = [name_MEX filesep name_MEX];
name_model_c           = [name_MEX filesep name_MEX '.c'];
name_model_mex         = [name_MEX filesep name_MEX '.' mexext];
name_create_alfile_m   = [name_MEX filesep 'create_alfile_' name_MEX '.m'];
hdl_create_alfile      = str2func(['@create_alfile_' name_MEX]);
name_alfile            = [name_MEX filesep name_MEX '.al'];
name_alfile_raw_c      = [name_MEX filesep name_MEX '_raw.c'];
name_alfile_raw_in     = [name_MEX filesep name_MEX '_raw.in'];
name_alfile_c          = [name_MEX filesep name_MEX '_al.c'];
name_alfile_o          = [name_MEX filesep name_MEX '_al.' objext];
if strcmp(name_MEX, 'gait2d_osim') || strcmp(name_MEX, 'gait10dof18musc') % For this model, separate FK files% For this model, separate FK files
    name_alfile_FK            = [name_MEX filesep name_MEX '_FK.al'];
    name_alfile_FK_raw_c      = [name_MEX filesep name_MEX '_FK_raw.c'];
    name_alfile_FK_raw_in     = [name_MEX filesep name_MEX '_FK_raw.in'];
    name_alfile_FK_c          = [name_MEX filesep name_MEX '_FK_al.c'];
    name_alfile_FK_o          = [name_MEX filesep name_MEX '_FK_al.' objext];
    name_alfile_DynNoDer        = [name_MEX filesep name_MEX '_NoDer.al'];
    name_alfile_DynNoDer_raw_c  = [name_MEX filesep name_MEX '_NoDer_raw.c'];
    name_alfile_DynNoDer_raw_in = [name_MEX filesep name_MEX '_NoDer_raw.in'];
    name_alfile_DynNoDer_c      = [name_MEX filesep name_MEX '_NoDer_al.c'];
    name_alfile_DynNoDer_o      = [name_MEX filesep name_MEX '_NoDer_al.' objext];
end
name_clean_alfile      = ['clean_' name_MEX];
name_clean_alfile_c    = [name_MEX filesep name_clean_alfile '.c'];
name_clean_alfile_mex  = [name_MEX filesep name_clean_alfile '.' mexext];
name_header_multibody  = [name_MEX filesep name_MEX '_multibody.h'];
name_contact       = 'contact_2d';
folder_contact         = ['contact' filesep];
name_contact_alfile    = [folder_contact name_contact '.al'];
name_contact_alfile_c  = [folder_contact 'contact_2d_al.c'];
name_contact_raw_c     = [folder_contact 'contact_2d_raw.c'];
name_contact_raw_in    = [folder_contact 'contact_2d_raw.in'];
name_contact_alfile_o  = [folder_contact 'contact_2d_al.' objext];
name_clean_contact     = 'clean_contact_2d';
name_clean_contact_c   = [folder_contact name_clean_contact '.c'];
name_clean_contact_mex = [folder_contact name_clean_contact '.' mexext];
name_header_contact    = [folder_contact 'gait2d_osim_contact.h'];


% if there is an input argument, just clean up
if nargin > 1 && clean == 1
    fprintf('Cleaning up...');
    warning('off', 'MATLAB:DELETE:FileNotFound');	% don't warn me if the file does not exist
    delete('alTmp.*');
    delete(name_model_mex);
    delete(name_alfile);
    delete(name_alfile_raw_c);
    delete(name_alfile_raw_in);
    delete(name_alfile_c);
    delete(name_alfile_o);
    if strcmp(name_MEX, 'gait2d_osim') || strcmp(name_MEX, 'gait10dof18musc')
        delete(name_alfile_FK);
        delete(name_alfile_FK_raw_c);
        delete(name_alfile_FK_raw_in);
        delete(name_alfile_FK_c);
        delete(name_alfile_FK_o);
        delete(name_alfile_DynNoDer);
        delete(name_alfile_DynNoDer_raw_c);
        delete(name_alfile_DynNoDer_raw_in);
        delete(name_alfile_DynNoDer_c);
        delete(name_alfile_DynNoDer_o);
    end
    delete(name_clean_alfile_mex);
    delete(name_contact_raw_c);
    delete(name_contact_raw_in);
    delete(name_contact_alfile_c);
    delete(name_contact_alfile_o);
    delete(name_clean_contact_mex);
    fprintf('Done.\n');
    warning('on', 'MATLAB:DELETE:FileNotFound');	% turn the warnings back on
    cd(currentFolder); % change back
    return
end

% If the Opensim model has changed since we changed the according autolev code, give a warning
% Eventually, we will just read the most recent Opensim model and generate the gait2d_osim.al file from that
if exist([name_MEX, '.m'],'file')
    if mdate(osimName) > mdate([name_MEX,'.m'])
        disp('Warning: the osim file has changed since you last changed this program.');
        disp('You may need to change create_alfile.m and/or associated files.');
        disp('Hit ENTER if you proceed with model build, or CTRL-C to quit');
        pause
    end
end

%Commented out because Autolev is not running on Windows
%Run Autolev, if needed
if needs_rebuild(name_alfile_raw_c, {name_create_alfile_m})
   hdl_create_alfile(name_alfile);
   runautolev(name_model_path, 12);		% second input is estimated completion time in minutes
end

if rebuild_contact || needs_rebuild(name_contact_raw_c, {name_contact_alfile})
   cd(folder_contact); % Why do we need this, but do not need it to run autolev for the model?
   runautolev(name_contact, 1);
   cd('../')
end

% Compile the C code generated by Autolev if needed
if needs_rebuild(name_alfile_o,{name_alfile_raw_c, name_clean_alfile_c, name_header_multibody})
    fprintf('Recompiling %s...\n', name_clean_alfile);
    mex(name_clean_alfile_c, '-outdir', name_MEX);
    fprintf('Done.\n');
    fprintf('Cleaning the C code...\n');
    feval(name_clean_alfile);
    fprintf('Done.\n');
    if optimize
        if strcmp(mexext,'mexw64')
            completion_time = 'unknown time';
        else
            completion_time = datestr(addtodate(now, 10, 'minute'));
        end
    else
        if strcmp(mexext,'mexw64')
            completion_time = datestr(addtodate(now, 10, 'minute'));
        else
            completion_time = datestr(addtodate(now, 10, 'minute'));
        end
    end
    fprintf('Compiling C code, estimated completion at %s...', completion_time);
    if isunix && ~optimize
        mex(mexoptUNIX, mexoptOpt,'-largeArrayDims','-c',name_alfile_c, '-outdir', name_MEX);
    elseif isunix && optimize
        mex(mexoptUNIX,'-largeArrayDims','-c',name_alfile_c, '-outdir', name_MEX);
    elseif ~isunix && ~optimize
        mex(mexoptOpt,'-largeArrayDims','-c',name_alfile_c, '-outdir', name_MEX);
    elseif ~isunix && optimize
        mex('-largeArrayDims','-c',name_alfile_c, '-outdir', name_MEX);
    end
    if strcmp(name_MEX, 'gait2d_osim') || strcmp(name_MEX, 'gait10dof18musc') % For this model, separate FK files
        if isunix && ~optimize
            mex(mexoptUNIX, mexoptOpt,'-largeArrayDims','-c',name_alfile_FK_c, '-outdir', name_MEX);
            mex(mexoptUNIX, mexoptOpt,'-largeArrayDims','-c',name_alfile_DynNoDer_c, '-outdir', name_MEX);
        elseif isunix && optimize
            mex(mexoptUNIX,'-largeArrayDims','-c',name_alfile_FK_c, '-outdir', name_MEX);
            mex(mexoptUNIX,'-largeArrayDims','-c',name_alfile_DynNoDer_c, '-outdir', name_MEX);
        elseif ~isunix && ~optimize
            mex(mexoptOpt,'-largeArrayDims','-c',name_alfile_FK_c, '-outdir', name_MEX);
            mex(mexoptOpt,'-largeArrayDims','-c',name_alfile_DynNoDer_c, '-outdir', name_MEX);
        elseif ~isunix && optimize
            mex('-largeArrayDims','-c',name_alfile_FK_c, '-outdir', name_MEX);
            mex('-largeArrayDims','-c',name_alfile_DynNoDer_c, '-outdir', name_MEX);
        end
    end
    fprintf('Done.\n');
end

% Rebuild the contact model C code if needed
if needs_rebuild(name_contact_alfile_o, {name_contact_raw_c, name_clean_contact_c, name_header_contact})
    fprintf('Recompiling %s...', name_clean_contact);
    mex(name_clean_contact_c, '-outdir', folder_contact);
    fprintf('Done.\n');
    fprintf('Cleaning the C code...');
    feval(name_clean_contact);
    fprintf('Done.\n');
    fprintf('Compiling contact_al.c...');
    if isunix && ~optimize
        mex(mexoptUNIX, mexoptOpt,'-largeArrayDims','-c',name_contact_alfile_c, '-outdir', folder_contact);
    elseif isunix && optimize
        mex(mexoptUNIX,'-largeArrayDims','-c',name_contact_alfile_c, '-outdir', folder_contact);
    elseif ~isunix && ~optimize
        mex(mexoptOpt,'-largeArrayDims','-c',name_contact_alfile_c, '-outdir', folder_contact);
    elseif ~isunix && optimize
        mex('-largeArrayDims','-c',name_contact_alfile_c, '-outdir', folder_contact);
    end
    fprintf('Done.');
end

% Rebuild the MEX binary if needed
if needs_rebuild(name_model_mex, {name_model_c, name_alfile_o, name_contact_alfile_o, name_header_contact, name_header_multibody})
    fprintf('Compiling %s...',name_model_c);
    if strcmp(name_MEX, 'gait2d_osim') || strcmp(name_MEX, 'gait10dof18musc')% For this model, separate FK files
        if isunix && ~optimize
            mex(mexoptUNIX, mexoptOpt,'-largeArrayDims',name_model_c,name_alfile_o,name_alfile_DynNoDer_o,name_alfile_FK_o,name_contact_alfile_o, '-outdir', name_MEX);
        elseif isunix && optimize

            
            mex(mexoptUNIX,'-largeArrayDims',name_model_c,name_alfile_o,name_alfile_DynNoDer_o,name_alfile_FK_o,name_contact_alfile_o, '-outdir', name_MEX);
        elseif ~isunix && ~optimize
            mex(mexoptOpt,'-largeArrayDims',name_model_c,name_alfile_o,name_alfile_DynNoDer_o,name_alfile_FK_o,name_contact_alfile_o, '-outdir', name_MEX);
        elseif ~isunix && optimize
            mex('-largeArrayDims',name_model_c,name_alfile_o,name_alfile_DynNoDer_o,name_alfile_FK_o,name_contact_alfile_o, '-outdir', name_MEX);
        end
    else
        if isunix && ~optimize
            mex(mexoptUNIX, mexoptOpt,'-largeArrayDims',name_model_c,name_alfile_o,name_contact_alfile_o, '-outdir', name_MEX);
        elseif isunix && optimize
            mex(mexoptUNIX,'-largeArrayDims',name_model_c,name_alfile_o,name_contact_alfile_o, '-outdir', name_MEX);
        elseif ~isunix && ~optimize
            mex(mexoptOpt,'-largeArrayDims',name_model_c,name_alfile_o,name_contact_alfile_o, '-outdir', name_MEX);
        elseif ~isunix && optimize
            mex('-largeArrayDims',name_model_c,name_alfile_o,name_contact_alfile_o, '-outdir', name_MEX);
        end
    end
    fprintf('Done.\n');
end

cd(currentFolder); % change back



%> @cond DO_NOT_DOCUMENT
    function [rebuild] = needs_rebuild(file, dependencies)
        % check if the file needs to be rebuilt, based on dependencies
        %
        % Input:
        %	file..................(string) Name of file we are checking for
        %   dependencies..........(cell array of strings) files that this new file depends on
        %
        % Output:
        %	rebuild.........(logical) true if file must be rebuilt
        
        % if file does not exist, it has to be rebuilt
        if ~exist(file)
            rebuild = true;
            return
        end
        
        % if any of the parent files have changed, rebuild
        rebuild = false;
        for i = 1:numel(dependencies)
            if mdate(dependencies{i}) > mdate(file)
                rebuild = true;
                fprintf('\nRebuilding %s because of changes in %s\n', file, dependencies{i});
            end
        end
    end

    function filedate = mdate(filename)
        % determines the date when the file was last modified
        if ~exist(filename)
            error(['filedate: ' filename ' does not exist.']);
        end
        file_info = dir(which(filename));
        filedate = file_info.datenum;
    end

    function runautolev(filename, minutes)

        % Autolev location depends on computer
        [~,computername] = system('hostname');
        if strfind(computername,'LRI-102855')
            autolev = '/home/Public/Autolev/al.exe';%'C:\Program Files\Autolev\al.exe';
        elseif strfind(computername,'mad-simulation')
            autolev = '/home/Public/Autolev/al.exe';
        else
            autolev = [];
        end
        
        filename_out = filename;
            
        % Run Autolev
        if ~exist(autolev)
            disp(['Error: Autolev needs to generate ' filename_out '_raw.c, but you do not have Autolev installed.']);
            disp('Ask Ton to run this part of the build process.');
            disp(['Hit ENTER to continue with existing ' filename_out '_raw.c, or CTRL-C to quit']);
            pause
        else
            warning('off', 'MATLAB:DELETE:FileNotFound');	% don't warn me if the file does not exist
            delete([filename_out '_raw.c']);		% we need to delete the .c and .in files so Autolev won't ask for overwrite permission
            delete([filename_out '_FK_raw.c']);		% we need to delete the .c and .in files so Autolev won't ask for overwrite permission
            delete([filename_out '_NoDer_raw.c']);		% we need to delete the .c and .in files so Autolev won't ask for overwrite permission
            delete([filename_out '_raw.in']);
            delete([filename_out '_FK_raw.in']);		% we need to delete the .c and .in files so Autolev won't ask for overwrite permission
            delete([filename_out '_NoDer_raw.in']);		% we need to delete the .c and .in files so Autolev won't ask for overwrite permission
            warning('on', 'MATLAB:DELETE:FileNotFound');
            
            fprintf('Autolev is generating %s\n',[filename_out '_raw.c']);
            fprintf('Estimated completion at %s...\n', datestr(addtodate(now, minutes, 'minute')));
            
            if ispc
                system(['"' autolev '" ' filename '.al > nul']);		% double quotes needed because autolev has spaces in its path
            else
                system(['wine "' autolev '" ' filename '.al > nul']);		% double quotes needed because autolev has spaces in its path
            end
            % check if Autolev was successful
            if ~exist([filename_out '_raw.c'])
                warning(['Autolev error in ' filename_out '.al']);
            end
            delete([filename_out '_raw.in']);
            fprintf('Done.\n');
        end
    end
%> @endcond

end
