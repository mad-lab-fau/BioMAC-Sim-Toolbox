%======================================================================
%> @file @Gait2dc/getMexFiles.m
%> @brief Gait2dc function to make MEX functions
%> @details
%> Details: Gait2dc::getMexFiles()
%>
%> @author Ton van den Bogert
%> @date January, 2018
%======================================================================

%======================================================================
%> @brief Function to make MEX functions
%> @static
%> @public
%>
%> @details
%> Generates Autolev source file for model and builds the MEX function
%>
%> @param  mexUpdate  (optional) Bool: If true, MEX will be updated using
%>                    the Makefile (default: 0)
%======================================================================
function getMexFiles(mexUpdate)

clear mex						% to make sure mex files are not loaded, so we can delete/overwrite

% change directory
currentFolder = pwd;
methodFolder = strrep(mfilename('fullpath'), [filesep , mfilename], '');
cfilesFolder = strrep(methodFolder, '@Gait2dc', 'c_files');
cd(cfilesFolder);

if isdeployed % don't compile if function is running in deployed mode (e.g. on HPC cluster)
    return
end

if nargin > 0 && mexUpdate % update the mex file
    
    if ismac
        if (system('make -f Makefile.mac') ~= 0)
            error('Error in compilation of MEX function.');
        end
    elseif ispc
        if (system('@echo off & call "vcvars32" & nmake /nologo /F Makefile.win') ~= 0)
            %cd ..
            error('Error in compilation of MEX function.');
        end
    else
        error('Ubuntu is currently not supported');
    end
    
else  % create new mex
    mexfile = ['gait2dc.' mexext];
    % Try to run Autolev if the raw C files are outdated
    if needs_rebuild('contact_al_raw.c',{'contact.al'})
        runautolev('contact');
    end
    if needs_rebuild('gait2dc_dynamics_al_raw.c',{'gait2dc_dynamics.al', 'gait2dc_FK.al'})
        runautolev('gait2dc_dynamics');
    end
    if needs_rebuild('gait2dc_stick_al_raw.c',{'gait2dc_stick.al', 'gait2dc_FK.al'})
        runautolev('gait2dc_stick');
    end
    
    % Clean the raw C files
    % They are in the repository, for developers who do not have Autolev
    if needs_rebuild('autolevclean.c',{}) || needs_rebuild('contact_al.c',{'contact_al_raw.c'}) || needs_rebuild('gait2dc_dynamics_al.c',{'gait2dc_dynamics_al_raw.c'}) || needs_rebuild('gait2dc_FK_al.c',{'gait2dc_FK_al_raw.c'})
        mex autolevclean.c
        autolevclean;
        delete(['autolevclean.' mexext]);		% clean up
    end
    
    % Compile and link the C code
    if needs_rebuild(mexfile, {'gait2dc.c','gait2dc_dynamics_al_raw.c','gait2dc_FK_al_raw.c', 'contact_al_raw.c', 'acc_al.c'})
        disp('Compiling C code, this takes approx. 2 minutes...');
        mex -largeArrayDims gait2dc.c gait2dc_dynamics_al.c gait2dc_FK_al.c gait2dc_stick_al.c contact_al.c acc_al.c
        
        % Completion message
        disp(['gait2dc.' mexext ' is ready.']);
    end
    
end

cd(currentFolder); % change back



%> @cond DO_NOT_DOCUMENT
    function runautolev(filename)
       
        % Autolev location depends on computer
        [~,computername] = system('hostname');
        if strfind(computername,'LRI-102855')
            autolev = '/home/Public/Autolev/al.exe';
        else
            autolev = [];
        end
            
        % Run Autolev
        if ~exist(autolev)
            disp(['Error: Autolev needs to generate ' filename '_raw.c, but you do not have Autolev installed.']);
            disp('Ask Ton to run this part of the build process.');
            disp(['Hit ENTER to continue with existing ' filename '_al_raw.c, or CTRL-C to quit']);
            pause
        else
            warning('off', 'MATLAB:DELETE:FileNotFound');	% don't warn me if the file does not exist
            delete([filename '_al_raw.c']);		% we need to delete the .c and .in files so Autolev won't ask for overwrite permission
            delete([filename '_al_raw.in']);
            warning('on', 'MATLAB:DELETE:FileNotFound');
            fprintf('Autolev is generating %s...\n',[filename '_raw.c']);
            
            if ispc
                system(['"' autolev '" ' filename '.al > nul']);		% double quotes needed because autolev has spaces in its path
            else
                system(['wine "' autolev '" ' filename '.al > nul']);		% double quotes needed because autolev has spaces in its path
            end            
            
            % check if Autolev was successful
            if ~exist([filename '_al_raw.c'])
                error(['Autolev error in ' filename '.al']);
            end
            delete([filename '_al_raw.in']);
            fprintf('Done.\n');
        end
    end

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
        file_info = dir(filename);
        filedate = file_info.datenum;
    end
%> @endcond

end
