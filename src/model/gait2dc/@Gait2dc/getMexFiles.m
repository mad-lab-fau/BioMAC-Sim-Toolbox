%======================================================================
%> @file @Gait2dc/getMexFiles.m
%> @brief Gait2dc function to make MEX functions
%> @details
%> Details: Gait2dc::getMexFiles()
%>
%> @author Ton van den Bogert
%> @date January, 2025
%======================================================================

%======================================================================
%> @brief Function to make MEX functions
%> @static
%> @public
%>
%> @details
%> Builds the MEX function from Autolev source code, and from gait2dc.c MEX wrapper.
%======================================================================
function getMexFiles()

    clear mex % to make sure mex files are not loaded, so we can delete/overwrite if needed
    
    % change directory
    currentFolder = pwd;
    methodFolder = strrep(mfilename('fullpath'), [filesep , mfilename], '');
    cfilesFolder = strrep(methodFolder, '@Gait2dc', 'c_files');
    cd(cfilesFolder);
    
    if isdeployed % don't compile if function is running in deployed mode (e.g. on HPC cluster)
        return
    end
    
    % Run Autolev (if installed) to rebuild any C files that are outdated
    % (Note, we don't need to run gait2dc_FK.al here, it will be done within
    % gait2dc_dyn.al)
    clean_needed = false;
    if needs_rebuild('contact_al.c',{'contact.al', 'extforce.a'})
        runautolev('contact');
        clean_needed = true;
    end

    if needs_rebuild('gait2dc_dyn_al.c',{'gait2dc_dyn.al', 'gait2dc_FK.al', 'gait2dc_kin.al'})
        runautolev('gait2dc_dyn');
        clean_needed = true;
    end
    if needs_rebuild('gait2dc_acc_al.c',{'gait2dc_acc.al', 'accelerometer.a', 'gait2dc_kin.al'})
        runautolev('gait2dc_acc');
        clean_needed = true;
    end
    if needs_rebuild('gait2dc_stick_al.c',{'gait2dc_stick.al', 'gait2dc_kin.al'})
        runautolev('gait2dc_stick');
        clean_needed = true;
    end
    
    % Clean the raw C files
    % They are in the repository, for developers who do not have Autolev
    if (clean_needed)
        fprintf('Cleaning the raw C code.\n')
        mex autolevclean.c  % compile the cleaning code into a MEX binary
        autolevclean;
        delete(['autolevclean.' mexext]);		% remove the MEX binary
        fprintf('Done.\n');

        % after generating the clean C code, the raw C code can be deleted,
        % and also the .in files
        delete('*_al_raw.c');
        delete('*_al_raw.in');
    end
    
    % Compile and link the C code for the gait2dc mex function
    mexfile = ['gait2dc.' mexext];
    if needs_rebuild(mexfile, {'gait2dc.c','gait2dc_dyn_al.c','gait2dc_FK_al.c', 'gait2dc_stick_al.c', 'contact_al.c', 'gait2dc_acc_al.c'})
        disp('Compiling C code, this takes approx. 30 seconds...');
        % mex -largeArrayDims gait2dc.c gait2dc_dynamics_al.c gait2dc_FK_al.c gait2dc_stick_al.c contact_al.c acc_al.c
        mex gait2dc.c gait2dc_dyn_al.c gait2dc_FK_al.c gait2dc_stick_al.c contact_al.c gait2dc_acc_al.c
        
        % Completion message
        disp(['gait2dc.' mexext ' is ready.']);
    end
    
    cd(currentFolder); % change back

end

%> @cond DO_NOT_DOCUMENT
function runautolev(filename)

    % look for Autolev in a few known places
    places = {'C:\Users\Ton\Documents\Software\Autolev\al.exe', '/home/Public/Autolev/al.exe'};
    autolev = [];
    for i = 1:numel(places)
        if exist(places{i})
            autolev = places{i};
            break;
        end
    end
              
    % Run Autolev
    if ~exist(autolev)
        disp(['Error: Autolev needs to generate a new version of ' filename '_al.c, but you do not have Autolev installed.']);
        disp('Ask Ton to run this part of the build process.');
        disp(['Hit ENTER to continue with existing ' filename '_al.c, or CTRL-C to quit']);
        pause
    else
        warning('off', 'MATLAB:DELETE:FileNotFound');	% don't warn me if the file does not exist
        delete([filename '_al_raw.c']);		% we need to delete some temporary files, if they exist, so Autolev won't ask for overwrite permission
        delete([filename '_al_raw.in']);
        warning('on', 'MATLAB:DELETE:FileNotFound');
        fprintf('Autolev is generating %s...\n',[filename '_al_raw.c']);
        
        if ispc
            system(['"' autolev '" ' filename '.al > nul']);		% double quotes needed because autolev has spaces in its path
        else
            system(['wine "' autolev '" ' filename '.al > nul']);		% double quotes needed because autolev has spaces in its path
        end            
        
        % check if Autolev was successful
        if ~exist([filename '_al_raw.c'])
            error(['Autolev error in ' filename '.al']);
        end
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
    
    % if file does not exist in the current directory, it has to be rebuilt
    if ~exist(['./' file])
        rebuild = true;
        fprintf('\nRebuilding %s because it does not exist.\n', file);
        return
    end
    
    % if any of the parent files have changed, rebuild
    rebuild = false;
    for i = 1:numel(dependencies)
        if mdate(dependencies{i}) > mdate(file)
            rebuild = true;
            fprintf('\nRebuilding %s because of changes in one of these: %s', file);
            for j = 1:numel(dependencies)
                fprintf(' %s',dependencies{j});
            end
            fprintf('\n');
            return
        end
    end
end

function filedate = mdate(filename)
    % determines the date when the file was last modified
    file_info = dir(filename);
    if numel(file_info) == 0
        warning(['filedate: ' filename ' does not exist.']);
        filedate = -Inf;
    else
        filedate = file_info.datenum;
    end
end
%> @endcond

