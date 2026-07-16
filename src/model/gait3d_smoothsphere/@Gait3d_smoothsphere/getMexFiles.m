function name_MEX = getMexFiles(osimName, rebuild_contact, debugMode)
    if nargin < 2
        rebuild_contact = 0;
    end
    if nargin < 3
        debugMode = 0;
    end

    % 1. Compile standard model dependencies first and get parent MEX name
    parent_name_MEX = Gait3d.getMexFiles(osimName);

    % 2. Set directory to C files folder
    currentFolder = pwd;
    methodFolder = strrep(mfilename('fullpath'), [filesep , mfilename], '');
    cfilesFolder = fullfile(fileparts(fileparts(methodFolder)), 'gait3d', 'c_files');
    cd(cfilesFolder);
    
    % Name of our smoothsphere MEX
    name_MEX = [parent_name_MEX '_smoothsphere'];
    
    % Check Unix/Windows options
    if ispc
        objext = 'obj';
    elseif isunix
        objext = 'o';
        mexoptUNIX = 'GCC=''/usr/bin/gcc''';
    else
        error('System is unknown');
    end
    
    % 3. Run SymPy generator to build/update contact_3d_smoothsphere_al.c
    if ~exist('contact/contact_3d_smoothsphere_al.c', 'file') || rebuild_contact
        fprintf('Generating smoothsphere contact C code using SymPy...\n');
        % run using uv env
        [status, cmdout] = system('uv run --with sympy python3 contact/generate_contact_3d.py');
        if status ~= 0
            error('Failed to generate contact C code: %s', cmdout);
        end
    end
    
    % 4. Compile contact_3d_smoothsphere_al.c to object file
    contact_c_file = 'contact/contact_3d_smoothsphere_al.c';
    contact_obj_file = ['contact/contact_3d_smoothsphere_al.' objext];
    
    fprintf('Compiling custom contact model object...\n');
    if isunix
        mex(mexoptUNIX, '-largeArrayDims', '-c', contact_c_file, '-outdir', 'contact');
    else
        mex('-largeArrayDims', '-c', contact_c_file, '-outdir', 'contact');
    end
    
    % 5. Compile the new MEX binary linking the model and objects
    model_c_file = [parent_name_MEX '/' parent_name_MEX '_smoothsphere.c'];
    
    % Check if smoothsphere.c exists, if not copy and adapt it
    if ~exist(model_c_file, 'file')
        % Copy original c file
        copyfile([parent_name_MEX '/' parent_name_MEX '.c'], model_c_file);
        % Read and replace inclusion of .h with _smoothsphere.h
        content = fileread(model_c_file);
        content = strrep(content, ['#include "' parent_name_MEX '.h"'], ['#include "' parent_name_MEX '_smoothsphere.h"']);
        
        % Insert get_bodyweight function below parameters declaration
        target = 'static param_struct parameters;';
        replacement = sprintf('static param_struct parameters;\n\ndouble get_bodyweight() {\n\treturn parameters.bodyweight;\n}\n');
        content = strrep(content, target, replacement);
        
        fid = fopen(model_c_file, 'w');
        fwrite(fid, content);
        fclose(fid);
    end
    
    obj_al = [parent_name_MEX '/' parent_name_MEX '_al.' objext];
    
    % Check if separate FK and NoDer files exist (like for pelvis213 models)
    obj_noder = [parent_name_MEX '/' parent_name_MEX '_NoDer_al.' objext];
    obj_fk = [parent_name_MEX '/' parent_name_MEX '_FK_al.' objext];
    
    fprintf('Compiling new smoothsphere model MEX...\n');
    if exist(obj_noder, 'file') && exist(obj_fk, 'file')
        if isunix
            mex(mexoptUNIX, '-largeArrayDims', model_c_file, obj_al, obj_noder, obj_fk, contact_obj_file, '-output', [parent_name_MEX '/' name_MEX]);
        else
            mex('-largeArrayDims', model_c_file, obj_al, obj_noder, obj_fk, contact_obj_file, '-output', [parent_name_MEX '/' name_MEX]);
        end
    else
        if isunix
            mex(mexoptUNIX, '-largeArrayDims', model_c_file, obj_al, contact_obj_file, '-output', [parent_name_MEX '/' name_MEX]);
        else
            mex('-largeArrayDims', model_c_file, obj_al, contact_obj_file, '-output', [parent_name_MEX '/' name_MEX]);
        end
    end
    
    cd(currentFolder);
end
