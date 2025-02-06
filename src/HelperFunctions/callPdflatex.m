%======================================================================
%> @file HelperFunctions/callPdflatex.m
%> @brief Function to call pdflatex
%>
%> @author Marlies Nitscke
%> @date January, 2019
%======================================================================

% ======================================================================
%> @brief Function to call pdflatex
%>
%> @details
%> pdflatex has to be installed. 
%>
%> @param  folder    String: Name of the folder containing the .tex file
%> @param  filename  String: Name of the .tex file without the extension ".tex"
%> @retval status    String: Exit status of the command. When the command is successful, 
%>                   status is 0. Otherwise, status is a nonzero integer.
%> @retval cmdout    String: Output of the command
% ======================================================================
function [status, cmdout] = callPdflatex(folder, filename)

    % Change current directory
    currentFolder = pwd;
    cd(folder);
    
    % Create command string
    command = sprintf('pdflatex --jobname=%s %s.tex', filename, filename);
    
    % Execute the command
    [status, cmdout] = system(command, '-echo');
    
    % delete .aux, .log, .out files if callPdflatex ran without error
    if status == 0
        filelist = [dir(cat(2,filename,'.aux')); ...
            dir(cat(2,filename,'.log')); ...
            dir(cat(2,filename,'.out'))];
        x = 1:numel(filelist);
        arrayfun(@(x) delete(filelist(x).name), x);
    end
    
    % Change directory back
    cd(currentFolder);
    
end