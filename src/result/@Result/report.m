%======================================================================
%> @file @Result/report.m
%> @brief Result function to create graphical report from a result
%> @details
%> Details: Result::report()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Matlab function to create graphical and textual report from a result
%>
%> @details
%> Please note that no PDF is generated, if no resultFilename is given as input argument!
%>
%> @param    obj             Result object
%> @param    settings        (optional) Struct: Settings defining what to report. It combines
%>                           the fields specified in Collocation.extractData() and Collocation.report().
%>                           (use empty to skip)
%> @param    style           (optional) Struct: Settings defining the style for reporting.
%>                           See plotVarType() for details. (use empty to skip)
%> @param    resultFilename  (optional) String: Filename to save the report in a LaTeX 
%>                           document. If the filename is an empty string, the result.filename 
%>                           will be used to save the LaTeX document. If there is no input for 
%>                           the filename at all, the report will be printed in the console
%>                           and nothing will be saved. 
%> @param    saveSeparatly   (optional) Boolean: If true, save the solver report will be saved in a 
%>                           separate .tex and all figures will be saved as .fig and tikz
%>                           standalone. Figures will be saved as .fig and tikz standalone. (default: 0)
%> @retval simVarTable       Table: Summarizing all data which was reported
%>                           (see Collocation.extractData() for details)
%======================================================================
function simVarTable = report(obj, settings, style, resultFilename, saveSeparatly)

% Set default
if nargin < 2
    settings = struct();
end
if nargin < 3
    style = struct();
end

% Check if we have to save the results
if nargin > 3
    saveIt = 1;
    % Use filename in result if string is empty
    if isempty(resultFilename)
        resultFilename = result.filename; % Use this if no resultFilename is given
    end
else
    saveIt = 0;
end

% Check if we have to save the results of solver and problem separatly
if nargin < 5
    saveSeparatly = 0;
end

% Warn if it did not converged
if obj.converged ~= 1
    warning('Result:report', 'This optimization did not converge.');
end


% Get report of the Solver
doReportSettings = 1;
if saveSeparatly
    % Input also filename such that the results will be saved separately
    [conStringSolver, texStringSolver] = obj.solver.report(doReportSettings, obj.info, resultFilename);
else
    [conStringSolver, texStringSolver] = obj.solver.report(doReportSettings, obj.info);
end


% Get report of the Problem
if saveSeparatly
    % Input also filename such that the results will be saved separately
    [conStringProblem, texStringProblem, simVarTable] = obj.problem.report(obj.X, settings, style, resultFilename);
else
    [conStringProblem, texStringProblem, simVarTable] = obj.problem.report(obj.X, settings, style);
end


% Return report
if saveIt
    % Combine the two latex reports into one document
    texString = [texStringSolver, texStringProblem];
    % Save the report 
    docTitle = sprintf('Report of file %s', strrep(strrep(resultFilename,'_','\_'), '^', '$\hat{}$'));
    saveLaTeXDocument(texString, [resultFilename '_ResultInfo'], docTitle);
    % Create the pdf from tex file
    [folderName,fileName]=fileparts(resultFilename);
    fileName=[fileName,'_ResultInfo'];
    % check whether the tex file has no errors so the pdf can be created
    % for sure
    try
        callPdflatex (folderName,fileName);
    catch
        errorMessage = sprintf('Could not make the pdf since the associated tex file has errors. Check the file:\n %s_ResultInfo.tex\n', resultFilename);
        warning (errorMessage);
        
    end
    
else
    % Display the result in the console
    disp(conStringSolver);
    disp(conStringProblem);
end


end
