General Usage
-------------
The documentation can be built using doxygen (http://www.stack.nl/~dimitri/doxygen/). 

The perl script m2cpp.pl is used to extract comments from Matlab .m files (http://www.mathworks.com/matlabcentral/fileexchange/25925-using-doxygen-with-matlab).




Installation on Windows
-----------------------
1. Install doxygen, e.g. from http://sourceforge.net/projects/doxygen/

2. Run cmd as administrator and type
	> assoc .pl=PerlScript
	> ftype PerlScript=C:\Program Files\MATLAB\R2014b\sys\perl\win32\bin\perl.exe %1 %*

Note: You need to replace the perl path accordingly, which can be found by e.g. typing "which perl" in the command window (MinGW/Cygwin). Perl is installed with e.g. Matlab. This should allow to run a perl script by typing "*.pl" instead of "perl *.pl." If 2. does not work, you have to make changes in the windows registry. Open "regedit" by searching it in the windows start menu. Open "HKEY_CLASSES_ROOT>Applications>perl.exe>shell>open>command" and add %* to the value, e.g. change the value from "C:\Program Files\MATLAB\R2014b\sys\perl\win32\bin\perl.exe" "%1" to "C:\Program Files\MATLAB\R2014b\sys\perl\win32\bin\perl.exe" "%1" %* . The location in the windows registry could be different.
 
3. Delete the first line of m2cpp.pl file, i.e. delete: #!/usr/bin/perl.exe

4. Open command-line interface, change directory to the documentation folder (docs/Doxygen) and run 
		> doxygen Doxyfile

5.  Open html documentation: documentation>html>index.html



Installation on MAC OS
-------------------------------
1. Install doxygen, e.g. from git repository http://www.stack.nl/~dimitri/doxygen/download.html

2.Set variable PERL_PATH in Doxyfile, which can be found by typing "which perl" in the command window, e.g. PERL_PATH = /usr/bin/perl

3. Change the first line of m2cpp.pl file to your perl path, e.g.:  #!/usr/bin/perl

4. In Doxyfile, add path of m2cpp.pl file to FILTER_PATTERNS = *.m=path_to_file/m2cpp.pl

5. Open command-line interface, change directory to the documentation folder (docs/Doxygen) and run 
		> doxygen Doxyfile

6. Open html documentation: Doxygen_HTML/index.html



Installation on UNIX
-------------------------------
1. Install doxygen, e.g. from git repository http://www.stack.nl/~dimitri/doxygen/download.html

2.Set variable PERL_PATH in Doxyfile, which can be found by typing "which perl" in the command window, e.g. PERL_PATH = /usr/bin/perl

3. Change the first line of m2cpp.pl file to your perl path, e.g.:  #!/usr/bin/perl

4. Replace FILTER_PATTERNS = *m=m2cpp.pl in Doxyfile by FILTER_PATTERNS = *.m="perl m2cpp.pl"

5. Open command-line interface, change directory to the documentation folder (docs/Doxygen) and run 
		> doxygen Doxyfile

6. Open html documentation: Doxygen_HTML/index.html




Documentation Guidelines
-----------------------
In this section, the guidlines for commenting Matlab .m files is described.

The following header should be inserted at the beginning of each .m file:
%======================================================================
%> @file <filename> 
%> @brief A brief description of the file.
%>
%> @author X Y
%> @date Month day, year
%>
%> @bug Bugs can be listed here.
%>
%> @todo Tasks which should be done are listed here.
%>
%> @details
%> A more detailed description of the file follows here.
%======================================================================
Note: The parts @bug, @todo, @details are optionally. You can also start a new paragraph (blank line) then the @details command is not necessary.



The following header should be inserted before each function:
%======================================================================
%> @brief A brief description of the function
%>
%> @details
%> More detailed description here. 
%>
%> @param input_param1 Description of input parameter 1
%> @param input_param2 Description of input parameter 2
%>
%> @retval return_variable1 Description of return variable 1
%> @retval return_variable2 Description of return variable 2
%======================================================================
Note: The part @details is optionally. You can also start a new paragraph (blank line) then the @details command is not necessary.




