# BioMAC-Sim-Toolbox: Coding Guidelines

This coding guidelines explain applied naming conventions and implementation principles. At the moment not all guidelines were applied in the code. Feel free to change and improve the code ;)

## Folders and Files
The folder structure of the protect is currently the following:

````
|-- project home directory
|   |-- docs
|   |   |-- Doxygen
|   |   |-- files
|   |-- src
|   |   |-- HelperFunctions
|   |   |-- model
|   |   |   |-- gait2dc
|   |   |   |   |-- ...
|   |   |   |-- gait2d_osim
|   |   |   |   |-- ...
|   |   |   |-- gait3d
|   |   |   |   |-- ...
|   |   |   |-- Model.m
|   |   |-- problem
|   |   |   |-- @Collocation
|   |   |   |-- Problem.m
|   |   |   |-- ...
|   |   |-- result
|   |   |   |-- ...
|   |   |-- scripts
|   |   |   |-- ...
|   |   |-- solver
|   |   |   |-- IPOPT.m
|   |   |   |-- Solver.m
|   |   |-- tests
|   |   |   |-- Gait2dcTest.m
|   |   |   |-- gait3dTest.m
|   |   |   |-- ModelTest.m
|   |   |   |-- ...
|   |   |-- Toolboxes
|   |   |-- trackingData
|   |   |   |-- @TrackingData
|   |   |-- ...
|	|-- ...
````


In Matlab, classes can have a own class folder named `@<ClassName>`. This folder contains the file defining the model class (e.g.: Gait2dc.m) and other functions of this class (e.g.: simuAccGyro.m). This other functions have to be declared in the class file but the content is implemented in a separate file.


## Classes
- Starting with uppercase
- E.g.: `Model`, `Problem`, `Collocation`, ...

## Functions
- Starting with lowercase
- New word starts with uppercase
- `initModel()`, `initMex()`, `solve()`, ...
- Exception: Constructors

## Properties
- Starting with lowercase
- Numbers:
	- `n<Name>`
	- E.g.: Number of `dofs` is `nDofs`
- Single index:
	- Used in loops
	- `i<Name>`
	- E.g.: 
		```matlab
		for iDof = 1 : obj.nDofs
			disp(num2str(iDof));
		end
		```
- Multiple indices in an array:
	- Parameter
	- `idx<Name>`
	- E.g.: `idxArmdof`
- Hidden properties:
	- Used to save the content of dependent variables
	- `h<Name>`
	- E.g.: `hidxArmdof` for `idxArmdof`
- Constants:
	- Uppercase
	- E.g.: `GRAVITY` in the Gait3d-Class

## Default Values
- Should be set in property definition:
	```matlab
	properties (SetObservable, AbortSet, SetAccess = protected)
		%> Bodyheight in m (default: 1.80)
		bodyheight = 1.80
	end
	```
- Or it should be set in the constructor:
	```matlab
	function [obj] = IPOPT()
		% Set default options for termination
		obj.solverOptions.tol = 0.0002; % Desired convergence tolerance
		obj.solverOptions.max_iter = 4000; % Maximum number of iterations 
		% ...
	end
	```
- Do not change them in the code!
- Use setter function to change the value after object creation. If needed, write an own setter function.

## Inheritance
- Use inheritance instead of changing the functionality of the basic class.
- E.g.: 
	````
	handle
	|--Model
		|-- Gait2dc
		|-- Gait3d
		|-- Gait2d_osim
	````

# How to Use Git Within this Project
First of all use git with all it's advantages. You can find a good introduction (with links to good tutorials) and a good tutorial for advanced users here:
[Git-Tutorial](https://www.atlassian.com/git/tutorials)


## Access
Our project can be found in the MaD GitHub which is hosted on our own server: [Git repository](https://github.com/mad-lab-fau/BioMAC-Sim-Toolbox/tree/main)

## Commits

Please, commit each atomic change and do not combine several changes in one commit. An atomic change is for example a change in a calculation, or a bug fix. 

The commit message should contain this important points:

- One summary line (corresponds to a header in the webview)
- What you changed (not code-wise, but content-wise)
- Why you changed this
- If there were some changes in usage etc.

## Issues and Milestones
Use issues to discuss about a specific topic. It offers a perfect possibility to keep track on the discussion or to come back at a later point in time.

Issues can be assigned to milestones. A milestone could be a publication or the implementation of a new feature.

## Line Endings
Git handles line endings differently on each operation system. To avoid differences produced by line endings please set up your git:

- On Windows: `git config --global core.autocrlf true`
- On Linux/Mac: `git config --global core.autocrlf input`

# How to Document with Doxygen
The new Gait3d-project is documented using [Doxygen](http://www.stack.nl/~dimitri/doxygen/).

It will use all files in the folder src for documentation.

## Doxygen and Matlab
The perl script m2cpp.pl is used to extract comments from Matlab .m files. (http://www.mathworks.com/matlabcentral/fileexchange/25925-using-doxygen-with-matlab)

## Installation on Windows
- Install Doxygen, e.g. from [Sourceforge](http://sourceforge.net/projects/doxygen/)
- Run cmd as administrator and type:
	```bash
	assoc .pl=PerlScript
	ftype PerlScript=C:\Program Files\MATLAB\R2014b\sys\perl\win32\bin\perl.exe %1 %*
	```
	Note: You need to replace the perl path accordingly, which can be found by e.g. typing `which perl` in the command window (MinGW/Cygwin). Perl is installed with e.g. Matlab. This should allow to run a perl script by typing `*.pl` instead of `perl *.pl`. 

	If this does not work, you have to make changes in the windows registry. Open "regedit" by searching it in the windows start menu. 

	Open HKEY_CLASSES_ROOT $\rightarrow$ Applications $\rightarrow$ perl.exe $\rightarrow$ shell $\rightarrow$ open $\rightarrow$ command  and add 
	`\%*` 
	to the value. 
	E.g. change the value from 
	```bash
	"C:\Program Files\MATLAB\R2014b\sys\perl\win32\bin\perl.exe" "%1" 
	```
	to 
	```bash
	"C:\Program Files\MATLAB\R2014b\sys\perl\win32\bin\perl.exe" "%1" %*
	```

	The location in the windows registry could be different.
- Delete the first line of m2cpp.pl file, i.e. delete 
	```bash
	#!/usr/bin/perl.exe
	```

## Installation on UNIX (including Mac OS X)
- Install Doxygen, e.g. from [Sourceforge](http://sourceforge.net/projects/doxygen/)
- Set variable PERL\_PATH in Doxyfile, which can be found by typing `which perl` in the command window, e.g. `PERL_PATH = /usr/bin/perl`
- Change the first line of m2cpp.pl file to your perl path and delete `.exe`:  
	```bash
	#!/usr/bin/perl
	```

## Usage
- If you are using Doxygen on a Linux machine you have to adapt the file docs/Doxygen/Doxyfile. There you have to replace 
	```bash
	FILTER_PATTERNS = *m=m2cpp.pl
	```
	with
	```bash
	FILTER_PATTERNS = *.m="perl m2cpp.pl"
	```
	Please do not commit this changes.
- Open the command-line interface, change directory to the folder docs/Doxygen and run:
	```bash
	doxygen Doxyfile
	```
- Open html documentation: docs $\rightarrow$ Doxygen $\rightarrow$ Doxygen\_HTML $\rightarrow$ index.html

# How to Document

In Matlab files, Doxygen will evaluate only lines starting with `%>`.
Each command starts with `@`.
If you type file, class or function names they will be linked to the respective description.

## Each File
Each file (classes, functions, etc.) should contain first a header describing really shortly the content, the authors and the date:

```matlab
%====================================================================
%> @file @Gait2dc/Gait2dc.m
%> @brief Matlab class describing the Gait2dc model
%>
%> @author Ton, Eva, Marlies, Anne
%> @date October, 2017
%=====================================================================
```

This brief description will be shown in the tab "Files".
A click on it will show the full information.

## Classes
Before each class definition: 
```matlab
%======================================================================
%> @brief The class describes the Gait2dc model
%> @details
%> - The lower body is modeled in the sagittal plane.
%> - The model is configured using an excel file. This code relies on the
%>   structure of this file (hardcoded indices for each section in the file).
%>   Please don't change the structure. Thus the order of the table entries
%>   should not be changed!!!
%> - For the list with the CPs, additional CPs can be added.
%=====================================================================

classdef Gait2dc < Model
	%...
end
```
\medskip
This should contain at least a
- short description marked with `@brief`and starting with an uppercase and
- detailed description marked with `@details` and starting in a new line with an uppercase.

## Properties
The description of each property must be above the property declaration:
```matlab
properties (SetObservable, AbortSet)
	%> Array with lambda values (default: [0.01 0.1]) (1x2)
	lambda = [0.01 0.1]
	%> Gravity in y in m/(s^2)
	gravity
end
````
This description should start with uppercase and should contain
- data type (e.g.: struct, string, array, etc.),
- a description,
- used default value,
- unit,
- size specified in brackets (use variables to define it!).

## Functions
Before each function definition: 
```matlab
%======================================================================
%> @brief Function returning muscle forces for the system in state x
%>
%> @details
%> This function calls the mex file of gait2dc.c:
%> [forces] = gait2dc('Muscleforces', x);
%>
%> @param   obj         Gait2dc class object
%> @param   x           State of the model (obj.nStates x 1)
%>
%> @retval  forces      Muscle forces (in N) (obj.nMus x 1)
%======================================================================
function [muscleForces] = getMuscleforces(obj, x)
	...
end
```
\medskip
This description should contain
- a brief description marked with `@brief` 
	- one line
	- starting with uppercase in the same line
- a detailed description marked with `@details` 
	- can be multiple lines and paragraphs
	- starting with uppercase in the next line
- all function and return parameters marked with `@param` and `@retval`.
	- The formatting is not preserved but recommended for readability in the matlab file.
	- The description should start with an uppercase and should contain
		- data type (e.g.: struct, string, array, etc.),
		- a description,
		- used default value,
		- unit, 
		- size specified in brackets (use variables to define it!).
		- Mark optional parameters with `(optional)` at the beginning of the description.

## Functions in Separate Files
If a function is save in a separate file, there should be no Doxygen comment above the declaration in the class file. Doxygen would link this to the wrong function.

```matlab
% Declared function to simulate acceleration and gyroscope signals
[s, ds_dq, ds_dqd, ds_dqdd] = simuAccGyro(obj, data, q, qd, qdd)
```

In the file containing the implementation, a file header as well as a function header must be used:
```matlab
%======================================================================
%> @file @Gait2dc/simuAccGyro.m
%> @brief Gait2dc function to simulate acceleration and gyroscope signals
%> @details
%> Details: Gait2dc::simuAccGyro()
%> 
%> @author Ton, Eva, Iris, Ann-Kristin, Marlies
%> @date July, 2017
%======================================================================

%======================================================================
%> @brief Function to simulate acceleration and gyroscope signals
%> @details
%> Do not change the "obj" or "data" struct in this function!  This will allow
%> Matlab to "pass by reference" and avoid function call overhead.
%> 
%> @param obj		Gait2dc object
%> @param data		Struct containing the accelerometer and gyroscope data
%> @param q			Generalized coordinates (obj.nDofs x 1)
%> @param qd 		First derivatives of generalized coordinates (obj.nDofs x 1)
%> @param qdd		Second derivatives of generalized coordinates (obj.nDofs x 1)
%> 
%> @retval s 		Simulated sensor signals (data.nVars.all x 1)
%> @retval ds_dq 	Sparse Jacobian ds/dq (data.nVars.all x obj.nDofs)
%> @retval ds_dqd 	Sparse Jacobian ds/dqd (data.nVars.all x obj.nDofs)
%> @retval ds_dqdd 	Sparse Jacobian ds/dqdd (data.nVars.all x obj.nDofs)
%======================================================================
function [s, ds_dq, ds_dqd, ds_dqdd] = simuAccGyro(obj, data, q, qd, qdd)
	...
end
```

The brief file description should state the class the method belongs to. Furthermore, the class should be directly linked in the detailed file description. In this way, you have in the file description the direct link to the function description. 

For 

For a correct assignment, it is recommended to add `@public`, `@protected`, or `@private` after the brief description of the function. This was not done in the example shown above. Consequently Doxygen marked simuAccGyro() as protected even if it was implemented as public function.

## Abstract Functions
Doxygen recognized abstract functions only if you put it in a specific method-block:

```matlab
methods (Abstract)
	% ======================================================================
	%> @brief Abstract function to show model as stick figure
	% ======================================================================
	showStick(obj,x)
end
```

## TODOs
You can found the todo list in the tab "Related Pages". There is a todo list for everybody of us and a general todo list. 

This are the commands to use the lists:
- General todo list: `@todo`
- Named todo list: `@todo<Name>`


It is important to write the todos into the header of a class or a function. Otherwise it will not be linked to the correct content!

Example:
```matlab
%======================================================================
%> @brief Function defining the state table
%>
%> @details
%> Table lists type, name, xmin, xmax and xneutral
%>
%> @todo 
%> Variables of contact points of neutral position shouldn't be zero.
%>
%> @param   obj     Gait2dc class object
%======================================================================
function update_states(obj)
	...
end
```

Even if there is a functionality in Doxygen to report bugs, bugs should be listed in the issues of the git repository.

## Formatting
- White space is not preserved.
- A new line in Matlab does not mean a new line in the documentation.
	- Use `/n` to enforce a new line.
	- Use a blank line to start a new paragraph.
		```matlab
		%> @details 
		%> First paragraph
		%>
		%> Second paragraph
		```
- Code snippets:
	- To show usage of the function.
	- E.g.:
		```matlab
		%======================================================================
		%> @brief Test class to test the class Gait2dc
		%>
		%> @details
		%> - run single test: 
		%> @code 
		%> testCase = Gait2dcTest;
		%> res = run(testCase, 'derivativetest_Dynamics');
		%> @endcode
		%> - run all tests:
		%> @code 
		%>   Tests = matlab.unittest.TestSuite.fromClass(?Gait2dcTest);
		%>   res = run(Tests);
		%> @endcode
		%======================================================================
		classdef Gait2dcTest < ModelTest
			...
		end
		```
- Lists
	- Useful for structs.
	- https://www.stack.nl/~dimitri/doxygen/manual/lists.html:
		- "By putting a number of column-aligned minus (-) signs at the start of a line, a bullet list will automatically be generated. Instead of the minus sign also plus (+) or asterisk (*) can be used." 
		- "Numbered lists can also be generated by using a minus followed by a hash or by using a number followed by a dot."
		- E.g.:
			```matlab
			%> Struct containing initialization parameters with the fields:
			%> - flag:  Boolean showing if the model is initialized. 
			%> - Nf:    Number of functions returned by the Dynamics function.
			%> - fmin:  Lower bounds for dynamics residuals f.
			%> - fmax:  Upper bounds for dynamics residuals f.
			```
- HTML Commands
	- For example for Greek letters.
	- E.g.: `&phi;` 	
