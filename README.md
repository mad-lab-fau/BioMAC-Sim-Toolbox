# BioMAC-Sim-Toolbox: The Simulation Toolbox for Biomechanical Movement Analysis and Creation

## Installation

### BioMAC-Sim-Toolbox

First, clone this repository. Then, open MatLab, navigate to the folder where the BioMAC-Sim-Toolbox is located and add all folders and subfolders from "src" to the path: 
```matlab
cd(what('BioMAC-Sim-Toolbox').path); % If not already in the BioMAC-Sim-Toolbox folder
addpath(genpath('src'));
savepath;
```
Alternatively, you can manually add the toolbox to the path by going to the "Current Folder Browser" and navigate to the folder where "src" is located. Right click, go to "add to path" and click "selected folders and subfolders".

When calling a model for the first time, it will be compiled automatically. If you would like to compile on different systems on the same machine, you need to manually delete the .o files after building the model files to ensure that you can build on the other system.

### Dependencies

**Compilers**

First, make sure that MatLab Compiler is installed. For **Windows**, mingw is needed. For instructions, see here: https://nl.mathworks.com/matlabcentral/answers/311290-faq-how-do-i-install-the-mingw-compiler. For Linux/MacOS, gcc should be installed by default.

**IPOPT**

We recommend using the obtaining the following version through the Add-On Explorer in MATLAB: https://github.com/ebertolazzi/mexIPOPT. For Mac Users using Apple Silicon native MatLab versions, follow these instructions: https://github.com/ebertolazzi/mexIPOPT/issues/24

**OpenSim** (optional)

When using the model gait3d, or gait2d_osim, it is recommended to install OpenSim and its MatLab API from https://simtk.org/frs/?group_id=91. We currently support versions 3.3, 4.0, 4.1, 4.3, and 4.5. For Mac Users using Apple Silicon native MatLab versions, OpenSim is not supported currently. 

You can use the model gait3d, or gait2d_osim, without an OpenSim installation: We have included .mat files that allow you to use the base models and those used in the ExampleScripts without an OpenSim installation. In that case, you should not recalculate the moment arms, which requires OpenSim.

## Using BioMAC-Sim-Toolbox

This is an instruction to start using the BioMAC-Sim-Toolbox. Extended documentation can be found here: https://mad-lab-fau.github.io/BioMAC-Sim-Toolbox/

### Models
The BioMAC-Sim-Toolbox contains three standard models: OpenSim-based `Gait2d_osim` and `gait3d`, as well as `gait2dc`, which is defined by a spreadsheet. A model can be instanciated using its definition file:
```matlab
model_3d = Gait3d('gait3d_pelvis213.osim'); %or
model_2d_osim = Gait2d_osim('gait2d.osim'); %or
model_2dc = Gait2dc('gait2dc_par.xls');
```

### Trajectory Optimization
In this example, a `Collocation` object is created for a simulation of 40 Nodes duration and backward Euler integration:

```matlab
nNodes = 40; Euler = 'BE'; logfile = 'example.log'; plotLog = 1;
problem = Collocation(model_2dc, nNodes, Euler, logfile, plotLog);
```
Next, we define which variables we want to optimize: states and controls. For all, we use the default lower and upper bounds that our model provides. 
Our target is a simulation at 3.5m/s, therefore, we fix the speed and set the duration loosely to an interval between 0.2s and 2s.
```matlab
xmin = repmat(model_2dc.states.xmin, 1, nNodes+1);
xmax = repmat(model_2dc.states.xmax, 1, nNodes+1);
xmax(model_2dc.extractState('q', 'pelvis_tx'), 1) = 0;
xmin(model_2dc.extractState('q', 'pelvis_tx'), 1) = 0; 
problem.addOptimVar('states', xmin, xmax);
problem.addOptimVar('controls',repmat(model_2dc.controls.xmin,1,nNodes+1), repmat(model_2dc.controls.xmax,1,nNodes+1));
problem.addOptimVar('dur',0.2, 2);
problem.addOptimVar('speed',3.5,3.5);
```

Then, we supply the initial guess to our problem. For simplicity, we set the initial guess to 'mid', which is the middle between upper and lower bound. Setting a good initial guess can speed up the simulation, see the ExampleScripts.
```matlab
problem.makeinitialguess('mid'); 
```

**Objectives**
In this simple example, we only optimize for muscle effort, therefore we set the weight to 1. An objective can be added the following way:
```matlab
Weffort = 1;
weightsType = 'equal'; 
exponent = 3; 
problem.addObjective(@effortTermMuscles, Weffort, weightsType, exponent);
```

**Constraints**
In our example, we make use of `@dynamicConstraints` for dynamically consistent output. We also set `@periodicityConstraint` for half-gait-cycle periodicity.
```matlab
problem.addConstraint(@dynamicConstraints,repmat(model_2dc.constraints.fmin,1,nNodes),repmat(model_2dc.constraints.fmax,1,nNodes))
problem.addConstraint(@periodicityConstraint,zeros(model_2dc.nStates+model_2dc.nControls,1),zeros(model_2dc.nStates+model_2dc.nControls,1),1)
```

**Calling IPOPT**
By adding objectives and constrained, the problem is now defined and can be solved using IPOPT in our default settings:
```matlab
solver = IPOPT();
result = solver.solve(problem);
```

**Inspecting Results**
To show a video of the resulting motion, call the `Collocation.writeMovie` function. 
```matlab
result.problem.writeMovie(result.X);
```
Note that the output of this example isn't quite realistic, as muscle effort alone is not sufficient to describe walking or running.

The above example is a simplified version of `script2D`, and we created introduction examples for all models in the toolbox:
- **gait2dc**: `IntroductionExamples.script2D` 
- **gait2d_osim**: `Treadmill.script2D`
- **gait3d**: `IntroductionExamples.script3D`

To run any of the introduction scripts, make sure that the folder "ExampleScripts" is also on your MATLAB path. The folders that are named starting with a + will then also be part of the path.


## Citation

If you are using our toolbox in your research, we would be glad for a citation. Our _paper_ is still under review, we will update the citation. Follow us on ... for updates.
