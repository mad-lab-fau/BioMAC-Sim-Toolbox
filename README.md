# BioMAC-Sim-Toolbox: The Simulation Toolbox for Biomechanical Movement Analysis and Creation

## Installation

### BioMAC-Sim-Toolbox

1. Clone this repository:
```bash
git clone https://github.com/MAD-Lab-FAU/BioMAC-Sim-Toolbox.git
```
3. Then, open MATLAB, navigate to the folder where the BioMAC-Sim-Toolbox is located and add all folders and subfolders from "src" to the path: 
```matlab
cd(what('BioMAC-Sim-Toolbox').path); % If not already in the BioMAC-Sim-Toolbox folder
addpath(genpath('src'));
savepath;
```
Alternatively, you can manually add the folder:

* In MATLAB's "Current Folder" browser, navigate to the src folder.
* Right-click -> Add to Path -> Selected Folders and Subfolders.

When calling a model for the first time, it will be compiled automatically. If you would like to compile on different systems on the same machine, you need to manually delete the .o files after building the model files to ensure that you can build on the other system.

### Dependencies

- **MATLAB Compiler**: Required. On Windows, install MinGW via [these instructions](https://nl.mathworks.com/matlabcentral/answers/311290-faq-how-do-i-install-the-mingw-compiler).
- **IPOPT**: Use [mexIPOPT](https://github.com/ebertolazzi/mexIPOPT). Mac users with Apple Silicon should follow [these instructions](https://github.com/ebertolazzi/mexIPOPT/issues/24). For Linux users using gcc versions > 9, use [our mexIPOPT fork](https://github.com/gambimar/mexIPOPT/)
- **OpenSim (Optional)**: Needed for `gait3d` and `gait2d_osim` models. Supported versions: 3.3, 4.0–4.5. For macOS, OpenSim is currently not supported when using native Apple Silicon MATLAB versions. 

You can use the model gait3d, or gait2d_osim, without an OpenSim installation: We have included .mat files that allow you to use the base models and those used in the ExampleScripts without an OpenSim installation. In that case, you should not recalculate the moment arms, which requires OpenSim.

## Using BioMAC-Sim-Toolbox

This is an instruction to start using the BioMAC-Sim-Toolbox. Extended documentation can be found here: https://mad-lab-fau.github.io/BioMAC-Sim-Toolbox/

### Model Loading
The BioMAC-Sim-Toolbox contains three standard models: OpenSim-based `Gait2d_osim` and `gait3d`, as well as `gait2dc`, which is defined by a spreadsheet. A model can be instantiated using its definition file:

```matlab
model_3d = Gait3d('gait3d_pelvis213.osim'); %or
model_2d_osim = Gait2d_osim('gait2d.osim'); %or
model_2dc = Gait2dc('gait2dc_par.xls');
```

### Trajectory Optimization

Set up the optimization problem:

```matlab
% In this example, a `Collocation` object is created for a simulation of 40 Nodes duration and backward Euler integration:
nNodes = 40; Euler = 'BE'; logfile = 'example.log'; plotLog = 1;
problem = Collocation(model_2dc, nNodes, Euler, logfile, plotLog);
```

Define variables: 

```matlab
% States, use default bounds
xmin = repmat(model_2dc.states.xmin, 1, nNodes+1);
xmax = repmat(model_2dc.states.xmax, 1, nNodes+1);
% Fix pelvis translation at initial node 
xmax(model_2dc.extractState('q', 'pelvis_tx'), 1) = 0;
xmin(model_2dc.extractState('q', 'pelvis_tx'), 1) = 0;
% Add states and controls to optimization variables
problem.addOptimVar('states', xmin, xmax);
problem.addOptimVar('controls', repmat(model_2dc.controls.xmin, 1, nNodes+1), ...
                              repmat(model_2dc.controls.xmax, 1, nNodes+1));
% Set duration and speed ranges, in this case fixed target speed
problem.addOptimVar('dur',0.2, 2);
problem.addOptimVar('speed',3.5,3.5);
```

Make initial guess:
```matlab
% Here, previous results could also be supplied
problem.makeinitialguess('mid'); 
```

Add objectives:
```matlab
% Effort term with cubic exponent and equal muscle weighting
problem.addObjective(@effortTermMuscles, 1, 'equal', 3);
```

Add constraints:
```matlab
problem.addConstraint(@dynamicConstraints, ...
                      repmat(model_2dc.constraints.fmin, 1, nNodes), ...
                      repmat(model_2dc.constraints.fmax, 1, nNodes));
problem.addConstraint(@periodicityConstraint, ...
                      zeros(model_2dc.nStates + model_2dc.nControls, 1), ...
                      zeros(model_2dc.nStates + model_2dc.nControls, 1), 1);
```

Solve:
```matlab
solver = IPOPT();
result = solver.solve(problem);
```

Visualize:
```matlab
result.problem.writeMovie(result.X);
```

Note that this example using only muscle-effort minimization yields unrealistic gait, as more objectives are often necessary.

The above example is a simplified version of `script2D`, and we created introduction examples for all models in the toolbox:
- **gait2dc**: `IntroductionExamples.script2D` 
- **gait2d_osim**: `Treadmill.script2D`
- **gait3d**: `IntroductionExamples.script3D`

To run any of the introduction scripts, make sure that the folder "ExampleScripts" is also on your MATLAB path. The folders that are named starting with a + will then also be part of the path.


## Citation

If you use BioMAC-Sim-Toolbox in your research, please cite us. A publication is under review and will be linked here soon.

Follow [Anne Koelewijn](https://www.linkedin.com/in/anne-koelewijn/) on LinkedIn or BlueSky for updates.

```bibtex
@misc{BioMACSimToolbox2025,
  author = {Koelewijn, Anne and others},
  title = {BioMAC-Sim-Toolbox: A Simulation Toolbox for Biomechanical Movement Analysis and Creation},
  year = {2025},
  note = {Manuscript in review},
  url = {https://github.com/MAD-Lab-FAU/BioMAC-Sim-Toolbox}
}
```
