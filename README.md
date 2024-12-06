Instruction to start using the toolbox

1. Make sure that you have compiler installed. For instructions, see here: https://nl.mathworks.com/matlabcentral/answers/311290-faq-how-do-i-install-the-mingw-compiler
2. Install IPOPT. For example, you can install it by installing the following version through the Add-On Explorer in MATLAB: https://github.com/ebertolazzi/mexIPOPT. To do so, go to the "Apps" tab, click "Get More Apps" in top left to open the Add-On Explorer, then search for "IPOPT" and click "Add" on the top right. 
3. When using the model gait3d, or gait2d_osim, it is recommended to install OpenSim from https://simtk.org/frs/?group_id=91. We currently support versions 3.3, 4.0, 4.1, 4.3, and 4.5. We have included .mat files that allow you to use the base models and those used in the ExampleScripts without an OpenSim installation. In that case, you should not recalculate the moment arms, which requires OpenSim
4. Add all folders and subfolders from "src" to the path. In the "Current Folder Browser", navigate to the folder where "src" is located. Right click, go to "add to path" and click "selected folders and subfolders".
5. To run any of the introduction scripts, make sure that the folder "ExampleScripts" is also on your MATLAB path. The folders that are named starting with a + will then also be part of the path.

If you are using a newer Apple Silicon Mac, please take a look at the following to get mexIPOPT running: https://github.com/ebertolazzi/mexIPOPT/issues/24

If you would like to compile on different systems on the same machine, you need to manually delete the .o files after building the model files to ensure that you can build on the other system.