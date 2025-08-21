---
title: 'BioMAC-Sim-Toolbox: A MATLAB toolbox for biomechanical motion analysis and creation through simulation'
tags:
  - MATLAB
  - Biomechanics
  - Movement Simulation
  - Movement Analysis
authors:
  - name: Anne D. Koelewijn  
    orcid: 0000-0002-1867-0374
    corresponding: true
    affiliation: "1, 2"
  - name: Marlies Nitschke
    orcid: 0000-0001-7318-9578
    affiliation: 1
  - name: Eva Dorschky
    orcid: 0000-0001-8708-0426
    affiliation: 1
  - name: Markus Gambietz
    orcid: 0000-0003-4096-2354
    affiliation: 1
  - name: Alexander Weiss
    orcid: 0009-0003-1604-6679
    affiliation: 1
  - name: Bjoern M. Eskofier
    orcid: 0000-0002-0417-0336
    affiliation: 1
  - name: Antonie J. van den Bogert
    orcid: 0000-0002-3791-3749
    affiliation: 2
affiliations:
 - name: Friedrich-Alexander-Universität Erlangen-Nürnberg, Erlangen, Germany
   index: 1
 - name: Cleveland State University, Cleveland, Ohio, United States of America
   index: 2
date: 13 August 2024
bibliography: paper.bib
---

# Summary
We present a MATLAB toolbox for human movement simulation and analysis. The toolbox's focus is on creating simulations of human walking and running by solving 
trajectory optimization problems. **In these problems, the muscle excitations  and initial state are found that minimize an objective related to energy expenditure. **
To solve the trajectory optimization problems, we have implemented a direct collocation framework.
Different objectives have been implemented that represent energy or effort minimization. We have also included 
different objectives for data tracking, including tracking of joint angles, ground reaction forces, accelerometer data, gyroscope data, and marker positions. We have implemented
two different two dimensional (2D) **musculoskeletal** models and one three dimensional (3D) musculoskeletal model. 
To show how the toolbox can be used, we include a **tutorial**, a 2D and a 3D introductory example, as well as different 
applications that are based on previous publications. The toolbox is versatile due 
to its object-oriented design, such that different dynamics models can be easily combined with different 
problem classes. In the future, we aim to implement different other 
problem classes as well, such as inverse kinematics and dynamics to process experimental data, static and dynamic optimizations to find muscle simulations from joint moments, while also 
forward shooting can be used to investigate neural control of gait. 

# Statement of Need
Movements of humans and other animals are extremely versatile, while our efficiency is unmatched by human-made 
machines without having accurate or efficient controllers. Therefore, a better understanding of human movement can have a large impact for persons with movement disabilities, sports performance optimization, 
and design of human-made devices. Movement simulations are a key tool for creating this understanding.
On the one hand, we can use so-called predictive simulations [@ackermann:2010; @falisse:2019;@koelewijn:2018] to investigate unseen movements by replicating the optimization in the central nervous system [@zarrugh:1974] in a computer optimization. **Such simulations can help us decipher fundamentals of movement **
**in humans and other animals.**
On the other hand, we can use so-called reconstructive simulations to estimate variables that cannot or were not measured directly by minimizing a data tracking error [@dorscky:2019a;@nitschke:2023;@nitschke:2024]. **These reconstructive simulations do not necessarily require laboratory equipment, **
**instead allowing biomechanical variables to be estimated accurately from wearable or contactless sensors, such as inertial sensors. Wide-spread measurement of movements and their biomechanical analysis are vital to be able to **
**reach a deep understanding** of human movement.

Here, we present the `BioMAC-Sim-Toolbox`, a MATLAB toolbox that can be used to simulate
and analyse human movements. **This toolbox is released under an Apache-2.0 license and can be accessed on GitHub: <https://github.com/mad-lab-fau/BioMAC-Sim-Toolbox>.** The main functionality of the toolbox is that it solves trajectory optimization problems, or optimal control problems, for human musculoskeletal dynamics models. These dynamics models combine multibody 
dynamics to model the movements of the skeleton with muscle dynamics models to model the dynamics of muscular 
contraction and activation. In the trajectory optimization problem, the muscle **excitations and initial state, as well as possibly other gait variables such as speed or duration,** are found that minimize an objective. This objective generally includes terms to minimize muscular effort and to minimize a tracking error
with respect to measured movement data [@koelewijn:2016;@dorscky:2019a;@dorschky:2019b;@koelewijn:2022;@nitschke:2023;@nitschke:2024], but different other objectives have been implemented as well, such as 
minimization of metabolic energy expenditure [@koelewijn:2018]. The optimization problem can be described mathematically as follows:

$$ \underset{\mathbf{x}(0),\mathbf{u}(t), T} {\text{minimize}} \quad J(\mathbf{x}, \mathbf{u}) =  \sum_{j=1}^{N_W} W_{j} \int_{t=0}^T c_{j}(x(t),u(t)) \mathrm{d} t + c_f(x(T),u(T)) $$

$$ \text{subject to} \quad  \mathbf{f}(\mathbf{x}(t), \mathbf{\dot{x}}(t), \mathbf{u}(t)) = \mathbf{0} \ \text{for }0 \leq t \leq T \ \text{(system dynamics)} $$

$$ \mathbf{x_L} \leq \mathbf{x}(t) \leq \mathbf{x_U} \ \text{for }0 \leq t \leq T \quad \text{(bounds on states)} $$

$$ \mathbf{u_L} \leq \mathbf{u}(t) \leq \mathbf{u_U} \ \text{for }0 \leq t \leq T \quad \text{(bounds on controls)} $$

for a musculoskeletal model with state $x$ and input $u$. Here, the goal is to find the initial state, $\mathbf{x}(0)$, control inputs, $\mathbf{u}(t)$,
and duration, $T$, that minimize the cost function, $J(\mathbf{x}, \mathbf{u})$. The cost function consists of $c_{j}(\mathbf{x}(t),\mathbf{u}(t))$, which 
describes the cost functions that are evaluated over time, which are multiplied with associated **weight** $W_{j}$, and of a final cost function $c_f(\mathbf{x}(T),\mathbf{u}(T))$, which is only evaluated at the end of the trajectory, at time $T$. We normally do not include a final cost function. The subscript $L$ defines the lower bounds, while $U$ defines the upper bounds. When speed or duration is known, we commonly still include them into the optimization variables, while making the lower and upper bound equal to the desired value. 

For walking and running, we include a periodicity constraint to the full gait cycle or a single step, in case of a symmetric motion. This constraint can be defined for straight and curved running and walking [@nitschke:2020] through the displacement in the horizontal plane. This periodicity constraint ensures that the states are the same at the end of the trajectory as at the beginning, while applying a horizontal displacement to the states in $\mathbf{x}_{hor}$:

$$ \mathbf{x}(T) = \mathbf{x}(0) + v T \mathbf{x}_{hor} \quad \text{(periodicity)} $$

So far, we have created the problem class `Collocation`, which creates a trajectory optimization problem using 
direct collocation [@betts:2010]. In direct collocation, the numerical integration of the differential equations that represent 
the dynamics are solved as equality constraints during the optimization. To do so, collocation nodes, or time nodes, 
need to be predefined, and the integration method should be selected. We have currently implemented a backward Euler and
a midpoint Euler integration method. Backward Euler has been more stable and has therefore been used for publications. A 
direct collocation approach creates a large-scale nonlinear optimization problem, for which the derivatives of the 
objectives and all constraint can be derived analytically. An optimization problem including a periodicity constraint could be transcribed as follows after applying direct collocation with backward Euler integration:

$$ \underset{\mathbf{x}(i),\mathbf{u}(i), v, T} {\text{minimize}} \quad J(\mathbf{x}, \mathbf{u}) =  \sum_{j=1}^{N_W} W_{j} \sum_{i=0}^N c_{j}(\mathbf{x}(i),\mathbf{u}(i)) + c_f(\mathbf{x}(N),\mathbf{u}(N)) $$

$$ \text{subject to} \quad  \mathbf{f}(\mathbf{x}(i+1), \frac{\mathbf{{x}}(i+1)-\mathbf{x}(i)}{\Delta t}, \mathbf{u}(i)) = \mathbf{0} \ \text{for }i=1,\dots,N \quad \text{(system dynamics)} $$

$$ \mathbf{x_L} \leq \mathbf{x}(i) \leq \mathbf{x_U} \ \text{for }i=1,\dots,N \quad \text{(bounds on states)} $$

$$ \mathbf{u_L} \leq \mathbf{u}(i) \leq \mathbf{u_U} \ \text{for }i=1,\dots,N \quad \text{(bounds on controls)} $$

$$ \mathbf{x}(N+1) = \mathbf{x}(0) + v T \mathbf{x}_{hor} \quad \text{(periodicity)} $$

for $N$ collocation points $i$. **Compared to the description of the optimization problem in continuous time, there are only subtile changes. **
**The integral in the objective has changed into a sum, and we have used the backward Euler transcription in the dynamics. Furthermore, the bounds **
**and the periodicity constraint are now applied to discrete variables.**

Problems can be solved using the solvers that are implemented in the class `solver`. So far, we have implemented **an interface to** IPOPT [@Wachter:2006] to solve 
the resulting optimization problem. This algorithm can solve optimal control problems transcribed with direct collocation in less than 10 minutes 
for a 2D model, **for example on a computer with an Intel Core i5-3210M CPU at 2.5 GHz clock speed** [@koelewijn:2016; @koelewijn:2022] and less than 
one hour for a 3D model **on a single core of a workstation with a 3.2GHz Xeon E5-1660v4 processor** [@nitschke:2020]. Different solvers could easily 
be integrated into this solver class.

As most analysis is specific to the problem type, most code for analysis can be found in the respective problem class. For example, metabolic cost requires an integration over time, which is dependent on the problem type that is being used. In addition, general methods, e.g., to calculate a correlation or root mean squared error between two variables, or to write results into an OpenSim file format, can be found in the function `HelperFunctions` folder. 

We have implemented different human dynamics models (see Table 1). All implemented models are musculoskeletal models, but they can 
also be used as skeletal models, such that the input is generated using joint moments instead of muscle **excitations**. 
We have implemented a 3D model version (`gait3d`) and two 2D model versions in c, which are compiled as mex functions. The model dynamics can be 
tested using the test cases coded in the `tests` folder. The **model parameters of the** 3D model and one 2D model (`gait2d_osim`) are loaded from OpenSim, but use our own 
muscle dynamics model, which is described in [@nitschke:2020], while the other (`gait2dc`) is used in our own previous work 
[@koelewijn:2016; @koelewijn:2022; @dorscky:2019a; @dorschky:2019b]. The model parameters are defined in the .osim file for the OpenSim models, 
and defined in an Excel file for model gait2dc (gait2dc_par.xlsx). The models can be personalized by directly adjusting the parameters in the 
.osim model or Excel file, such that, for example, OpenSim scaling can be used [@seth:2018]. They can also still be adjusted in MATLAB, for 
example to investigate virtual participants [@dorscky:2019a; @koelewijn:2022]. The model dynamics are explained further for `gait2dc` in
@koelewijn:2022; @dorscky:2019a; @dorschky:2019b and for `gait3d` in @nitschke:2020. **The `gait2d_osim model` has been used in @gambietz:2024.** 
It is based on the `gait10dof18musc.osim` model [@seth:2018], and can be loaded in two ways: the original version, called `gait10dof18musc`, and a version with 
the lumbar joint locked, called `gait2d_osim`.

| Model Name      | Dimensionality | Degrees of freedom | Number of muscles | Arm torque actuators | Source | 
| --------------- | -------------  | ------------------ | ----------------  | -------------------- | ------ | 
| gait2dc         |      2D        |          9         |         16        |    -                 |[@koelewijn:2016; @koelewijn:2022; @dorscky:2019a; @dorschky:2019b] |
| gait2d_osim     |      2D        |          9         |         18        |    -                 |[@seth:2018]
| gait10dof18musc |      2D        |          10        |         18        |    -                 |[@seth:2018]
| gait3d          |      3D        |          33        |         92        |       10             |[@nitschke:2020; @nitschke:2023; @nitschke:2024] |
**Table 1:** Overview of dynamics models implemented in the BioMAC-Sim-Toolbox

**To get started, we recommend that users first look at the tutorial in the folder `Tutorial` and the folder `IntroductionExamples` in 
**`ExampleScripts`. The tutorial is created to help users become familiar with the inputs and outputs of the simulations, as well as create and run simulations.
**To perform the tutorial, users should start by looking at the file TutorialInTheWild.PDF. In addition, we have added an introductory example for 2D simulations (`script2D.m`) and
**3D simulations (`script3D.m`), which is based on @nitschke:2020, in the folder `IntroductionExamples`. The goal of these introductory examples is to show how the toolbox works, and
**highlight different options in the implementation. Therefore, these functions, and those called inside them, have a lot of comments. We recommend
**running them line by line to follow along with the example. 

**Besides the introductory examples, we have included several other examples in the folder `ExampleScripts`. These examples all run on their own through the function
**in the folder that has a starting with `script`. They should produce a relevant graph or movie at the end. These examples are added to highlight different ways that the code can be used. The code in the folder `ExoPaper` can be used to generate predictive simulations of walking with an exoskeleton, based on 
**the paper by @koelewijn:2022. The example is added to show how models can be extended by assistive devices, in this case an exoskeleton, and to show
**how predictive simulations can be generated with the `BioMAC-Sim-Toolbox`. The code in the folder `IMUTracking2D` can be used to generate tracking
**simulations based on inertial sensor data. This example is added to show how inertial sensor data can be tracked in a 2D musculoskeletal model,
**and is based on the publications by @dorscky:2019a and @dorschkynitschke:2025. The code in `MarkerTracking3D` can also be used to generate tracking
**simulations, but from marker and ground reaction force data. This example is added as an example of how the 3D musculoskeletal model has been used and
**is based on the publication by @nitschke:2024. The code in `Treadmill`  can be used to generate simulations of walking on the treadmill. This example
**is not based on published work, but was added to show two how a treadmill can be implemented in the ground contact model and how the `gait2d_osim` model can be used in the `BioMAC-Sim-Toolbox`.

Currently, the toolbox can be used to solve trajectory optimization problems to solve predictive and reconstructive simulations. We have implemented code to track marker positions, joint angles, translations, ground reaction forces, 
accelerations, angular velocities, movement 
duration, and movement speed. We have further implemented objectives to minimize muscular effort, metabolic cost, squared torque, joint 
accelerations, passive joint moments, ground reaction force impact, and head stability. By combining these different objectives, many different 
types of simulations can be generated and the movement kinematics and kinetics investigated. **To generate new simulations, it is necessary to write 
**scripts in MATLAB. We recommend to use functions, as done in the examples, to set up the different problems. Currently, it is not possible to 
**use the software in a GUI. Due to these limitations, the learning curve is steep for new users. To help new users, we have added slides that were presented 
**in a tutorial at the 2025 Conference of the International Society of Biomechanics to explain the concepts that are used in the toolbox. 

**Several software suits exist with similar functionality, such as OpenSim Moco [@dembia:2020], Bioptim [@michaud:2022], and PredSim [@dhondt:2024]. 
**In OpenSim Moco, more time is required to solve, e.g., a tracking problem required 130 minutes [@dembia:2020], where a similar problem takes 8 minutes with the 
**`BioMAC-Sim-Toolbox` [@koelewijn:2016]. Furthermore, Bioptim and PredSim do not include functionality to track experimental measurement 
**data and therefore are limited to predictive simulations, while the `BioMAC-Sim-Toolbox` can be used for both tracking and predictive simulations. 
**A disadvantage of the `BioMAC-Sim-Toolbox` is that we have, so far, only implemented one-step discretization schemes, which are not as 
**accurate as higher-order schemes, which are implemented in the other toolboxes. We are currently working on an implementation of higher order 
**collocation schemes as well, and plan to add methods up to the fifth order, which are two-step and three-step. 

In future, the toolbox could be expanded. **We invite users to suggest new features by creating issues and also to help implement them through pull requests. We encourage 
**users to add models, objective and constraint functions to the toolbox to expand the range of optimal control problems that can be solved.** Furthermore, new problem 
classes can be defined, for example for inverse kinematics and inverse dynamics, or to solve muscle activations and excitations statically or dynamically. It would also 
be possible to implement single shooting approach, which would allow for models with discontinuous dynamics or control to be investigated, such as a reflex model [@Geyer:2010]. 
**Through these expansions, the `BioMAC-Sim-Toolbox` would also allow users to implement and compare different movement simulation and analysis approaches and ensure that the 
**outcomes as as good as possible for their research questions.

# Acknowledgements

We acknowledge contributions from Marko Ackermann, Dieter Heinrich, Maria Eleni Athanasiadou, Chuyi Wang, Christopher Löffelmann, Linus Hötzel, Utkarsha Shukla, Arne Kuederle, Tobias Luckfiel, Heiko Schlarb. The development of the BioMAC-Sim-Toolbox was supported by adidas AG and by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) through Project-ID442419336, SFB 1483 - EmpkinS and through Project-ID520189992

# References
