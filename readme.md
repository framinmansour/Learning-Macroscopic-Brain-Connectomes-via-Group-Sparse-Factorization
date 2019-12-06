# Learning-Macroscopic-Brain-Connectomes-via-Group-Sparse-Factorization

![alt tag](https://github.com/framinmansour/Learning-Macroscopic-Brain-Connectomes-via-Group-Sparse-Factorization/blob/master/pictures/Selection_355.png)
![alt tag](https://github.com/framinmansour/Learning-Macroscopic-Brain-Connectomes-via-Group-Sparse-Factorization/blob/master/pictures/Selection_356.png)

# About
In this work, we explore a framework that facilitates applying learning algorithms to automatically extract brain connectomes. Using a tensor encoding, we design an objective with a group-regularizer that prefers biologically plausible fascicle structure. 

The resuts of this software can have a long list of potential applications:
* Normal brain development and aging;
* Congenital anomalies, leukodystrophies;
* Tumors and preoperative planning;
* Ischemia and stroke;
* Encephalopathies (toxic, metabolic, infectious);
* Traumatic brain injury;
* Psychiatric disorders, dementia, depression;
* Functional connectivity mapping, cognitive neuroscience;

We provide demos to explain how to:
 (1) Stage 1: initialize the three dimnesional tensor of a brain connectome by assigning an orientation candidate set to each voxel using either of the following greedy forward selection strategies: Orthogonal Matching Persuit (OMP), or our proposed algorithm called GreedyOrientation.
 (2) Stage 2: Build and optimize the objective function consists of the proposed group reguarizer to enforce biolagical plausiblity of the fascicles.
 (3) Visualization: 
 (4) Evaluation: 
 

## [License](TBA).
#### Copyright (2019), [Martha White](https://webdocs.cs.ualberta.ca/~whitem/), whitem@ualberta.ca, [Franco Pestilli](http://francopestilli.com/), frakkopesto@gmail.com, [Cesar Caiafa](http://web.fi.uba.ar/~ccaiafa), ccaiafa@gmail.com, [Russel Greiner](https://rgreiner6.wixsite.com/greiner), [Andrew Patterson](https://www.linkedin.com/in/andy-patterson-1940b068), ap3@ualberta.ca, [Lei Le](http://homes.sice.indiana.edu/leile/), leile@iu.edu, rgreiner@ualberta.ca, [Farzane Aminmansour](https://www.linkedin.com/in/faminmansour/), aminmans@ualberta.ca
 
## [Documentation](TBA).

## [Stable code release](TBA).

## [How to cite the software](TBA).

## Funding.
This research was funded by NSERC, Amii and CIFAR. Computing was generously provided by Compute Canada and Cybera. 
F.P. was supported by NSF IIS-1636893, NSF BCS-1734853, NSF AOC 1916518, NIH NCATS UL1TR002529, a Microsoft Research Award, Google Cloud Platform, and the Indiana University Areas of Emergent Research initiative “Learning: Brains, Machines, Children.

## Installation.
1. Download (https://github.com/framinmansour/Learning-Macroscopic-Brain-Connectomes-via-Group-Sparse-Factorization).
2. [Start MatLab](http://www.mathworks.com/help/matlab/startup-and-shutdown.html).
3. Add repository to the [matlab search path](http://www.mathworks.com/help/matlab/ref/addpath.html).
4. Use **mex** function of __*matlab*__ to compile __*forLoop*__.c in **mex** folder. For more information, please check this [link](https://www.mathworks.com/help/matlab/ref/mex.html). Remember to move the compiled file from **mex** folder to **codes** folder.
5. If not exist yet, please create three folders named as **subsets**, **newsubsets** and **stage1** under **data** folder

## Dependencies.
* [MatLab](http://www.mathworks.com/products/matlab/).
* [C/C++ compiler]()

## Getting started.

### 1. [Download the repository](https://github.com/brain-life/encode).
* Download the repository from the TAR/ZIP files linked [here](https://github.com/framinmansour/Learning-Macroscopic-Brain-Connectomes-via-Group-Sparse-Factorization/archive/master.zip).
* UNZIP/UNTAR the file. 

### 2. [Run the demo_learning_brain_connectome code](/code/demo_learning_brain_connectome.m).
Here you will learn about extracting prerequisit matrices into the folders **subsets** and **newsubsets** in order to run the main script to initialize and optimize brain tensors. Then the demo will visualize fascicles in the resulting brain tensors.
We need to pass two arguments to this function: the first argument determines which dataset to be processed where `1` means Arcuate and `2` means ARC-SLF. The second argument determines how do we want the predicted brain tensor Phi to be initialized, `1` initializes the predicted_Phi with the expert_Phi, `2` initalizes it with the output of OMP algorithm, and `3` will use GreedyOrientation method.
```
  >>  demo_learning_brain_connectome(#dataIndex, #stage1)
```
__*brain*__.m,  __*stageOne*__.m and __*stageTwo*__.m are main files to do optimization and deliver results. Among them, __*brain*__.m is the one users will call the other two matlab files and save results; __*stageOne*__.m is to do jobs at stage 1 (initialize Phi with different strategies as warm start in the optimization at stage2) and to call __*stage2*__.m; __*stageTwo*__.m is to do jobs at stage 2 (optimize Phi). To know more about each file, please check their comments

### 3. [Run other scripts to reproduce evaluation results](/code/figures).
The code files in **figures** folder are to draw figures on the results to evaluate OMP and GreedyOrientation in the **Screening stage** or **stage1**.
