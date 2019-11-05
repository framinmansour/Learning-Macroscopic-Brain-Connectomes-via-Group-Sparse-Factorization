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
#### Copyright (2019), [Martha White](https://webdocs.cs.ualberta.ca/~whitem/), whitem@ualberta.ca, [Franco Pestilli](http://francopestilli.com/), frakkopesto@gmail.com, [Cesar Caiafa](http://web.fi.uba.ar/~ccaiafa), ccaiafa@gmail.com, [Russel Greiner](https://rgreiner6.wixsite.com/greiner), russel rgreiner@ualberta.ca, [Farzane Aminmansour](https://www.linkedin.com/in/faminmansour/), aminmans@ualberta.ca
 
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
4. 

## Dependencies.
* [MatLab](http://www.mathworks.com/products/matlab/).

## Getting started.

### 1. [Download the repository](https://github.com/brain-life/encode).
* Download the repository from the TAR/ZIP files linked [here](https://github.com/framinmansour/Learning-Macroscopic-Brain-Connectomes-via-Group-Sparse-Factorization/archive/master.zip).
* UNZIP/UNTAR the file. 

### 2. [Run the demo_connectome_encoding code](/scripts/demos/demo_connectome_encoding.m).
Here you will learn about creating the tensor representation of a connectoms and perform basic operations such as identifying fascicles having a particular spatial orientation in a small voxel area. 
```
  >>  demo_connectome_encoding.m
```
### 6. [Run the demo_connectome_data_comparison code](/scripts/demos/demo_connectome_data_comparison.m).
This code reproduce Fig. 3 of the paper "Multidimensional encoding of brain connectomes", by C. Caiafa and F. Pestilli. 
```
  >>  demo_connectome_data_comparison.m
```
### 7. [Run the demo_virtual_lesion code](/scripts/demos/demo_virtual_lesion.m).
This code allows you to compute virtual lesions on a particular brain dataset and visualize particular major tracts together with their path-neighborhood, i.e. fascicles sharing same voxels. 
```
  >>  demo_virtual_lesion.m
```
### 8. [Run the demo_LiFE code](/scripts/demos/demo_LiFE.m).
This code allows you to compute compute the fascicles weights for two different tractography methods, probabilistic and deterministic tractographies, on a same brain. This is similar to the original LiFE demo in  https://github.com/francopestilli/life but here a full brain dataset is used. The optimization (fitting fascicles weights) runs in about 3 hours on a modern Intel processor with 8GB of RAM. This code has been tested with MatLab 2015b on Ubuntu 15+ and Mac OSX 10.11.
```
  >>  demo_LiFE.m
```


# How to Run these Codes?
1. Download datasets. For more information, please check the file __*save_compact_matrices*__.m in **data** folder (here is the [link](./data/save_compact_matrices.m "save_compact_matrices"))
2. If not exist yet, please create three folders named as **subsets**, **newsubsets** and **stage1** under **data** folder
3. Run  __*save_compact_matrices*__.m after choosing desired dataset. Check the comments in __*save_compact_matrices*__.m to know how to do so. Note that for each dataset this file is to be run only once unless the data is modified, as this file is used to both transfer the origin data to the data we will actually use and store the transfered data 
4. Use **mex** function of __*matlab*__ to compile __*forLoop*__.c in **mex** folder. For more information, please check this [link](https://www.mathworks.com/help/matlab/ref/mex.html). Remember to move the compiled file from **mex** folder to **codes** folder. 
5. __*brain*__.m,  __*stageOne*__.m and __*stageTwo*__.m are main files to do optimization and deliver results. Among them, __*brain*__.m is the one users will call the other two matlab files and save results; __*stageOne*__.m is to do jobs at stage 1 (initialize Phi with different strategies as warm start in the optimization at stage2) and to call __*stage2*__.m; __*stageTwo*__.m is to do jobs at stage 2 (optimize Phi). To know more about each file, please check their comments
6. The code files in **figures** folder are to draw figures on results. For now, let's just ignore them. 
