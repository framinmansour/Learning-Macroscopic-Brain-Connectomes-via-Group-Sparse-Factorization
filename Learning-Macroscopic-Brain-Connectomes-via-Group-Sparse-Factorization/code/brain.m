function [pars] = brain(dataIndex, Nsv, k, stg1)
% Nsv is the maximum number of voxels to be processed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataIndex is the parameter to denote which dataset is to be processed
% 1 denotes Arcuate dataset
% 2 denotes ARC_SFL dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Of course, we need to add tensor toolbox
addpath tensor_toolbox

% For the data used in the folder ../data/subsets
% (those are compact data from Arcuate dataset)
% full Nsv is 11823
numberArcuate = 11823;

% For the data used in the folder ../data/newsubsets
% (those are compact data from ARC_SLF dataset)
% full Nsv is 15033
numberARC_SLF = 15033;

% if do not offer the parameters, set default parameters
% this is only for safety, and please always offer dataIndex and Nsv
if nargin < 1
    dataIndex = 1;
    Nsv = numberArcuate;
elseif nargin < 2
    if dataIndex == 1
        Nsv = numberArcuate;
    elseif dataIndex == 2
        Nsv = numberARC_SLF;
    else
        dataIndex = 1;
        Nsv = numberArcuate;
    end
end
%% Settings of Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calStage1 decides if re-calculate stage1
calStage1 = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stage1 decides the method used in stage 1
% 1 is Model 1 (just use the expert Phi as initialization for Stage2)
% 2 is OMP for the selection of orientations
% 3 is GreedyDirections (main method in the paper)
stage1 = stg1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 not go through stage 2
% 1 go through stage 2
stage2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximum # of orientations set nonzero for each fascicle at stage1
% Only valid for stage1 == 2 or 3
numOrient = k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is used as the regularization
% parameter for assigning the warm start
regStage1 = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the max iterations for optimization at stage2
maxIter = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the max iterations for optimization at stage2
% the relative tolerance for difference of last two objective value
% this is for one ending conditions for optimization at stage2
relTolerance = 1e-10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is stepsize for optimization at stage2
stepSize = 1e-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambdaL1 is the regularization parameter for L1 norm of optimization at stage2
% lambdaGroup is regularization parameter for group regularization of optimization at stage2
lambdaL1 = 10;
lambdaGroup = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename of storing results
saveNameStage1 = ['../data/stage1/Phi_' int2str(Nsv) '_' int2str(dataIndex) '_' int2str(stage1) '_' int2str(numOrient) '_' num2str(regStage1) '.mat'];
saveNamePhi1 = ['../data/subsets/Phi_' int2str(Nsv) '_' int2str(dataIndex) '_' int2str(stage1) '_' int2str(numOrient) '_' num2str(regStage1) '_' num2str(lambdaL1) '_' num2str(lambdaGroup) '.mat'];
saveNamePhi2 = ['../data/newsubsets/Phi_' int2str(Nsv) '_' int2str(dataIndex) '_' int2str(stage1) '_' int2str(numOrient) '_' num2str(regStage1) '_' num2str(lambdaL1) '_' num2str(lambdaGroup) '.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% truncate denotes if truncating small values to be zero
% relThres sets the relative threshold values for truncation
relThres = 0.001;
truncate = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximal operator for l1 norm part
% the whole optimization is like proximal gradient descent by using sub-gradient descent
% for group regularization part
useProximal = 1;

% sparse coding parameters
pars.dataIndex = dataIndex;
pars.calStage1 = calStage1;
pars.stage1 = stage1;
pars.stage2 = stage2;
pars.numOrient = numOrient;
pars.regStage1 = regStage1;
pars.maxIter = maxIter;
pars.relTolerance = relTolerance;
pars.stepSize = stepSize;
pars.lambdaL1 = lambdaL1;
pars.lambdaGroup = lambdaGroup;
pars.saveNameStage1 = saveNameStage1;
pars.saveNamePhi1 = saveNamePhi1;
pars.saveNamePhi2 = saveNamePhi2;
pars.relThres = relThres;
pars.truncate = truncate;
pars.useProximal = useProximal;
%% Load Data
% load B (Dictionary)
if pars.dataIndex == 1
    load('../data/subsets/B.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/B.mat');
end

% load w (linear parameters)
if pars.dataIndex == 1
    load('../data/subsets/weights.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/weights.mat');
end
% since Phi absorbs w
% hence, w is a vector of ones 
% (only used to sum up the dimension)
w = ones(size(w));

% load Y
if pars.dataIndex == 1
    load('../data/subsets/compact_Y.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/compact_Y.mat');
end

% if Nsv is set even larger than the full number of voxels
% then set it to the full number
if (Nsv > size(Y,2))
    Nsv = size(Y,2);
end

% load Phi absorbing w
if pars.dataIndex == 1
    load('../data/subsets/compact_Phi_withw.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/compact_Phi_withw.mat');
end
% Phi_sp will store the returned Phi
% size is the compact Phi
Phi_sp = sptensor([], [], size(Phi));

% neighborhood index for voxels
if pars.dataIndex == 1
    load('../data/subsets/voxel_vicinity.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/voxel_vicinity.mat');
end

% neighborhood index for orientations
if pars.dataIndex == 1
    load('../data/subsets/atom_vicinity.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/atom_vicinity.mat');
end

% voxel index (of necessary voxels) from original dataset
if pars.dataIndex == 1
    load('../data/subsets/voxel_indices.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/voxel_indices.mat');
end


%% Call Stage1 (and Stage2)

% let's go to stage1: initialize Phi for optimization
% Note that stage2 (the optimization for Phi) will be called inside stage1
% so returned Phi_sp will be already the final one through these two stages
[Phi_sp, ss] = stageOne(Y, B, Phi_sp, Phi, w, voxel_indices, voxel_vicinity, atom_vicinity, Nsv, pars);

% not go through stage 2
if stage2 == 0
    return;
end

% store returned Phi
if pars.dataIndex == 1
    %filename = ['../data/subsets/Phi_' int2str(Nsv) '_' char(datetime('now','Format','yyyy_MM_dd_HHmmss')) '.mat'];
    filename = saveNamePhi1;
elseif pars.dataIndex == 2
    %filename = ['../data/newsubsets/Phi_' int2str(Nsv) '_' char(datetime('now','Format','yyyy_MM_dd_HHmmss')) '.mat'];
    filename = saveNamePhi2;
end
% save the result as Phi_sp
save(filename, 'Phi_sp', 'ss', '-v7.3');
% 
% %% visualization
% load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')
% load('../data/stage1/Phi_50_1_3_5_0.01.mat')
% load('../data/subsets/compact_Phi_withw.mat');
% hs = visualize_Phi_with_fiber(fe, Phi, Phi_sp, pars);
% hs