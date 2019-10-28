%%
% choose the dataset
% 1: Arcuate Fasciculus
% 2: ARC-SLF
dataIndex = 1;

%%
% load the fe_struct to get the vecotrs of orientations. 
if dataIndex == 1
    load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
    orients = fe.life.M.Atoms.orient;
elseif dataIndex == 2
    load('fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
    orients = fe.life.M.Atoms.orient;
end


%%
% load Phi_pred that you want to evaluate
% format of the file name: 
% Phi_' int2str(Nsv) '_' int2str(dataIndex) '_' int2str(stage1) '_' 
% int2str(numOrient) '_' num2str(regStage1) '.mat']
%Phi_filename = '../experiments/dataset_one/Phi_11823_1_1_5_0.01_1_1.mat';
%Phi_filename = '../data/stage1/Phi_11823_1_3_5_0.01.mat'
Phi_filename = '../data/stage1/Phi_OMP_all.mat'
%Phi_init_filename = '../data/stage1/Phi_11823_1_1_5_0.01.mat';
%load(Phi_init_filename)
load(Phi_filename);
Phi_sp = Phi;
%%
% load Phi absorbing w. This is our expert Phi which we assumed is the
% ground truth.
if dataIndex == 1
    load('../data/subsets/compact_Phi_withw.mat');
elseif dataIndex == 2
    load('../data/newsubsets/compact_Phi_withw.mat');
end

%%
[wangdist, angdist, countdist] = orientation_evaluation(Phi_sp(:,:,:), Phi(:,:,:), orients);

%%
% Get the average of weighted angular distance over all of the activated 
% voxels. 
[Nv, ~] = max(size(unique(Phi_sp.subs(:,2))));
fprintf('Num voxels: %d',Nv)
mean_wangdist_per_voxel = wangdist / Nv;
mean_wangdist_per_voxel

% Get the average of weighted angular distance over all of the activated 
% fascicles. 
[Nf, ~] = max(size(unique(Phi_sp.subs(:,3))));
fprintf('Num fascicles: %d', Nf)
mean_wangdist_per_fascicle = wangdist / Nf;
mean_wangdist_per_fascicle

%%
% Get the average of non-weighted angular distance over all of the activated 
% voxels.
mean_angdist_per_voxel = angdist / Nv;
mean_angdist_per_voxel

% Get the average of weighted angular distance over all of the activated 
% fascicles.
mean_angdist_per_fascicle = angdist / Nf;
mean_angdist_per_fascicle

%%
countdist
countdist / Nv