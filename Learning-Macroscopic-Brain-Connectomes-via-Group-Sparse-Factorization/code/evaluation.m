function evaluation(dataIndex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataIndex is the parameter to denote which dataset is to be processed
% 1 denotes Arcuate dataset
% 2 denotes ARC_SFL dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Of course, we need to add tensor toolbox
addpath tensor_toolbox

if nargin < 1
    dataIndex = 2;
end

%%
% load predicted phi for each dataset
if dataIndex == 1
    load('../experiments/dataset_one/Phi_11823_1_1_5_0.01_1000_1000.mat');
elseif dataIndex == 2
    load('../experiments/dataset_two/Phi_15033_2_1_5_0.01_1000_1000.mat');
else
    print('There is no dataset for the given data index!');
end

size(Phi_sp);

%%
% load Phi absorbing w. This is our expert Phi which we assumed is the
% ground truth.
if dataIndex == 1
    load('../data/subsets/compact_Phi_withw.mat');
elseif dataIndex == 2
    load('../data/newsubsets/compact_Phi_withw.mat');
end

size(Phi);


%%
% voxel_indices a list of indexes for compact Phi's voxels (e.g. [52,53,65,..])
% of course, if learning full voxels at stage 1, then voxels = voxel_indices
% voxels and voxel_indices could be the same if learning all voxels at
% stage1
if dataIndex == 1
    load('../data/subsets/voxel_indices.mat');
elseif dataIndex == 2
    load('../data/newsubsets/voxel_indices.mat');
end

size(voxel_indices);

%% Dimensions
% # of total orientations
Na = size(Phi, 1);
% # of total voxels (~190k)
Nv = size(Phi, 2);
% # of total fascicles
Nf = size(Phi, 3);
% # of selected voxels
%Nvb = size(Y, 2);

% err_not_sorted = 0;
% for v=1:Nv
%    err_per_v = immse(double(Phi(:,v,:)), double(Phi_sp(:,v,:)));
%    err_not_sorted = err_not_sorted + err_per_v;
% end
% err_not_sorted

% load Y
if dataIndex == 1
    load('../data/subsets/compact_Y.mat');
elseif dataIndex == 2
    load('../data/newsubsets/compact_Y.mat');
end


% load B (Dictionary)
if dataIndex == 1
    load('../data/subsets/B.mat');
elseif dataIndex == 2
    load('../data/newsubsets/B.mat');
end

Phi(:,3,:)
Phi_sp(:,3,:)


Phi_diff = Phi_sp - Phi;
Phi_diff(:,2,:);
size(Y);
Y(:,2);

ones_f = ones(Nf, 1);

Y_pred = ttm(ttv(Phi_sp, ones_f, 3), B, 1);
Y_exp = ttm(ttv(Phi, ones_f, 3), B, 1);

%Y_pred = B * Phi_pred_hat;
%Y_exp = B * Phi_exp_hat;

Y_pred(:,3)
Y_exp(:,3)


%% Sorting
f_ind_phi = [];
a_ind_phi = [];
for voxel=1:1
    index = find(Phi(:,voxel,:));
    f_ind_phi = [f_ind_phi index(:,2)];
    a_ind_phi = [a_ind_phi index(:,1)];
    
end

