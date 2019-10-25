function save_compact_matrices(dataIndex)
% dataIndex is the parameter to denote which dataset is to be processed
% 1 denotes Arcuate dataset
% 2 denotes ARC_SFL dataset

% Of course, we need to add tensor toolbox
addpath ../code/tensor_toolbox

% If have not done so, please download Arcuate data set from 
% https://iu.app.box.com/file/171262214861
% and ARC_SFL dataset from
% https://iu.app.box.com/file/189177566190
% and move them to the data folder. 

% Load dataset
if dataIndex == 1
    % this is to load Arcuate dataset
    load ../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat
elseif dataIndex == 2
    % this is to load ARC_SFL dataset
    load ../data/fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat
end

Y = fe.life.diffusion_signal_img';
ly = fe.life.predicted_signal_demean;
Y = reshape(ly, size(Y)); % get the predicted Y of the correct size

% Dimension of Phi = # of atoms x # of voxels x # of fascicles
Phi = fe.life.M.Phi;
Phi = sptensor(Phi.subs, Phi.vals, Phi.size);
w = fe.life.fit.weights;

PhiTilde = ttv(Phi, w, 3);
PhiTilde = spmatrix(PhiTilde);
flag = PhiTilde ~= 0;

% if zero, then the whole column is zero (this is voxel)
% then that voxel is not necessary
% in both datasets, we find lots of unnecessary voxels
sc = sum(flag, 1);

% if zero, then the whole row is zero (this is orientation)
% then that orientation is not necessary
% in either dataset, we do not find any unnecessary orientations
sr = sum(flag, 2);

% the index of voxels which are not all zeros (namely necessary) in fitted Y
% we only care about those voxels
voxel_indices = find(sc);
% Store the indices of necessary voxels (nonzero ones)
if dataIndex == 1
    % store it in the folder ../data/subsets/
    save('../data/subsets/voxel_indices.mat', 'voxel_indices');
elseif dataIndex == 2
    % store it in the folder ../data/newsubsets/
    save('../data/newsubsets/voxel_indices.mat', 'voxel_indices');
end

Y = Y(:,voxel_indices);
% Store compact Y
if dataIndex == 1 
    % Store compact data from Arcuate dataset in the folder ../data/subsets/
    save('../data/subsets/compact_Y.mat', 'Y', 'voxel_indices');
elseif dataIndex == 2
    % Store compact data from ARC_SLF dataset in the folder ../data/newsubsets/
    save('../data/newsubsets/compact_Y.mat', 'Y', 'voxel_indices');
end

% if some w is 0, corresponding fascicle is not necessary
% like orientations, in both datasets, we find no fascicle is unnecessary
f = w ~= 0;
f = find(f);

Phi = Phi(:,voxel_indices,f);
if dataIndex == 1
    save('../data/subsets/compact_Phi.mat', 'Phi');
elseif dataIndex == 2
    save('../data/newsubsets/compact_Phi.mat', 'Phi', 'voxel_indices');
end

w = w(f);
if dataIndex == 1
    save('../data/subsets/weights.mat', 'w');
elseif dataIndex == 2
    save('../data/newsubsets/weights.mat', 'w', 'voxel_indices');
end

B = fe.life.M.DictSig;
if dataIndex == 1
    save('../data/subsets/B.mat', 'B');
elseif dataIndex == 2
    save('../data/newsubsets/B.mat', 'B');
end

voxel_vicinity = fe.life.M.Voxels.vicinity;
if dataIndex == 1
    save('../data/subsets/voxel_vicinity', 'voxel_vicinity');
elseif dataIndex == 2
    save('../data/newsubsets/voxel_vicinity', 'voxel_vicinity');
end

atom_vicinity = fe.life.M.Atoms.vicinity;
if dataIndex == 1
    save('../data/subsets/atom_vicinity', 'atom_vicinity');
elseif dataIndex == 2
    save('../data/newsubsets/atom_vicinity', 'atom_vicinity');
end
 
% absorb weights into Phi
for i = 1 : size(Phi,3)
    Phi(:,:,i) = Phi(:,:,i) * w(i);
end
if dataIndex == 1
    save('../data/subsets/compact_Phi_withw.mat', 'Phi');
elseif dataIndex == 2
    save('../data/newsubsets/compact_Phi_withw.mat', 'Phi');
end
