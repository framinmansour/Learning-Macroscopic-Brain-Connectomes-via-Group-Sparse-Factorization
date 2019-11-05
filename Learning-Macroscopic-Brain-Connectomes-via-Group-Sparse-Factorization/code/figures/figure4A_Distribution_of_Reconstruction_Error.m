%% Demonstrates a comparison between the distributions of reconstructio error of signal matrix Y after optimization 
%% for Phi_predicts initialized with each of the OMP vs Greedy vs GroundTruth in stage 1. 
% Set parameters: 
% bins: number of bins in histogram
% dataIndex: datset 1 = Arcuate, dataser 2 = ARC-SLF
% opt: if true, considers Phi after optimization. If false, uses
% initialized Phi
bins = 50;
dataIndex = 2;
opt = true;

if dataIndex == 1
    load('../data/subsets/voxel_indices.mat');
    load('../data/subsets/B.mat');
    load('../data/subsets/compact_Y.mat');
    load('../data/subsets/weights.mat');
    load('../data/subsets/voxel_vicinity.mat');
elseif dataIndex == 2
    load('../data/newsubsets/voxel_indices.mat');
    load('../data/newsubsets/B.mat');
    load('../data/newsubsets/compact_Y.mat');
    load('../data/newsubsets/weights.mat');
    load('../data/newsubsets/voxel_vicinity.mat');
end


%% For OMP:
if dataIndex == 1
    load('../data/stage1/Phi_11823_1_2_5_0.01.mat');
    if opt
        load('../experiments/dataset_one_new/Phi_11823_1_2_5_0.01_10_10.mat');
    end
elseif dataIndex == 2
    load('../data/stage1/Phi_15033_2_2_5_0.01.mat');
    if opt
        load('../experiments/dataset_two_new/Phi_15033_2_2_5_0.01_10_10.mat');
    end
end

a = zeros(1,size(voxel_vicinity,2));
a(voxel_indices) = true;
voxels = logical(a);
v_compact_ind = zeros(size(voxels));
v_compact_ind(voxels) = 1:nnz(voxels);

vind = v_compact_ind(vlist);

Phi = Phi_sp(:,vind,:);

ones_f = ones(size(w,1),1);
temp = ttv(Phi,ones_f,3);
Y_hat = ttm(temp,B,1);% get the fitted Y
%temp = ttm(Phi,B,1);
%Y_hat = ttv(temp,ones_f,3);% get the fitted Y
Y_hat = Y_hat.data;
Y_exp = Y(:,vind);

r = (Y_hat-Y_exp);

r = r.*r;

res = sum(r,1);

res = sqrt(res);

histogram(res,bins,'Normalization','probability','DisplayStyle','stairs','LineWidth',1.5);
set(gca,'xscale','log')

hold on

%% For Greedy

if dataIndex == 1
    load('../data/stage1/Phi_11823_1_3_5_0.01_gs_0.01.mat');
    if opt
        load('../experiments/dataset_one_new/Phi_11823_1_3_5_0.01_10_10_gs_0.01.mat');
    end
elseif dataIndex == 2
    load('../data/stage1/Phi_15033_2_3_5_0.01_gs_0.01.mat');
    if opt
        load('../experiments/dataset_two_new/Phi_15033_2_3_5_0.01_10_10_gs_0.01.mat');
    end
end

a = zeros(1,size(voxel_vicinity,2));
a(voxel_indices) = true;
voxels = logical(a);
v_compact_ind = zeros(size(voxels));
v_compact_ind(voxels) = 1:nnz(voxels);

vind = v_compact_ind(vlist);

Phi = Phi_sp(:,vind,:);

ones_f = ones(size(w,1),1);
temp = ttv(Phi,ones_f,3);
Y_hat = ttm(temp,B,1);% get the fitted Y
%temp = ttm(Phi,B,1);
%Y_hat = ttv(temp,ones_f,3);% get the fitted Y
Y_hat = Y_hat.data;
Y_exp = Y(:,vind);

r = (Y_hat-Y_exp);

r = r.*r;

res = sum(r,1);

res = sqrt(res);

histogram(res,bins,'Normalization','probability','DisplayStyle','stairs','LineWidth',1.25);
set(gca,'xscale','log')
hold on

if dataIndex == 1
    load('../data/stage1/Phi_11823_1_1_5_0.01.mat');
    if opt
        load('../experiments/dataset_one_new/Phi_11823_1_1_5_0.01_10_10.mat');
    end
elseif dataIndex == 2
    load('../data/stage1/Phi_15033_2_1_5_0.01.mat');
    if opt
        load('../experiments/dataset_two/Phi_15033_2_1_5_0.01_10_10.mat');
    end
end

a = zeros(1,size(voxel_vicinity,2));
a(voxel_indices) = true;
voxels = logical(a);
v_compact_ind = zeros(size(voxels));
v_compact_ind(voxels) = 1:nnz(voxels);

vind = v_compact_ind(vlist);

Phi = Phi_sp(:,vind,:);

ones_f = ones(size(w,1),1);
temp = ttv(Phi,ones_f,3);
Y_hat = ttm(temp,B,1);% get the fitted Y
%temp = ttm(Phi,B,1);
%Y_hat = ttv(temp,ones_f,3);% get the fitted Y
Y_hat = Y_hat.data;
Y_exp = Y(:,vind);

r = (Y_hat-Y_exp);

r = r.*r;

res = sum(r,1);

res = sqrt(res);

histogram(res,bins,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
set(gca,'xscale','log')
if dataIndex == 1
    xlabel('Model Error in reconstructing Y after optimization for Arcuate');
elseif dataIndex == 2
    xlabel('Model Error in reconstructing Y after optimization for ARC-SLF');
end
ylabel('P(Reconstruction Error)');
legend('OMP', 'Greedy', 'GroundTruth');
%legend('Stage1(GD)', 'Stage2(GD)', 'Stage2(GD) more iterations');