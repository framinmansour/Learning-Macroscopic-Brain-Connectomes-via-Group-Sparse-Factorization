%load('../data/newsubsets/voxel_indices.mat');
%load('../data/subsets/B.mat');
%load('../data/subsets/compact_Y.mat');
%load('../data/subsets/weights.mat');
%load('../data/subsets/voxel_vicinity.mat');

load('../data/newsubsets/voxel_indices.mat');
load('../data/newsubsets/B.mat');
load('../data/newsubsets/compact_Y.mat');
load('../data/newsubsets/weights.mat');
load('../data/newsubsets/voxel_vicinity.mat');

load('../data/stage1/Phi_GD_Large.mat');
Phi_sp = Phi;

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

histogram(res,100,'Normalization','probability','DisplayStyle','stairs');
hold on

load('../data/newsubsets/Phi_final_15033_lei_Large.mat');

Phi_sp = Phi_sp;

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

histogram(res,100,'Normalization','probability','DisplayStyle','stairs');
load('../data/newsubsets/Phi_final_15033_lei_Large_5.mat');

Phi_sp = Phi_sp;

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

histogram(res,100,'Normalization','probability','DisplayStyle','stairs');
xlabel('Model Error in constructing Y');
ylabel('P(Error)');
%legend('Model1','Model2(OMP)', 'Model3(Greedy)');
legend('Stage1(GD)', 'Stage2(GD)', 'Stage2(GD) more iterations');





