finalRes = [];
finalAng = [];
%load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')
load ('../data/fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33')

%load('../data/subsets/voxel_indices.mat');
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

%error
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
Y_hat = Y_hat.data;
Y_exp = Y(:,vind);
r = (Y_hat-Y_exp);
r = r.*r;
res = sum(r,1);
res = sum(sqrt(res))/length(voxel_indices);

finalRes = [finalRes,res];

% angle
%load('../data/subsets/compact_Phi_withw.mat');
load('../data/newsubsets/compact_Phi_withw.mat');

orient = fe.life.M.Atoms.orient';
orient_t = sptensor(orient);
angs_exp = ttt(Phi, orient_t, 1, 1);

angs_norm_exp = ttv(angs_exp .* angs_exp, [1,1,1]', 3);
[subs, vals] = find(angs_norm_exp);
angs_norm_exp = sptensor(subs, sqrt(vals), size(angs_norm_exp));


%get weighted average direction for each voxels
angs_our = ttt(Phi_sp, orient_t, 1, 1);

%get inner product between each directions that we want to compare
angs_prod = ttv(angs_our .* angs_exp, [1,1,1]', 3);

%norm of each vectors
angs_norm_our = ttv(angs_our .* angs_our, [1,1,1]', 3);
[subs, vals] = find(angs_norm_our);
angs_norm_our = sptensor(subs, sqrt(vals), size(angs_norm_our));

%calculate the angle with arccos
angs_diff = (angs_norm_our .* angs_norm_exp);
[subs, vals] = find(angs_diff);
angs_diff = sptensor(subs, acos(angs_prod(subs) ./ vals), size(angs_diff));

%if the angle is larger than 90 degrees, change to the complement angle.
flags = angs_diff > pi/2;
flags = find(flags);
if (sum(flags)~=0)
    angs_diff(flags) = pi - angs_diff(flags);
end;

angs_diff = sum(angs_diff.values)/(size(Phi,2)*size(Phi,3));
finalAng = [finalAng, angs_diff];

load('../data/newsubsets/Phi_final_15033_Large.mat');
%Phi_sp = Phi;

%error
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
Y_hat = Y_hat.data;
Y_exp = Y(:,vind);
r = (Y_hat-Y_exp);
r = r.*r;
res = sum(r,1);
res = sum(sqrt(res))/length(voxel_indices);
finalRes = [finalRes,res];

% angle
%load('../data/subsets/compact_Phi_withw.mat');
load('../data/newsubsets/compact_Phi_withw.mat');

orient = fe.life.M.Atoms.orient';
orient_t = sptensor(orient);
angs_exp = ttt(Phi, orient_t, 1, 1);

angs_norm_exp = ttv(angs_exp .* angs_exp, [1,1,1]', 3);
[subs, vals] = find(angs_norm_exp);
angs_norm_exp = sptensor(subs, sqrt(vals), size(angs_norm_exp));


%get weighted average direction for each voxels
angs_our = ttt(Phi_sp, orient_t, 1, 1);

%get inner product between each directions that we want to compare
angs_prod = ttv(angs_our .* angs_exp, [1,1,1]', 3);

%norm of each vectors
angs_norm_our = ttv(angs_our .* angs_our, [1,1,1]', 3);
[subs, vals] = find(angs_norm_our);
angs_norm_our = sptensor(subs, sqrt(vals), size(angs_norm_our));

%calculate the angle with arccos
angs_diff = (angs_norm_our .* angs_norm_exp);
[subs, vals] = find(angs_diff);
angs_diff = sptensor(subs, acos(angs_prod(subs) ./ vals), size(angs_diff));

%if the angle is larger than 90 degrees, change to the complement angle.
flags = angs_diff > pi/2;
flags = find(flags);
if (sum(flags)~=0)
    angs_diff(flags) = pi - angs_diff(flags);
end;

angs_diff = sum(angs_diff.values)/(size(Phi,2)*size(Phi,3));
finalAng = [finalAng, angs_diff];

load('../data/newsubsets/Phi_final_15033_Large_fabs.mat');
%Phi_sp = Phi;

%error
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
Y_hat = Y_hat.data;
Y_exp = Y(:,vind);
r = (Y_hat-Y_exp);
r = r.*r;
res = sum(r,1);
res = sum(sqrt(res))/length(voxel_indices);
finalRes = [finalRes,res];

% angle
%load('../data/subsets/compact_Phi_withw.mat');
load('../data/newsubsets/compact_Phi_withw.mat');

orient = fe.life.M.Atoms.orient';
orient_t = sptensor(orient);
angs_exp = ttt(Phi, orient_t, 1, 1);

angs_norm_exp = ttv(angs_exp .* angs_exp, [1,1,1]', 3);
[subs, vals] = find(angs_norm_exp);
angs_norm_exp = sptensor(subs, sqrt(vals), size(angs_norm_exp));


%get weighted average direction for each voxels
angs_our = ttt(Phi_sp, orient_t, 1, 1);

%get inner product between each directions that we want to compare
angs_prod = ttv(angs_our .* angs_exp, [1,1,1]', 3);

%norm of each vectors
angs_norm_our = ttv(angs_our .* angs_our, [1,1,1]', 3);
[subs, vals] = find(angs_norm_our);
angs_norm_our = sptensor(subs, sqrt(vals), size(angs_norm_our));

%calculate the angle with arccos
angs_diff = (angs_norm_our .* angs_norm_exp);
[subs, vals] = find(angs_diff);
angs_diff = sptensor(subs, acos(angs_prod(subs) ./ vals), size(angs_diff));

%if the angle is larger than 90 degrees, change to the complement angle.
flags = angs_diff > pi/2;
flags = find(flags);
if (sum(flags)~=0)
    angs_diff(flags) = pi - angs_diff(flags);
end;

angs_diff = sum(angs_diff.values)/(size(Phi,2)*size(Phi,3));
finalAng = [finalAng, angs_diff];

plot([finalRes],[finalAng]);

finalRes
finalAng

ylabel('Average Angle Differences of Fascicle Directions');
xlabel('Average Error from Expert Phi');
%legend('Model2(OMP)', 'Model2(Greedy)');