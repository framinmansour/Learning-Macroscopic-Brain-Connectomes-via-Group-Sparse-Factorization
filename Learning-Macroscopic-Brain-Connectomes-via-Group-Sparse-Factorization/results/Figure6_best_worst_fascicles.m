%% Calls displayPhi function and passes the indices of best and worst fascicles based on average angular difference of them with ground truth. 
% 1 is Arcuate and 2 is ARC-SLF
dataIndex = 2;
% 2 is OMP and 3 is Greedy
stg1 = 2;
%before or after optimization
opt = true;

f_num = 5;


if dataIndex == 1
    load('../codes/data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')
    load('../codes/data/subsets/compact_Phi.mat');
    if ~opt
        if stg1 == 2
            load('../codes/data/stage1/Phi_11823_1_2_5_0.01.mat');
        elseif stg1 == 3
            load('../codes/data/stage1/Phi_11823_1_3_5_0.01_gs0.01.mat');
        end
    else
        if stg1 == 2
            load('../codes/experiments/dataset_one_new/Phi_11823_1_2_5_0.01_10_10.mat');
        elseif stg1 == 3
            %load('../codes/experiments/dataset_one_new/Phi_11823_1_1_5_0.01_0_0.1_5extra_l0_g0.1_1e3.mat');
            load('../codes/experiments/dataset_one_new/Phi_11823_1_5_1_0.01_0_0.0001.mat');
        end
    end
elseif dataIndex == 2
    load '../codes/data/fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33'
    load('../codes/data/newsubsets/compact_Phi.mat');
    if ~opt
        if stg1 == 2
            load('../codes/data/stage1/Phi_15033_2_2_5_0.01.mat');
        elseif stg1 == 3
            load('../codes/data/stage1/Phi_15033_2_3_5_0.01_gs_0.01.mat');
        end
    else
        if stg1 == 2
            load('../codes/experiments/dataset_two_new/Phi_15033_2_2_5_0.01_10_10.mat');
        elseif stg1 == 3
            load('../codes/experiments/dataset_two_new/Phi_15033_2_1_5_0.01_0_0.0001.mat');
        end
    end
end

orient = fe.life.M.Atoms.orient';

orient_t = sptensor(orient);
angs_exp = ttt(Phi, orient_t, 1, 1);
angs_norm_exp = ttv(angs_exp .* angs_exp, [1,1,1]', 3);
[subs, vals] = find(angs_norm_exp);
angs_norm_exp = sptensor(subs, sqrt(vals), size(angs_norm_exp));

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
end

Na = size(Phi_sp,1);
Nv = size(Phi_sp,2);
Nf = size(Phi_sp,3);
mean_angdiff_per_f = mean(double(angs_diff(:,:)));
[angdiff, aindx] = sort(mean_angdiff_per_f);

best_f = aindx(1:f_num);
best_f
worst_f = aindx(Nf-f_num+1:Nf);
worst_f

lineSpec1 = '-';
displayPhi(Phi_sp(:,vind,best_f),fe.life.M.Atoms.orient,fe.roi.coords(vlist,:), lineSpec1, false);
%displayPhi(Phi_sp(:,vind,worst_f),fe.life.M.Atoms.orient,fe.roi.coords(vlist,:), lineSpec1,false)

if dataIndex == 1
    load('../codes/data/stage1/Phi_11823_1_1_5_0.01.mat');
elseif dataIndex == 2 
    
    load('../codes/data/stage1/Phi_15033_2_1_5_0.01.mat');
end

lineSpec2 = '--';
displayPhi(Phi_sp(:,vind,best_f),fe.life.M.Atoms.orient,fe.roi.coords(vlist,:),lineSpec2, false);
%displayPhi(Phi_sp(:,vind,worst_f),fe.life.M.Atoms.orient,fe.roi.coords(vlist,:), lineSpec2,false)

