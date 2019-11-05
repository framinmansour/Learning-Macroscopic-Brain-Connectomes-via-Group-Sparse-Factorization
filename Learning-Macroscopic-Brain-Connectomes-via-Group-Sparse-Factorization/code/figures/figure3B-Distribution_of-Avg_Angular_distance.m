%% Figure shows comparison of distribution of average angular distance of OMP vs Greedy

dataIndex = 1;
opt = 1;

if dataIndex == 1
    load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')
    if opt == 0
        load('../data/stage1/Phi_11823_1_2_5_0.01.mat');
        Phi_sp1 = Phi_sp;
        load('../data/stage1/Phi_11823_1_3_5_0.01_gs_0.01.mat');
        Phi_sp2 = Phi_sp;
    elseif opt == 1
        load('../experiments/dataset_one_new/Phi_11823_1_2_5_0.01_10_10.mat');
        Phi_sp1 = Phi_sp;
        load('../experiments/dataset_one_new/Phi_11823_1_3_5_0.01_10_10_gs_0.01.mat');
        Phi_sp2 = Phi_sp;
    end
    load('../data/subsets/compact_Phi.mat');
elseif dataIndex == 2
    load '../data/fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33'
    if opt == 0
        load('../data/stage1/Phi_15033_2_2_5_0.01.mat');
        Phi_sp1 = Phi_sp;
        load('../data/stage1/Phi_15033_2_3_5_0.01.mat');
        Phi_sp2 = Phi_sp;
    elseif opt == 1
        load('../experiments/dataset_one_new/Phi_11823_1_2_5_0.01_10_10.mat');
        Phi_sp1 = Phi_sp;
        load('../experiments/dataset_two_new/Phi_15033_1_3_5_0.01_10_10_gs_0.01.mat');
        Phi_sp2 = Phi_sp;
    end
    load('../data/newsubsets/compact_Phi.mat');
end


orient = fe.life.M.Atoms.orient';
orient_t = sptensor(orient);
angs_exp = ttt(Phi, orient_t, 1, 1);

angs_norm_exp = ttv(angs_exp .* angs_exp, [1,1,1]', 3);
[subs, vals] = find(angs_norm_exp);
angs_norm_exp = sptensor(subs, sqrt(vals), size(angs_norm_exp));

%get weighted average direction for each voxels
angs_our = ttt(Phi_sp1, orient_t, 1, 1);

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
%cdfplot(real(angs_diff.vals));
histogram(real(angs_diff.vals),'Normalization','probability','DisplayStyle','stairs');
hold on

%get weighted average direction for each voxels
angs_our = ttt(Phi_sp2, orient_t, 1, 1);

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

angs_diff.vals(imag(angs_diff.vals)~=0)
%cdfplot(real(angs_diff.vals));
histogram(real(angs_diff.vals),'Normalization','probability','DisplayStyle','stairs');

if dataIndex == 1
    xlabel('Angular Differences of Fascicle Orientations for Arcuate');
elseif dataIndex == 2
    xlabel('Angular Differences of Fascicle Orientations for ARC-SF');
end

ylabel('P(Angle Differences of Fascicle Orientations)');
legend('OMP', 'Greedy');