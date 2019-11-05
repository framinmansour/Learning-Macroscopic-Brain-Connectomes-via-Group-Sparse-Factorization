%% This code counts the average of number of orientations missing in the candidate set per voxel for initialized Phi over different values of k (candidate set size). 
% % Of course, we need to add tensor toolbox
% addpath tensor_toolbox
% 
% %%
% % choose the dataset
% % 1: Arcuate Fasciculus
% % 2: ARC-SLF
% dataIndex = 1;
% max_orient = 62;
% steps = 10;
% 
% %%
% % load the fe_struct to get the vecotrs of orientations. 
% if dataIndex == 1
%     load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
%     orients = fe.life.M.Atoms.orient;
%     Nsv = 11823;
% elseif dataIndex == 2
%     load('../data/fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
%     orients = fe.life.M.Atoms.orient;
%     Nsv = 15033;
% end
% 
% 
% %%
% % load Phi absorbing w. This is our expert Phi which we assumed is the
% % ground truth.
% if dataIndex == 1
%     load('../data/subsets/compact_Phi_withw.mat');
% elseif dataIndex == 2
%     load('../data/newsubsets/compact_Phi_withw.mat');
% end
% 
% Phi_exp = Phi(:,:,:);
% Na_exp = size(Phi_exp, 1);
% Nf_exp = size(Phi_exp, 3);
% 
% % all of the atoms activated in Phi_exp
% a_exp = unique(Phi_exp.subs(:,1));
% 
% ones_f = ones(Nf_exp, 1);
% 
% 
% %% For Greedy
% stg1 = 3;
% GD_angdiff = [];
% for k = 2:steps:max_orient
%     brain(dataIndex, Nsv, k, stg1);
%     % load Phi_pred that you want to evaluate
%     % format of the file name: 
%     % Phi_' int2str(Nsv) '_' int2str(dataIndex) '_' int2str(stage1) '_' 
%     % int2str(numOrient) '_' num2str(regStage1) '.mat']
%     if dataIndex == 1
%     load(['../data/stage1/Phi_11823_1_3_', num2str(k), '_0.01.mat']);
%     elseif dataIndex == 2
%     load(['../data/stage1/Phi_15033_2_3_', num2str(k), '_0.01.mat']);
%     end
%     Phi_pred = Phi_sp;
%     Na_pred = size(Phi_pred, 1);
%     Nv = size(Phi_pred, 2);
%     Nf_pred = size(Phi_pred, 3);
%     % all of the atoms activated in Phi_pred
%     a_pred = unique(Phi_pred.subs(:, 1));    
%     
%     countdist = 0.0;
%     Total_num_orients_in_Phiexp = 0.0;
%     for vp = 1:Nv
%         if mod(vp,10000) == 0
%             fprintf("Voxel %d of %d\n", vp, Nv);
%         end
% 
%         % Find non-zero oriantations in each voxel
%         Phi_pred_nz_atoms = find(Phi_pred(:,vp,:));
% 
%         % Make an orientation set for each voxel. The goal is to find different
%         % types of orientations in each voxel for Phi_pred.
%         apredset_per_v = unique(Phi_pred_nz_atoms(:,1));
%         %fprintf('\n atoms in voxel %d of pred', vp);
%         %apredset_per_v
% 
%         % Find the non-zero orientations in voxel vp of Phi_exp
%         Phi_exp_nz_atoms = find(Phi_exp(:,vp,:));
% 
%         % Make a set of active orientations for voxel vp
%         aexpset_per_v = unique(Phi_exp_nz_atoms(:,1));
%         %fprintf('\n atoms in voxel %d of exp', vp);
%         %aexpset_per_v
%         Total_num_orients_in_Phiexp = Total_num_orients_in_Phiexp + size(aexpset_per_v);
% 
%         % The intersection of orientations in Phi_exp and Phi_pred for the same
%         % voxel shows that how close could we potentially get to the ground 
%         % truth.
%         mutual_orients = intersect(apredset_per_v, aexpset_per_v);
% 
%         % We want the distance, so we need the numeber of orientations which
%         % has not been included in Phi_pre(:,vp,:) but included in Phi_exp(:,vp,:)
%         % The ideal case is that the length of mutual orientations be exactly
%         % the same as the activated orientations in Phi_exp(:,vp,:)
%         countdist = countdist + (length(aexpset_per_v) - length(mutual_orients));
%     end
%     GD_angdiff = [GD_angdiff, countdist/Nv];
%     GD_angdiff
% end
% 
% %% For OMP
% stg1 = 2;
% OMP_angdiff = [];
% for k = 2:steps:max_orient
%     brain(dataIndex, Nsv, k, stg1);
%     % load Phi_pred that you want to evaluate
%     % format of the file name: 
%     % Phi_' int2str(Nsv) '_' int2str(dataIndex) '_' int2str(stage1) '_' 
%     % int2str(numOrient) '_' num2str(regStage1) '.mat']
%     if dataIndex == 1
%         load(['../data/stage1/Phi_11823_1_2_', num2str(k),'_0.01.mat']);
%     elseif dataIndex == 2
%         load(['../data/stage1/Phi_15033_2_2_', num2str(k), '_0.01.mat']);
%     end
%     Phi_pred = Phi_sp;
%     Na_pred = size(Phi_pred, 1);
%     Nv = size(Phi_pred, 2);
%     Nf_pred = size(Phi_pred, 3);
%     % all of the atoms activated in Phi_pred
%     a_pred = unique(Phi_pred.subs(:, 1));
% 
%     countdist = 0.0;
%     Total_num_orients_in_Phiexp = 0.0;
%     for vp = 1:Nv
%         if mod(vp,10000) == 0
%             fprintf("Voxel %d of %d\n", vp, Nv);
%         end
% 
%         % Find non-zero oriantations in each voxel
%         Phi_pred_nz_atoms = find(Phi_pred(:,vp,:));
% 
%         % Make an orientation set for each voxel. The goal is to find different
%         % types of orientations in each voxel for Phi_pred.
%         apredset_per_v = unique(Phi_pred_nz_atoms(:,1));
%         %fprintf('\n atoms in voxel %d of pred', vp);
%         %apredset_per_v
% 
%         % Find the non-zero orientations in voxel vp of Phi_exp
%         Phi_exp_nz_atoms = find(Phi_exp(:,vp,:));
% 
%         % Make a set of active orientations for voxel vp
%         aexpset_per_v = unique(Phi_exp_nz_atoms(:,1));
%         %fprintf('\n atoms in voxel %d of exp', vp);
%         %aexpset_per_v
%         Total_num_orients_in_Phiexp = Total_num_orients_in_Phiexp + size(aexpset_per_v);
% 
%         % The intersection of orientations in Phi_exp and Phi_pred for the same
%         % voxel shows that how close could we potentially get to the ground 
%         % truth.
%         mutual_orients = intersect(apredset_per_v, aexpset_per_v);
% 
%         % We want the distance, so we need the numeber of orientations which
%         % has not been included in Phi_pre(:,vp,:) but included in Phi_exp(:,vp,:)
%         % The ideal case is that the length of mutual orientations be exactly
%         % the same as the activated orientations in Phi_exp(:,vp,:)
%         countdist = countdist + (length(aexpset_per_v) - length(mutual_orients));
%     end
%     OMP_angdiff = [OMP_angdiff, countdist/Nv];
%     OMP_angdiff
% end
% 
% 
% x = 2:steps:max_orient;
% if dataIndex == 1
%     save('../data/figures_mats/Arcuate_num_mismatch_orients', 'x', 'OMP_angdiff', 'GD_angdiff')
% elseif dataIndex == 2
%     save('../data/figures_mats/ARC-SLF_num_mismatch_orients', 'x', 'OMP_angdiff', 'GD_angdiff')
% end

%% After you ran the first part, uncomment this part and comment last part to get the graph
%load('../data/figures_mats/Arcuate.mat')
%base = 0.0086;
dataIndex = 1;
if dataIndex == 1
    load('../data/figures_mats/Arcuate_num_mismatch_orients_gs1.mat')
elseif dataIndex == 2
    load('../data/figures_mats/ARC-SLF_num_mismatch_orients_gs1.mat')
end

% For Arcuate
%GD_angdiff = [5.2593, 3.1960, 2.1199, 1.5493, 1.2397, 1.0638, 0.9590];
GD_angdiff = [5.1298, 3.6061, 2.6846, 2.1215, 1.7667, 1.5484, 1.4022];

% For ARC-SL
%GD_angdiff = [5.1471, 3.6299, 2.9789, 2.0570, 1.6350, 1.3939, 1.1758]

%bar(x, OMP_angdiff);
semilogy(x, OMP_angdiff, 'LineWidth',2);
hold on
semilogy(x, GD_angdiff, 'LineWidth',2);
%bar(x, GD_angdiff);

xlabel('Number of selected orientations');
ylabel('Average Num of Missing Orientations per Voxel - Arcuate');
legend('OMP', 'Greedy');