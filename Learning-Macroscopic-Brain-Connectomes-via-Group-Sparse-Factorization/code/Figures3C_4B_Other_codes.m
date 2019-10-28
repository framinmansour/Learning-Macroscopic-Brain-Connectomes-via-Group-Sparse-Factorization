% %%%
% % choose the dataset
% % 1: Arcuate Fasciculus
% % 2: ARC-SLF
% dataIndex = 1;
% 
% %%
% % load the fe_struct to get the vecotrs of orientations. 
% if dataIndex == 1
%     load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
%     orients = fe.life.M.Atoms.orient;
% elseif dataIndex == 2
%     load('fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
%     orients = fe.life.M.Atoms.orient;
% end
% 
% %%
% % load Phi absorbing w. This is our expert Phi which we assumed is the
% % ground truth.
% if dataIndex == 1
%     load('../data/subsets/compact_Phi_withw.mat');
% elseif dataIndex == 2
%     load('../data/newsubsets/compact_Phi_withw.mat');
% end
% %nnz(Phi)
% 
% Phi_filename = '../data/stage1/Phi_11823_1_3_5_0.01.mat';
% load(Phi_filename);
% 
% Phi_pred = Phi_sp(:,vind,:);
% Phi_exp = Phi(:,vind,:);
% Nv = size(Phi_pred, 2);
% countdist = 0.0;
% count_diff = 0.0;
% count_captured = 0.0;
% Total_num_orients_in_Phiexp = 0.0;
% for vp = 1:Nv
%     if mod(vp,100) == 0
%         fprintf("Voxel %d of %d\n", vp, Nv);
%     end
%     
%     % Find non-zero oriantations in each voxel
%     Phi_pred_nz_atoms = find(Phi_pred(:,vp,:));
%     
%     % Make an orientation set for each voxel. The goal is to find different
%     % types of orientations in each voxel for Phi_pred.
%     apredset_per_v = unique(Phi_pred_nz_atoms(:,1));    
%     
%     % Find the non-zero orientations in voxel vp of Phi_exp
%     Phi_exp_nz_atoms = find(Phi_exp(:,vp,:));
%     
%     % Make a set of active orientations for voxel vp
%     aexpset_per_v = unique(Phi_exp_nz_atoms(:,1));
%     
%     Total_num_orients_in_Phiexp = Total_num_orients_in_Phiexp + size(Phi_exp_nz_atoms);
%     
%     % The intersection of orientations in Phi_exp and Phi_pred for the same
%     % voxel shows that how close could we potentially get to the ground 
%     % truth.
%     mutual_orients = intersect(apredset_per_v, aexpset_per_v);
%     count_captured = count_captured + length(mutual_orients);
%     
%     diff = length(setdiff(aexpset_per_v, apredset_per_v));
%     diff;
%     % We want the distance, so we need the numeber of orientations which
%     % has not been included in Phi_pre(:,vp,:) but included in Phi_exp(:,vp,:)
%     % The ideal case is that the length of mutual orientations be exactly
%     % the same as the activated orientations in Phi_exp(:,vp,:)
%     countdist = countdist + (length(aexpset_per_v) - length(mutual_orients));
%     count_diff = count_diff + diff;
% end
% countdist
% count_diff
% count_captured
% fprintf('Total number of active ground truth orientations %d \n', Total_num_orients_in_Phiexp);




% %%
% % load Phi_pred that you want to evaluate
% % format of the file name: 
% % Phi_' int2str(Nsv) '_' int2str(dataIndex) '_' int2str(stage1) '_' 
% % int2str(numOrient) '_' num2str(regStage1) '.mat']
% %Phi_filename = '../experiments/dataset_one/Phi_11823_1_1_5_0.01_1_1.mat';
% Phi_filename = '../data/stage1/Phi_200_1_2_5_0.01.mat'
% %Phi_init_filename = '../data/stage1/Phi_11823_1_1_5_0.01.mat';
% %load(Phi_init_filename)
% load(Phi_filename);
% Phi_omp = Phi_sp(:,vind,:);
% Phi_filename = '../data/stage1/Phi_200_1_3_5_0.01.mat'
% %Phi_init_filename = '../data/stage1/Phi_11823_1_1_5_0.01.mat';
% %load(Phi_init_filename)
% load(Phi_filename);
% Phi_greedy = Phi_sp(:,vind,:);



% Nv = size(Phi_greedy, 2);
% for vp = 1:Nv
%     if mod(vp,100) == 0
%         fprintf("*");
%     end
%         
%     % Find non-zero oriantations in each voxel
%     Phi_omp_nz_atoms = find(Phi_omp(:,vp,:));
% 
%     % Make an orientation set for each voxel. The goal is to find different
%     % types of orientations in each voxel for Phi_pred.
%     aompset_per_v = unique(Phi_omp_nz_atoms(:,1));
%     %aompset_per_v
%     % Find the non-zero orientations in voxel vp of Phi_exp
%     Phi_greedy_nz_atoms = find(Phi_greedy(:,vp,:));
%     
%     % Make a set of active orientations for voxel vp
%     agreedyset_per_v = unique(Phi_greedy_nz_atoms(:,1));
%     %agreedyset_per_v
%     
%     diff = setdiff(aompset_per_v, agreedyset_per_v);
%     if ~isempty(diff)
%         %aompset_per_v agreedyset_per_v diff
%     end
% end


%%
% x = 2:5:42;
% GD_angdiff = [0.0079, 0.0077, 0.0076, 0.0076, 0.0076, 0.0076, 0.0075, 0.0075, 0.0075];
% OMP_angdiff = [0.1282, 0.1230, 0.1204, 0.1190, 0.1172, 0.1161, 0.1156, 0.1154, 0.1153];
% 
% bar(x, OMP_angdiff);
% hold on
% bar(x, GD_angdiff);
% 
% xlabel('Number of selected orientations');
% ylabel('Average Angular Differences of Fascicle Directions)');
% legend('OMP', 'Greedy');

%%
dataIndex = 2;
if dataIndex == 1
    X = 0:15;
    Greedy =[480470, 94864, 31480.3, 25510.8, 17358.8, 11090.7, 6692.51, 3976.73, ...
    2539.5, 1860.75, 1612.04, 1463.07, 1381.2, 1288.1, 1214.9, 1131.38];

    OMP = [478410, 478212, 478013, 477815, 477616, 477418, 477220, 477021, 476823, ...
    476625, 476428, 476230, 476032, 475834, 475637, 475439];
elseif dataIndex == 2
    X = 0:15;
    Greedy = [514630, 1951.22, 1872.21, 1761.85, 1708.85, 1590.13, 1538.91, ...
        1419.56, 1376.86, 1262.26, 1218.41, 1118.91, 1095.75, 1012.2, 1002.23, 994.517];
    OMP = [512770, 512554, 512336, 512118, 511900, 511682, 511464, 511247, 511029, ...
    510812, 510594, 510377, 510160, 509943, 509726, 509509];
end
%plot(X, OMP);

semilogy(X, OMP, 'LineWidth',3);
hold on
%plot(X, Greedy);
semilogy(X, Greedy, 'LineWidth',3);
xlabel('Optimization Step');
if dataIndex == 1
    ylabel('Reconstruction error of Y for Arcuate');
elseif dataIndex
    ylabel('Reconstruction Error');
end
legend('OMP', 'Greedy');


%% Plotting the average of precise ang diff, Martha suggested
% dataIndex = 2;
% 
% if dataIndex == 1
%     X = 2:10:62;
%     Greedy = [48.5095, 14.7677, 15.6316, 15.1871, 15.6394, 15.1613, 15.5747];
%     OMP = [68.8883, 178.5572  245.7969  304.7122  349.2896, 373.0730  388.8001];
% elseif dataIndex == 2
%     X = 2:8:52;
%     Greedy = [59.3420, 26.4703, 26.0034, 28.1007, 28.2257, 29.8276, 30.8977];
%     OMP = [69.8499, 157.1341, 218.0036, 266.0400, 310.8500, 343.8009, 362.6211];
% end
% 
% semilogy(X, OMP, 'LineWidth',3);
% hold on
% %plot(X, Greedy);
% semilogy(X, Greedy, 'LineWidth',3);
% xlabel('Iteration Number During Optimization');
% if dataIndex == 1
%     ylabel('Average of Angular Difference over all voxels - Arcuate');
% elseif dataIndex == 2
%     ylabel('Average of Angular Difference over all voxels - ARC-SLF');
% end
% legend('OMP', 'Greedy');



