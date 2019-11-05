%% First part, It compute the precise angular difference for  based on the best solution of orientations' weights.
% dataIndex = 2;
% opt = false;
% 
% stage1 = 3;
% numOrient = 5;
% regStage1 = 0.01;
% 
% %% Load Data
% if dataIndex == 1
%     Nv = 11823;
%     load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
%     orients = fe.life.M.Atoms.orient;
%     load('../data/subsets/B.mat');
%     load('../data/subsets/weights.mat');
%     load('../data/subsets/compact_Y.mat');
%     load('../data/subsets/compact_Phi_withw.mat');
%     load('../data/stage1/Phi_11823_1_1_5_0.01.mat');
%     Phi_gt = Phi_sp;
%     if opt
%         if stage1 == 3
%             load('../experiments/dataset_one_new/Phi_11823_1_3_5_0.01_10_10_gs_0.01.mat');
%         elseif stage1 == 2
%             load('../experiments/dataset_one_new/Phi_11823_1_2_5_0.01_10_10.mat');
% 
%         end
%         Phi_pred = Phi_sp;
%     else
%         if stage1 == 3
%             load('../data/stage1/Phi_11823_1_3_5_0.01.mat');
%         elseif stage1 == 2
%             load('../data/stage1/Phi_11823_1_2_5_0.01.mat');
%         end
%         Phi_pred = Phi_sp;
%     end
%     saveName = ['../data/stage1/Angdiff_' int2str(Nv) '_' int2str(dataIndex) '_' int2str(stage1) '_' int2str(numOrient),'.mat'];
% elseif dataIndex == 2
%     Nv = 15033;
%     load('fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat');
%     orients = fe.life.M.Atoms.orient;
%     load('../data/newsubsets/B.mat');
%     load('../data/newsubsets/weights.mat');
%     load('../data/newsubsets/compact_Y.mat');
%     load('../data/newsubsets/compact_Phi_withw.mat');
%     load('../data/stage1/Phi_15033_2_1_5_0.01.mat');
%     Phi_gt = Phi_sp;
%     if opt
%         if stage1 == 3
%             load('../experiments/dataset_two_new/Phi_15033_2_3_5_0.01_10_10 _gs_0.01_15extra.mat');
%         elseif stage1 == 2
%             load('../experiments/dataset_two_new/Phi_15033_2_2_5_0.01_10_10.mat');
%         end
%         Phi_pred = Phi_sp;
%     else
%         if stage1 == 3
%             load('../data/stage1/Phi_15033_2_3_5_0.01.mat');
%         elseif stage1 == 2
%             load('../data/stage1/Phi_15033_2_2_5_0.01.mat');
%         end
%         Phi_pred = Phi_sp;
%     end
%     saveName = ['../data/stage1/Angdiff' int2str(Nv) '_' int2str(dataIndex) '_' int2str(stage1) '_' int2str(numOrient),'.mat'];
% end
% 
% % load w (linear parameters)
% % since Phi absorbs w
% % hence, w is a vector of ones 
% % (only used to sum up the dimension)
% w = ones(size(w));
% Nf = size(w, 1);
% Na = size(B, 2);

% 
% %% Algorithm starts
% % sum over all fascicles in ground truth Phi_f_gt(Na * Nv)
% Phi_f_gt = ttv(Phi_gt, w, 3);
% % sum over all fascicles in ground truth Phi_f_pred(Na * Nv)
% Phi_f_pred = ttv(Phi_pred, w, 3);
% 
% % angular distance matrix (Na * Nv), contains the 
% ang_diff = sptensor([], [], size(Phi_f_gt));
% for v_i = 1:Nv
%     if mod(v_i, 100) == 1
%             fprintf(1, '*');
%     end
%     pred_a_ind = unique(find(Phi_f_pred(:, v_i)));
%     Bsmall = B(:,pred_a_ind);
%     
%     gt_a_ind = unique(find(Phi_f_gt(:, v_i)));
%     for a_i = 1:size(gt_a_ind)
%         % Compute the signal y of ground truth for each node.
%         y_a = Phi_f_gt(gt_a_ind(a_i), v_i) * B(:, gt_a_ind(a_i));
%         vals = pinv(Bsmall' * Bsmall) * Bsmall' * y_a;
%         
%         exp_a_vec = orients(:, gt_a_ind(a_i));
%         pred_a_vec = orients(:, pred_a_ind) * vals;
%         
%         angdif = AngularDistance(exp_a_vec, pred_a_vec);
%         
%         ang_diff(gt_a_ind(a_i), v_i) = angdif;
%         
%     end
% end
% 
% 
% save(saveName, 'ang_diff');
% 
% %% Angular distance
% % This fuction returns the angular distance of any two input vectors.
%     function angdist = AngularDistance(vec1, vec2)
%         angdist = abs(atan2d(norm(cross(vec1, vec2, 1)), dot(vec1, vec2, 1)));
%         if angdist > 90
%             angdist = 180 - angdist;
%         end
%     end
% 

%% After you ran the first part, uncomment this part and comment last part to get the graph
dataIndex = 2;
if dataIndex == 1
    datasetName = 'Arcuate';
    saveName1 = '../data/stage1/Angdiff11823_1_2_5.mat';
    saveName2 = '../data/stage1/Angdiff11823_1_3_5.mat';
elseif dataIndex == 2
    datasetName = 'ARC-SLF';
    saveName1 = '../data/stage1/Angdiff_15033_2_2_5.mat';
    saveName2 = '../data/stage1/Angdiff_15033_2_3_5.mat';
end
    
load(saveName1)
ang_diff.vals
ang_diff.vals(imag(ang_diff.vals)~=0)
%cdfplot(real(ang_diff.vals));
histogram(real(ang_diff.vals),'Normalization','probability','DisplayStyle','stairs', 'LineWidth',1);
hold on

load(saveName2)
ang_diff.vals
ang_diff.vals(imag(ang_diff.vals)~=0)
%cdfplot(real(ang_diff.vals));
histogram(real(ang_diff.vals),'Normalization','probability','DisplayStyle','stairs', 'LineWidth',1);
xlabel(['Angular Difference of Orientations - ', datasetName]);
ylabel('P(Angular Differences of Fascicle Orientations)');
legend('OMP', 'Greedy');
