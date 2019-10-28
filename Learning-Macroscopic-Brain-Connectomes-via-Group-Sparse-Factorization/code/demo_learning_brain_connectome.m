function demo_learning_brain_connectome(dataIndex)

addpath(genpath('tensor_toolbox'),genpath('../data'));	

%% Extracting sub tensors and matrices
% Extrcts the dataset into subsets and newsubsets. The resulting matrices
% are used later on in the code. Their data has beem computed and extracted
% from the original dataset.
if dataIndex == 1
    save_compact_matrices(1);
elseif dataIndex == 2
    save_compact_matrices(2);
else
    fprintf('There is no dataset with dataIndex %d', dataIndex);
end

%% Setup parameters
% Setup size of orientation candidate set: k
k = 5;
% setup whether we want stage1 of the algorithm to initialize Phi with
% expert Phi (stg1 = 1), with OMP (stg1 = 2), or with OrientationGreedy(stg1 = 3).
stg1 = 3;


%% Start optimisation procedure
if dataIndex == 1
    %Nsv = 11823;
    Nsv = 5;
    brain(dataIndex, Nsv, k, stg1);
elseif dataIndex == 2
    Nsv = 15033;
    brain(dataIndex, Nsv, k, stg1);
else
    fprintf('There is no dataset with dataIndex %d', dataIndex);
end


%% Visualization of the optimized Phi
if dataIndex == 1
    load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')
    load();
    load();
    
elseif dataIndex == 2
    load('../data/fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')
else
    fprintf('There is no dataset with dataIndex %d', dataIndex);
end