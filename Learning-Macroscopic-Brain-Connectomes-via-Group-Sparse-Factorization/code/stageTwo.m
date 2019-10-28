function [Phib, message] = stageTwo(Y, B, Phi, voxels, voxel_vicinity, atom_vicinity, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% key parameters explanations: 
% Phi is of the original size and with initialized values from stage 1
% voxels are a list of indexes for selected voxels at stage 1 (e.g. [52, 65, ....])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
% fill in necessary setting parameters   
maxIter = pars.maxIter;
relTolerance = pars.relTolerance;
stepSize = pars.stepSize;
lambdaL1 = pars.lambdaL1;
lambdaGroup = pars.lambdaGroup;
relThres = pars.relThres;
truncate = pars.truncate;
useProximal = pars.useProximal;
%%
% voxel_indices a list of indexes for compact Phi's voxels (e.g. [52,53,65,..])
% of course, if learning full voxels at stage 1, then voxels = voxel_indices
% voxels and voxel_indices could be the same if learning all voxels at
% stage1
if pars.dataIndex == 1
    load('../data/subsets/voxel_indices.mat');
elseif pars.dataIndex == 2
    load('../data/newsubsets/voxel_indices.mat');
end

% set it 1 if updating fascicle screen after each iteration
% this is still screen based on fascicles of current choice
% has nothing to do with neighbor information
% e.g. this is for f_screen rather than f_mask
% check these two screen/mask below
% I do not think this flag parameter is necessary
% optimize_fascicles = 0;


%% Dimensions
% # of total orientations
Na = size(Phi, 1);
% # of total voxels (~190k)
Nv = size(Phi, 2);
% # of total fascicles
Nf = size(Phi, 3);
% # of selected voxels
Nvb = size(Y, 2);

%% Ones vectors
ones_f = ones(Nf, 1);
ones_a = ones(Na, 1);
ones_vb = ones(Nvb, 1);

%% Initial state
% only those voxels dealt at stage 1; the values are initialized from stage 1
% notice Phi,Phi_ls have the size of the original Phi
% Phib and Phib_ls have the size of compact Phi
Phib = Phi(:,voxels,:);
Phi_ls = Phi;
Phib_ls = Phib;

%% Build a grad selection mask (spread atoms to neighbor voxels)
[~, cols, vals] = find(voxel_vicinity(:,voxels));
% Nv x Nvb, notice Nvb is the index within the selected voxels (1 to Nvb)
% while Nv is the index of the original Phi (1 to ~190k)
% Gv is a Nv x Nvb tensor, where Gv(i,j) = 1 denotes
% that voxel i is in the neighborhood of voxel j
% notice i and j are in difference index system!
Gv = sptensor([vals, cols], 1, [Nv, Nvb]);

% like Gv, but for orientations
% Ga is a Na x Na tensor, where Ga(i,j) = 1 denotes
% that orientation i is in the neighborhood of orientation j
% notice unlike Gv, i and j are in the same index system
[~, cols, vals] = find(atom_vicinity);
Ga = sptensor([vals, cols], 1, [Na, Na]);

% this construct a mask matrix that only the position
% that a voxel and an orientation belong to a selected
% voxel's neighbor and orientation's neighbor will be nonzero
% notice ttt is tensor times tensor, which changes the order of dimensions
% that's way using permute to make the result is still 
% oritentations x voxels x fascicles
% note that entry_masks do not care about if the voxels are in the compact
% Phi or not; it only cares about the neighbors of voxels in the compact
% Phi and the neighbors of orientations
% In other words, the nonzero could be in the position of voxels not in the
% compact Phi
entry_masks = ttt(ttt(Phib, Gv, 2, 2), Ga, 1, 2);
entry_masks = permute(entry_masks, [3,2,1]);
% Note that now we select the voxels in the masks, so we are still considering
% only the voxels selected at Stage 1 (among voxels in the compact Phi)
entry_masks = entry_masks(:,voxels,:);
entry_masks = ones(entry_masks);

% build atoms and fascicles filter masks

% a_mask is orientations x voxels
a_mask = ones(ttv(entry_masks, ones_f, 3));

% f_mask is fasciles x voxels
f_mask = ones(ttv(entry_masks, ones_a, 1));
f_mask = permute(f_mask, [2, 1]);

% build orientations and fascicles masks 
% for only possible values
% a_screen and f_screen do not contain neighbors information while a_mask and f_mask does;  

a_screen = ones(ttv(Phib, ones_f, 3));

f_screen = ones(ttv(Phib, ones_a, 1));
f_screen = permute(f_screen, [2,1]);

%% Calculate prediction and diff
Y_hat = ttm(ttv(Phib_ls, ones_f, 3), B, 1);
% if Y_hat is still sparse tensor, then convert to sparse matrix
if (strcmp(class(Y_hat), 'sptensor'))
    Y_hat = spmatrix(Y_hat);
else % 'dense tensor'
    Y_hat = Y_hat.data;
end
Y_diff = Y_hat - Y;

%% Calculate loss

% sign of Phib_ls
sign_Phib_ls = sptensor(Phib_ls.subs, sign(Phib_ls.values), size(Phib_ls));

% absolute value of Phi_ls and Phib_ls
abs_Phi_ls = sptensor(Phi_ls.subs, abs(Phi_ls.values), size(Phi_ls));
abs_Phib_ls = sptensor(Phib_ls.subs, abs(Phib_ls.values), size(Phib_ls));

% this is l1 norm of Phib_ls
% notice this is the same as 
% sum(abs_Phi_ls.values) as nonzero part 
% of Phib_ls and Phi_ls is the same
norm1 = sum(abs_Phib_ls.values);

% this is |\Phi| x_1 G_A^T in the paper
% notice in the paper G_A and G_V are matrices
% but in this implementation, Ga and Gv are both tensors
norm_g1_x = ttt(Ga, abs_Phi_ls, 1, 1);

% this is A^2 in the paper
% namely (|\Phi| x_1 G_A^T)^2 x_2 G_V^T
norm_g1_x2 = ttt(Gv, norm_g1_x .* norm_g1_x, [1], [2]); % Nvb x Na x Nf
% get squre root, namely norm_g1 is A in the paper
% but need to adjust dimension order since Ga and Gv are 
% tensors rather than matrices
if (strcmp(class(norm_g1_x2), 'sptensor'))
    norm_g1 = sptensor(norm_g1_x2.subs, sqrt(norm_g1_x2.vals), size(norm_g1_x2)); % Nvb x Na x Nf
    lpart3 = lambdaGroup * sum(norm_g1.values);% group regularization
else % 'tensor'
    norm_g1 = norm_g1_x2.^(1/2); % Nvb x Na x Nf
    lpart3 = lambdaGroup * sum(norm_g1);
end

% norm_g1 is the group regularization part inside the sum of f, Gv and Ga
% notice that this is also A in the paper's gradient denotion
norm_g1 = permute(norm_g1, [2, 1, 3]); % Na x Nvb x Nf

%lpart1 = norm(Y_diff)^2;%recovery loss
lpart1 = sum(sum(Y_diff .* Y_diff));%recovery loss
lpart2 = lambdaL1 * norm1;% l1 norm
lold = lpart1 + lpart2 + lpart3;

%% Optimization Loop
niter = 1;
%% decide if we are tuning stepsize, so we won't repeat the calculation 
%% of gradient (which is very costly) last time
skip = false;
while (true)
    %% Calculate the gradient
    % gradient for the prediction difference
    % grad_p1_x_t could be quite NOT sparse
    if (skip == false)
        grad_p1_x = 2 * (B)' * Y_diff;
        grad_p1_x_t = sptensor(grad_p1_x);

        % gradient for the L1 regularization
        if (useProximal)
            % if useProximal, we do not need subgradient of l1 regularizer
            g = sptensor([1,1,1], [0], size(sign_Phib_ls));
        else
            grad_l1 = sign_Phib_ls;
            g = lambdaL1 * grad_l1;
        end

        % 1 / A
        % Because norm_g1 could be tense tensor, so cannot use norm_g1.subs or
        % norm_g1.vals
        [subs, vals] = find(norm_g1);
        grad_g1_x2 = sptensor(subs, 1 ./ vals, size(norm_g1));
        size(subs)

        % the whole gradient
        fprintf('forLoop!\n')
        [resSub, resVal] = forLoop(lambdaGroup, a_mask.subs, a_mask.vals, ...
                                   f_mask.subs, f_mask.vals, entry_masks.subs, ...
                                   entry_masks.vals, Ga.subs, Ga.vals, ...
                                   Gv.subs, Gv.vals, grad_p1_x_t.subs, ...
                                   grad_p1_x_t.vals, grad_g1_x2.subs, ...
                                   grad_g1_x2.vals, Phib_ls.subs, Phib_ls.vals, ...
                                   voxels, g.subs, g.vals, a_screen.subs, ...
                                   a_screen.vals, f_screen.subs, ...
                                   f_screen.vals);

        grad = sptensor(resSub, resVal, g.size);
    else
        skip = false;
    end
    
    %% Update using GD
    Phib_tmp = Phib_ls;
    Phib_ls = Phib_ls - stepSize * grad;
    % if using proximal 
    % again, this is not correct 
    % since we have not found a correct way
    % to do so (so we keep use_proximal 0)
    if (useProximal)
        Phib_ls_vals_t = abs(Phib_ls.vals) - stepSize * lambdaL1;
        Phib_ls_vals_t(Phib_ls_vals_t < 0) = 0;
        Phib_ls = sptensor(Phib_ls.subs, sign(Phib_ls.vals).* ...
                 Phib_ls_vals_t, Phib_ls.size);       
    end
    
    % update results
    Phi_ls(:,voxels',:) = Phib_ls;
    
    %% Calculate prediction and diff
    Y_hat = ttm(ttv(Phib_ls, ones_f, 3), B, 1);
    % if Y_hat is still sparse tensor, then convert to sparse matrix
    if (strcmp(class(Y_hat), 'sptensor'))
        Y_hat = spmatrix(Y_hat);
    else % 'dense tensor'
        Y_hat = Y_hat.data;
    end
    Y_diff = Y_hat - Y;

    %% Calculate loss

    % sign of Phib_ls
    sign_Phib_ls = sptensor(Phib_ls.subs, sign(Phib_ls.values), size(Phib_ls));

    % absolute value of Phi_ls and Phib_ls
    abs_Phi_ls = sptensor(Phi_ls.subs, abs(Phi_ls.values), size(Phi_ls));
    abs_Phib_ls = sptensor(Phib_ls.subs, abs(Phib_ls.values), size(Phib_ls));

    % this is l1 norm of Phib_ls
    % notice this is the same as 
    % sum(abs_Phi_ls.values) as nonzero part 
    % of Phib_ls and Phi_ls is the same
    norm1 = sum(abs_Phib_ls.values);

    % this is |\Phi| x_1 G_A^T in the paper
    % notice in the paper G_A and G_V are matrices
    % but in this implementation, Ga and Gv are both tensors
    norm_g1_x = ttt(Ga, abs_Phi_ls, 1, 1);

    % this is A^2 in the paper
    % namely (|\Phi| x_1 G_A^T)^2 x_2 G_V^T
    norm_g1_x2 = ttt(Gv, norm_g1_x .* norm_g1_x, [1], [2]); % Nvb x Na x Nf
    % get squre root, namely norm_g1 is A in the paper
    % but need to adjust dimension order since Ga and Gv are 
    % tensors rather than matrices
    if (strcmp(class(norm_g1_x2), 'sptensor'))
        norm_g1 = sptensor(norm_g1_x2.subs, sqrt(norm_g1_x2.vals), size(norm_g1_x2)); % Nvb x Na x Nf
        lpart3 = lambdaGroup * sum(norm_g1.values);% group regularization
    else % 'tensor'
        norm_g1 = norm_g1_x2.^(1/2); % Nvb x Na x Nf
        lpart3 = lambdaGroup * sum(norm_g1);
    end

    % norm_g1 is the group regularization part inside the sum of f, Gv and Ga
    % notice that this is also A in the paper's gradient denotion
    norm_g1 = permute(norm_g1, [2, 1, 3]); % Na x Nvb x Nf

    %lpart1 = norm(Y_diff)^2;%recovery loss
    lpart1 = sum(sum(Y_diff .* Y_diff));%recovery loss
    lpart2 = lambdaL1 * norm1;% l1 norm
    lnew = lpart1 + lpart2 + lpart3;
    
    %% current results
    fprintf(1, '%g: old obj = %g, total obj = %g, reconstruction error = %g, l1 error = %g, group error = %g, stepsize = %g\n', niter, lold, lnew, lpart1, lpart2, lpart3, stepSize);
    
    %% termination condition
    if ((lold - lnew < lold * relTolerance) || (niter >= maxIter))
        message = ['Done optimization with lold = ' num2str(lold) ' and lnew = ' num2str(lnew) ', stepsize = ' num2str(stepSize)];
        fprintf(1, [message '\n']);    
        if (lold < lnew)
            Phib_ls = Phib_tmp;
            stepSize = stepSize / 10;
            skip = true;
        else
            break;
        end
    end

    %% Update screen for fascicles
    % If zeroed fascicles are free to be updated nonzero, 
    % or if nonzeroed fascicles are optimized to be nonzero during the
    % optimization,
    % then update screen of fascicle after each updation since it will change
    %if optimize_fascicles 
        %f_screen = ones(ttv(Phib_ls, ones_a, 1));
        %f_screen = permute(f_screen, [2,1]);
    %end
    %% Update loss value and iteration for next run
    if (lold >= lnew)
        lold = lnew;
        niter = niter + 1;
    end
end
%% truncate or not
%% set values directly 0 if their absolute values are tiny enough
if truncate
    Phi_sum_a = ttv(abs_Phib_ls, ones_a, 1);
    if (strcmp(class(Phi_sum_a), 'sptensor'))
        Phi_sum_a = spmatrix(Phi_sum_a);
    else % 'tensor'
        Phi_sum_a = Phi_sum_a.data; % Nvb x Na x Nf
    end

    v_max = max(Phi_sum_a, [], 2);
    [subs, vals] = find(abs_Phib_ls);
    % if absolute values are smaller than some 
    % portion of the largest ones 
    % (for a voxel, the largest sum of all orientations on one fascicle)
    % then set them 0
    flags = vals < relThres * v_max(subs(:,2));
    [subs, vals] = find(Phib_ls);
    vals(flags) = 0;
    Phib = sptensor(subs, vals, size(Phib));
else
    Phib = Phib_ls;
end
message = [message ' proximal ' num2str(useProximal)];
end