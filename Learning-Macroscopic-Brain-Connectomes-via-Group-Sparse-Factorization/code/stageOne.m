function [Phi_sp, message] = stageOne(Y, B, Phi_sp, PhiP, w, voxels, voxel_vicinity, atom_vicinity, Nsv, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% key parameters explanations: 
% Phi is an empty compact tensor, PhiP is the original compact tensor (absorting w already)
% voxels is voxel_indices, namely the indexes of necessary voxels
% Nsv is the maximum # of voxels to be processed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global Parameters
% fill in necessary setting parameters  
dataIndex = pars.dataIndex;
calStage1 = pars.calStage1;
stage1 = pars.stage1;
stage2 = pars.stage2;
numOrient = pars.numOrient;
regStage1 = pars.regStage1;
saveNameStage1 = pars.saveNameStage1;
saveNamePhi1 = pars.saveNamePhi1;
saveNamePhi2 = pars.saveNamePhi2;

message = [];
%% Local variables

% number of directions of atoms (namely orientations)
Na = size(B, 2);

% number of fasicles
Nf = size(w, 1);

% number of all voxels
Nv = size(voxel_vicinity, 2);

% start from the first voxels in compact_phi
% Note that vlist contains the index in the original phi
vlist = voxels(1);

% a logical array to denote if the voxel is added already
% notice size(voxel_vicinity,2) is the # of original voxels (~190k)
added = false(1, size(voxel_vicinity, 2));
added(vlist) = true;% now the first voxel is set true since we add it already

% transform voxels to a logical array to denote if the index of the
% original voxels is in the compact voxels
a = zeros(1, size(voxel_vicinity, 2));
a(voxels) = true;
voxels = logical(a);

% for the index in the original 190k voxels Phi, v_compact_ind will return
% the index in the compact voxels
v_compact_ind = zeros(size(voxels));
v_compact_ind(voxels) = 1 : nnz(voxels);

% notice that v_compact_ind and voxel_indices are mutually inverse encoders / decoders
% v_compact_ind maps index of the original voxels to index of the compact voxels
% voxel_indices maps index of the compact voxels to index of the original voxels


%% STAGE 1
%% At Stage 1, our goal is to initialize Phi;
% At the beginning, we will add the first voxel in the compact Phi in the candidate list;
% We will deal with the voxels in the candidate list sequentially while 
% adding the neighbor voxels of the currently iterated voxel to the candidate list
% if they are not in the list yet;
% We will initialize the selected fasicles (nonzeros across the whole orientations) and
% orientations of the slice matrix of the currently iterated voxel (different strategies to do so);
% for the other parts, just keep the same values as the expert phi

% print out what settings we are trying to produce or load
fprintf("Nsv : %d; dataIndex : %d; stage1 : %d; numOrient : %d; regStage1 : %f\n", Nsv, dataIndex, stage1, numOrient, regStage1);

if ~calStage1
    %% If not initialize Phi through stage1
    % we will load the Phi returned by stage1 last time
    try
        load(saveNameStage1);
    catch
        fprintf('Error: Even No Stage1 Results!');
        return;
    end
    if pars.dataIndex == 1
        % it is possible there are not results from previous experiments
        % then we can just use the result from stage1
        try
            load(saveNamePhi1);
        catch
            fprintf('Warn: Use Results from Stage 1 Directly');
        end
    elseif pars.dataIndex == 2
        try
            load(saveNamePhi2);
        catch
            fprintf('Warn: Use Results from Stage 1 Directly');
        end
    end
else
    %% Otherwise, let's process Stage1
    % now iterate the voxels
    fprintf('stage1!\n');
    for v_i = 1 : Nsv
        % get enough candidate voxels
        if (v_i > size(vlist, 1)) % No new neighbor voxels could be added in
            Nsv = size(vlist, 1);
            fprintf("This is new Nsv: %d\n", Nsv);
            break;
        end
        % print * for each 100 voxels finished
        % makes the screen not so boring 
        % as well as we can expect the time to finish
        % the job
        if mod(v_i, 100) == 1
            fprintf(1, '*');
        end
        
        v0 = vlist(v_i);% v0 is the index of vi in the original voxels
        v = v_compact_ind(v0);% v is the index of vi in the compact voxels

        y = Y(:,v);

        % find adjacent voxels of v0
        ad_voxels = voxel_vicinity(:,v0);
        % ad_voxels(1) is v0 itself, so we start from 2
        neig_voxels = ad_voxels(2:length(ad_voxels));
        % set the first element false, because that
        % denotes the current voxel
        % only among the 26 neighbors who are in compact voxel sets 
        % and not added yet are candidate voxels
        % next_voxels are the index for original voxels
        next_voxels = ad_voxels([false, voxels(neig_voxels) & ~added(neig_voxels)]);
        vlist = [vlist;next_voxels];
        added(next_voxels) = true;
        
        % For Model 1 (just use the expert Phi as initialization for Stage2)
        if stage1 == 1
            % Only initialize Phi by expert Phi; 
            % so nothing to do yet
        else
            % Note ones here is only makes the nonzero elements in PhiP be 1
            % Not like the usual ones
            % But do we really need ones?? (TODO)
            P = spmatrix(ones(PhiP(:,v,:)));
            % The reason why we need P, it is because sum cannot be used
            % directly for sptensor. So we use a sparse matrix P to store the
            % sptensor to be used by sum
            selection_f = sum(P, 1) ~= 0;% keep the fasicles in the expert Phi       
            selection_f_ind = find(selection_f);% index the kept fasicle
            
            % For OMP for the selection of orientations
            if stage1 == 2
                %% Farzane's implementation
                Phiv = Phi_sp(:,v,:);
                x = full(spmatrix(Phiv(:,:)));
                selection_a = false(1, Na); % at the beginning, no directions have been selected
                
                % get difference between target vector and the combination
                % of basis now we have
                r = y - B * x * w;
                
                while norm(r) > 1e-12 && nnz(selection_a) < numOrient
                    angle_cos = (r' * B)./(norm(r) * sqrt(sum(B .^ 2)));

                    permute = 1 : Na;  
                    B_sort = [permute;angle_cos;B];
                    B_sort = B_sort(:, ~selection_a);
                    B_sort = B_sort';
                    % sort according to the second column(angle)
                    B_sort = sortrows(B_sort, -2);
                    B_sort = B_sort';
                    
                    if (B_sort(2, 1) < 0.01) % the corelation is too small
                        Nsv = v_i - 1;
                        fprintf("This is new Nsv: %d, because reaching ending condition2 of OMP\n", Nsv);
                        break;
                    end
    
                    % update selections
                    % select top ``numOrient'' directions
                    selection_a(B_sort(1,1)) = true;
                    selection_a_ind = find(selection_a);
                    Nas = nnz(selection_a);
                    
                    Nfs = nnz(selection_f);

                    Bsmall = B(:,selection_a_ind);
                    cinit = pinv(Bsmall' * Bsmall + regStage1 * eye(size(Bsmall, 2))) * Bsmall' * y;
                    xsmall = repmat(cinit, 1, Nfs);
                    
                    x(selection_a_ind, selection_f_ind) = xsmall;
                    %cinit = Bsmall \ y;
                    %cinit = cinit/Nfs;
                    
                    r = y - Bsmall * x(selection_a_ind, :) * w;
                    
                    if norm(r) <= 10e-12
                        fprintf('norm(r) is: %d, reaching the ending condition of OPM, no more projection required!\n num of atom selected is: %d', norm(r), nnz(selection_a))
                        break
                    end
                    
                end
                
                selection_a_ind = find(selection_a);

                Nfs = nnz(selection_f);
                %fprintf('selection_a in 2 %d \n', selection_a_ind)

                Bsmall = B(:,selection_a_ind);
                cinit = pinv(Bsmall' * Bsmall + regStage1 * eye(size(Bsmall, 2))) * Bsmall' * y;
                %cinit = Bsmall \ y;
                %cinit = cinit/Nfs;

                Phib = sptensor(repmat(cinit, 1, Nfs));
                if (Nfs == 1)
                    Phi_sp(selection_a_ind,v,selection_f_ind) = Phib(:,1);
                else
                    Phi_sp(selection_a_ind,v,selection_f_ind) = Phib;
                end
                               
                %% Lei's implementation
%                 Phiv = Phi_sp(:,v,:);
%                 x = full(spmatrix(Phiv(:,:)));
%                 selection_a = false(1, Na); % at the beginning, no directions have been selected
%                 
%                 % get difference between target vector and the combination
%                 % of basis now we have
%                 r = y - B * x * w;
%                 
%                 % ending condition
%                 if (norm(r) / norm(y) < 0.1) % the magnitude of diff becomes small
%                     Nsv = v_i - 1;
%                     fprintf("This is new Nsv: %d, because reaching ending condition1 of OMP\n", Nsv);
%                     break;
%                 end
% 
%                 angle_cos = (r' * B)./(norm(r) * sqrt(sum(B .^ 2)));
%                 size(angle_cos)
% 
%                 permute = 1 : Na;
% 
%                 B_sort = [permute;angle_cos;B];
%                 size(B_sort)
%                 B_sort = B_sort(:, ~selection_a);
%                 B_sort = B_sort';
%                 % sort according to the second column(angle)
%                 B_sort = sortrows(B_sort, -2);
%                 B_sort = B_sort';
%                 permute = B_sort(1,:);
%                 if (B_sort(2, 1) < 0.01) % the corelation is too small
%                     Nsv = v_i - 1;
%                     fprintf("This is new Nsv: %d, because reaching ending condition2 of OMP\n", Nsv);
%                     break;
%                 end
% 
%                 % update selections
%                 % select top ``numOrient'' directions
%                 selection_a(permute(1 : numOrient)) = true;
%                 selection_a_ind = find(selection_a);
% 
%                 Nfs = nnz(selection_f);
%                 %fprintf('selection_a in 2 %d \n', selection_a_ind)
% 
%                 Bsmall = B(:,selection_a_ind);
%                 cinit = pinv(Bsmall' * Bsmall + regStage1 * eye(size(Bsmall, 2))) * Bsmall' * y;
%                 %cinit = Bsmall \ y;
%                 %cinit = cinit/Nfs;
% 
%                 Phib = sptensor(repmat(cinit, 1, Nfs));
%                 if (Nfs == 1)
%                     Phi_sp(selection_a_ind,v,selection_f_ind) = Phib(:,1);
%                 else
%                     Phi_sp(selection_a_ind,v,selection_f_ind) = Phib;
%                 end
                
            % GreedyDirections (main method in the paper)    
            elseif stage1 == 3
                BTB = B' * B;     
                b = B' * y;              
                C = diag(BTB);       
                myG = (b .^ 2) ./ C;       
                [~, aMax] = max(myG);             
                S = aMax;        
                cInv = 1 / C(aMax); 
                for i = 2 : numOrient
                    [gS, V] = ComputeGain(BTB, C, S, cInv, b);
                    gBar =  0.01 * gS + myG;
                    gBar(S) = -1e20;% do not select added ones
                    [~, aMax] = max(gBar);
                    cInv = [cInv + V(aMax) * cInv * BTB(S, aMax) * BTB(S, aMax)' * cInv, -V(aMax) * cInv * BTB(S, aMax); -V(aMax) * BTB(S, aMax)' * cInv, V(aMax)];
                    S = [S, aMax];
                end

                my_temp = 1 : Na;
                [selection_a, ~] = ismember(my_temp, S);

                % Nas is the # of selected orientations
                % Nfs is the # of selected fascicles
                Nas = nnz(selection_a);
                Nfs = nnz(selection_f);
                
                % this is the first way
                %if (Nas == Nfs == 1)
                %    Phi_sp(S,v,selection_f_ind) = 0.0001 * rand();
                    
                %    % Added for debugging
                %   %Phi_pred_nz_atoms = find(Phi_sp(:,v,:));
                %    %apredset_per_v = unique(Phi_pred_nz_atoms(:,1));
                %    %fprintf('1 Nas is %d and Nfs is %d', Nas, Nfs)
                %    %apredset_per_v
                %elseif (Nas == 1)
                %    randt = 0.0001 * sptenrand(Nfs, Nfs);
                %    Phi_sp(S,v,selection_f_ind) = randt;
                    
                %    % Added for debugging
                %    %Phi_pred_nz_atoms = find(Phi_sp(:,v,:));
                %    %apredset_per_v = unique(Phi_pred_nz_atoms(:,1));
                %    %fprintf('2 Nas is %d and Nfs is %d', Nas, Nfs)
                %    %apredset_per_v
                %elseif (Nfs == 1)
                %    randt = 0.0001 * sptenrand(Nas, Nas);
                %    Phi_sp(S,v,selection_f_ind) = randt;
                    
                %    % Added for debugging
                %    %Phi_pred_nz_atoms = find(Phi_sp(:,v,:));                    
                %    %apredset_per_v = unique(Phi_pred_nz_atoms(:,1));
                %    %fprintf('3: in voxel %d Nas is %d and Nfs is %d', v, Nas, Nfs)
                %    %S
                %    %randt
                %    %Phi_pred_nz_atoms
                %else
                %    %Phi(S,v,selection_f_ind) = sptensor(0.0001 * rand(Nas, Nfs));
                %    Bsmall = B(:,S);
                %    % closed form solution of ||D\phi - y||_2^2 + regStage1 * ||\phi||_2^2
                %    cinit = pinv(Bsmall' * Bsmall + regStage1 * eye(size(Bsmall,2))) * Bsmall' * y;
                %    % since there can be several fascicles, so 
                %    % assign evently to each of them
                %    cinit = cinit / Nfs;
                %    Phi_sp(S,v,selection_f_ind) = sptensor(repmat(cinit, 1, Nfs));
                    
                %    % Added for debugging
                %    %Phi_pred_nz_atoms = find(Phi_sp(:,v,:));
                %    %apredset_per_v = unique(Phi_pred_nz_atoms(:,1));
                %    %fprintf('4 Nas is %d and Nfs is %d', Nas, Nfs)
                %    %apredset_per_v
                %end
                

                % this is the second way
                Bsmall = B(:,S);
                cinit = pinv(Bsmall' * Bsmall + regStage1 * eye(size(Bsmall, 2))) * Bsmall' * y;
                %cinit = cinit / Nfs;
                mytemp = sptensor(repmat(cinit, 1, Nfs));       
                if (Nfs == 1)
                    Phi_sp(S,v,selection_f_ind) = sptensor(mytemp(:,1));
                else
                    Phi_sp(S,v,selection_f_ind) = sptensor(mytemp);
                end    
            end
        end
    end
    
    % vlist the index of voxels processed
    % if Nsv is full voxels, then vlist is
    % the same as voxel_indices
    vlist = vlist(1 : Nsv);
    vlist = sort(vlist);
    % vind : voxels with compact indices
    % e.g. 1st in compact Phi, 2nd in compact Phi,.etc. 
    vind = v_compact_ind(vlist);
    if stage1 == 1
        Phi_sp = PhiP;
        Phib = Phi_sp(:,vind,:);
        Phi_sp(:,:,:) = 0;
        Phi_sp(:,vind,:) = Phib;
    end
    % save Phi of stage 1
    save(saveNameStage1, 'Phi_sp', 'vlist', 'vind');
end

fprintf('\n');
% if NOT go through stage 2
if stage2 == 0
    return;
end

%% STAGE 2
%% call the Stage2 function after getting ready
fprintf('stage2!\n');


% Just to run the experiment on pruning section
% Phi_sp_1 = PhiP;
% Phib_1 = Phi_sp_1(:,vind,:);
% Phi_sp_1(:,:,:) = 0;
% Phi_sp_1(:,vind,:) = Phib_1;
% Phi_sp(:,vind,:) = Phi_sp(:,vind,:) + Phi_sp_1(:,vind,:);
% save([saveNameStage1, 'OMP_GroundTruth'], 'Phi_sp', 'vind', 'vlist');

Phi_all = sptensor([], [], [Na, Nv, Nf]);
Phi_all(:,vlist',:) = Phi_sp(:,vind,:);
% Phi_all is the Phi of original size 

% with only positions of compact voxels returned
% at the first stage nonzero
% (the values are the same as the returned Phi at the first stage)

% voxels_ind : sorted compact indices
% voxels_ind = find(voxels);
% we need voxels_ind here; in other words, we cannot use voxels directly 
% because sptensor type does not support logic objective to index elements
% directly (weird though)

% PhiP_sparse = sptensor([], [], [Na, Nv, Nf]);
% PhiP_sparse(:,voxels_ind,:) = PhiP(:,:,:);
% PhiP_sparse is the original Phi with only compact position of voxels
% nonzero(the values are the same as expert Phi of compact_voxel)
% call stage2 function

[Phib, message] = stageTwo(Y(:,vind), B, Phi_all, vlist, voxel_vicinity, atom_vicinity, pars);

Phi_sp(:,vind,:) = Phib;

end