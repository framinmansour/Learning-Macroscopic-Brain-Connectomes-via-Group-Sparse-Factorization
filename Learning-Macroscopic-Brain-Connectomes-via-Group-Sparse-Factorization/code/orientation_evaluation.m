%% This file is based on orientation evaluation script Farzaen suggested. 
% Note: What Martha suggested was better and more precise. 
% Farzane: I guess we don't need this for ang_diff part.
function [wangdist, angdist, countdist] = orientation_evaluation(Phi_pred, Phi_exp, orients)

% Of course, we need to add tensor toolbox
addpath tensor_toolbox

%%
% This fuction returns the angular distance of any two input vectors.
    function angdist = AngularDistance(vec1, vec2)
        angdist = abs(atan2d(norm(cross(vec1, vec2, 1)), dot(vec1, vec2, 1)));
        if angdist > 90
            angdist = 180 - angdist;
        end
    end

%%
% Returns the powerset of set S
%
% S is a cell array
% P is a cell array of cell arrays
    function [ P ] = PowerSet( S )
        n = numel(S);
        x = 1:n;
        P = cell(1,2^n);
        p_ix = 2;
        for nn = 1:n
            a = combnk(x,nn);
            for j=1:size(a,1)
                P{p_ix} = S(a(j,:));
                p_ix = p_ix + 1;
            end
        end
    end


%%
Na_exp = size(Phi_exp, 1);
Na_pred = size(Phi_pred, 1);
Nv = size(Phi_pred, 2);
Nf_exp = size(Phi_exp, 3);
Nf_pred = size(Phi_pred, 3);

% all of the atoms activated in Phi_exp and Phi_pred
a_exp = unique(Phi_exp.subs(:,1));
a_pred = unique(Phi_pred.subs(:, 1));

ones_f = ones(Nf_exp, 1);

%%
wangdist = 0.0;
angdist = 0.0;
countdist = 0.0;
Total_num_orients_in_Phiexp = 0.0;
for vp = 1:Nv
    if mod(vp,100) == 0
        fprintf("Voxel %d of %d\n", vp, Nv);
    end
    % This variable contains the total distance of orientations per voxel. 
    % We will add the distance of each orientation to it to get the 
    % distance per voxel.
    wdist_per_voxel = 0.0;
    dist_per_voxel = 0.0;
    
    % Find non-zero oriantations in each voxel
    Phi_pred_nz_atoms = find(Phi_pred(:,vp,:));
    
    % Make an orientation set for each voxel. The goal is to find different
    % types of orientations in each voxel for Phi_pred.
    apredset_per_v = unique(Phi_pred_nz_atoms(:,1));
    %fprintf('\n atoms in voxel %d of pred', vp);
    %apredset_per_v
    
    % Compute the powerset of all orientations in each voxel of Phi_pred.
    % We want to find the closeness of each orientation in the
    % corresponding voxel of Phi_exp with each of these subsets in the
    % powerset. We want to measure how far the estimation of Phi_pred would
    % be from the actual orientations in Phi_exp. P_pred is a cell array.
    % Access to elements like P_pred{n}.
    P_pred = PowerSet(apredset_per_v);
    
    % Find the non-zero orientations in voxel vp of Phi_exp
    Phi_exp_nz_atoms = find(Phi_exp(:,vp,:));
    
    % Make a set of active orientations for voxel vp
    aexpset_per_v = unique(Phi_exp_nz_atoms(:,1));
    %fprintf('\n atoms in voxel %d of exp', vp);
    %aexpset_per_v
    Total_num_orients_in_Phiexp = Total_num_orients_in_Phiexp + size(aexpset_per_v);
    
    % The intersection of orientations in Phi_exp and Phi_pred for the same
    % voxel shows that how close could we potentially get to the ground 
    % truth.
    mutual_orients = intersect(apredset_per_v, aexpset_per_v);
    
    % We want the distance, so we need the numeber of orientations which
    % has not been included in Phi_pre(:,vp,:) but included in Phi_exp(:,vp,:)
    % The ideal case is that the length of mutual orientations be exactly
    % the same as the activated orientations in Phi_exp(:,vp,:)
    countdist = countdist + (length(aexpset_per_v) - length(mutual_orients));
    
    % For each unique orientation in vp
    for aexp = 1: length(aexpset_per_v)
        % Find the orientation vector of it from fe_struct. 
        orient_vec_exp = orients(:, aexpset_per_v(aexp));
        
        % For non_weighted expert vector
        vec_exp = orient_vec_exp;
        % Compute the sum of weights over all fascicles which have the same
        % orientation in the same voxel. It shows that what portion of the
        % signal from this voxel is generated from this specific
        % orientation.
        orient_weights_exp = ttv(Phi_exp(aexpset_per_v(aexp),vp,:), ones_f,1);
        
        % Generate the expert's weighted vector of orientations in voxel vp
        wvec_exp = orient_vec_exp * orient_weights_exp;
        
        % It keeps the minimum distance of an orientation in Phi_exp with 
        % a weighted combination of subsets in the powerset of orientations 
        % in Phi_pred.
        min_wdist = Inf;
        min_dist = Inf;
        % For each subset in the set of powerset
        for s = P_pred
            % Excluding none
            if ~isempty(s{:})
                % Which dimension is the one that we want to get the sum of
                % weight over it
                [~,dim] = max(size(Phi_pred(s{:},vp,:)));
                % Fetch the orientation vectors (x,y,z)
                orient_vec_pred = orients(:, s{:});
                
                % To check which subset, use:
                % s{:}
                
                % Sum of weights over all fascicles having the same
                % orientation in a voxel.
                orient_weights_pred = ttv(Phi_pred(s{:},vp,:), ones_f,dim);
                
                % When we only have one orientation activated for getting
                % the vector sum of orientations orient_vec_pred.
                ones_a = ones(1, 1);
                
                % Just to match the sizes with orient_vec_pred
                sz = size(orient_weights_pred,1);
                if sz > 1
                    orient_weights_pred = reshape(orient_weights_pred,[sz,1]);
                    % To get the vector sum over all active orientations
                    % for non-weighted angular distance measurement.
                    ones_a = ones(sz, 1);
                end
                
                % For non-weighted predicted vector, vector sum over all
                % non-weighted orentations in the subset s{:}
                vec_pred = ttv(sptensor(orient_vec_pred), ones_a, 2);
                
                %ttv(sptensor(orient_vec_pred), ttv(Phi_pred(s{:},vp,:), ones_f,dim), 2)
                % Multiply wights by orientations
                wvec_pred = ttv(sptensor(orient_vec_pred), orient_weights_pred, 2);
                if (strcmp(class(wvec_pred), 'sptensor'))
                    wvec_pred = double(wvec_pred);
                else % 'dense tensor'
                    wvec_pred = wvec_pred.data;
                end
                
                % Find the weighted angular distance
                wdistance = AngularDistance(wvec_pred, wvec_exp);
                
                % Find the non-weighted angular distance
                if (strcmp(class(vec_pred), 'sptensor'))
                    vec_pred = double(vec_pred);
                else % 'dense tensor'
                    vec_pred = vec_pred.data;
                end
                distance = AngularDistance(vec_pred, vec_exp);
                
                if min_wdist > wdistance
                    min_wdist = wdistance;
                end
                if min_dist > distance
                    orient_set = s{:};
                    min_dist = distance;
                end
            end
        end
        %fprintf("Voxel %d with min %d\n", vp, min_dist);
        %orient_set
        %s{:}
        wdist_per_voxel = wdist_per_voxel + min_wdist;
        dist_per_voxel = dist_per_voxel + min_dist;
    end
    wangdist = wangdist + wdist_per_voxel;
    angdist = angdist + dist_per_voxel;

end
fprintf('Total number of active ground truth orientations %d \n', Total_num_orients_in_Phiexp);
end
