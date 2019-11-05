function hs = visualize_Phi_with_fiber(fe, Phi, Phi0, pars)
% Phi should be a tensor
% n number of frontal slices to be visualized
load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')    
load('../data/subsets/voxel_indices.mat')
    
    roi = fe.roi;
    orient = fe.life.M.Atoms.orient';
        
    Na = size(Phi,1);
    Nv = size(Phi,2);
    Nf = size(Phi,3);

    ones_f = ones(Nf,1);
    
    if ( ~exist('pars') )
        n = 10;
        fidx = 1:n;
        fig = figure;
        ax = axes('Parent', fig);
    else
        if (~isfield(pars,'n'))
            n = 10;
        else
            n = pars.n;
        end
        
        if (~isfield(pars,'f'))
            fidx = 1:n;
        else
            fidx = pars.f;
        end

        if (~isfield(pars,'axes'))
            fig = figure;
            ax = axes('Parent', fig);
        else
            ax = pars.axes;
            fig = get(ax,'Parent');
            cla(ax);
        end
    end
    
    tsubs.type = '()';
    
    % hold on;
    hold(ax, 'on');
    daspect(ax,[1 1 1])
    grid(ax,'on')
    
    [subs, vals] = find(Phi);
    
    ones_a = ones(size(Phi,1),1);
    if (Nf ~= 1)
        Phi_abs = sptensor(Phi.subs, abs(Phi.vals), size(Phi));
        % lei's change: without summary of atoms. Instead, 
        % use the max of atoms
        f_max = [];
        for i = 1:size(Phi_abs,3)
            f_max =[f_max;max(max(spmatrix(Phi_abs(:,:,i))))];
        end
        %vsum = spmatrix(ttv(Phi_abs, ones_a, 1));
%         v_max = max(vsum')'; %max among fascicles
        %f_max = max(vsum,[],1)';
        %v_max = max(vsum,[],2);
    else
        Phi_abs = abs(Phi);
        %vsum = ones_a' * spmatrix(Phi);
        f_max = Phi_abs';
        v_max = max(spmatrix(Phi_abs(:,:,1)),[],2);
        %v_max = vsum';
    end
    
    [subs, vals] = find(Phi);
%     vals = vals ./ v_max(subs(:,2));

    vals = vals ./ f_max(subs(:,3));
    Phi = sptensor(subs, vals, size(Phi));


    [subs, vals] = find(Phi0);
    
    ones_a = ones(size(Phi0,1),1);
    if (Nf ~= 1)
        Phi0_abs = sptensor(Phi0.subs, abs(Phi0.vals), size(Phi0));
        %lei's change: without summary of atoms. Instead, 
        % use the max of atoms
        f_max = [];
        for i = 1:size(Phi0_abs,3)
            f_max = [f_max;max(max(spmatrix(Phi0_abs(:,:,i))))];
        end
        %vsum = spmatrix(ttv(Phi0_abs, ones_a, 1));
%         v_max = max(vsum')'; %max among fascicles
        %f_max = max(vsum,[],1)';
        %v_max = max(vsum,[],2);
    else
        Phi0_abs = abs(Phi0);
        %vsum = ones_a' * spmatrix(Phi0);
        f_max = Phi0_abs';
        %v_max = vsum';
        v_max = max(spmatrix(Phi0_abs(:,:,1)),[],2);
    end
    
    [subs, vals] = find(Phi0);
%     vals = vals ./ v_max(subs(:,2));
    vals = vals ./ f_max(subs(:,3));
    Phi0 = sptensor(subs, vals, size(Phi0));
    
    %vcoords = roi.coords(subs(:,2), :);
    %angs = orient(subs(:,1), :);
    
    hs = {};
    
    vcenters = roi.coords(voxel_indices(1:size(Phi,2)), :);
    % find voxels that have non-zeros
    flag = ones(Na,1)' * spmatrix(ttv(Phi,ones(Nf,1),3));
    vcenters = vcenters(flag~=0,:);
%     hs.centers = scatter3(ax,vcenters(:,1),vcenters(:,2),vcenters(:,3), 1, '+');
    
    limits(2,:) = max(vcenters);
    limits(1,:) = min(vcenters);
    limits = limits';
    xlim(limits(1,:));
    ylim(limits(2,:));
    zlim(limits(3,:));
    
    for f = fidx
%         tsubs.subs = {':',':',f};
%         Phif = subsref(Phi, tsubs);
        
        Phif = Phi(:,:,f);
        Phi0f = Phi0(:,:,f);
        
        fprintf('.');
        if (mod(f, 100) == 0)
            fprintf('\n');
        end
        
        if (nnz(Phif) == 0)
            continue
        end

        [subs, vals] = find(Phif);
        
%         vcoords = roi.coords(subs(:,2), :);
        vcoords = roi.coords(voxel_indices(subs(:,2)), :);
        angs = (0.5*vals * ones(1,3)) .* orient(subs(:,1), :);
%         angs = 0.5*orient(subs(:,1), :);

%         quiver3(vcoords(:,1),vcoords(:,2),vcoords(:,3),angs(:,1),angs(:,2),angs(:,3),0,'b');
        h = quiver3(ax,vcoords(:,1),vcoords(:,2),vcoords(:,3),angs(:,1),angs(:,2),angs(:,3),0,'b','LineWidth',1','MaxHeadSize',0.8);
%         h = quiver3(ax,vcoords(:,1),vcoords(:,2),vcoords(:,3),angs(:,1),angs(:,2),angs(:,3),0.4,'b','LineWidth',1','MaxHeadSize',0.8);
        hs.fs(f,1) = h;
        
        [subs, vals] = find(Phi0f);
        
%         vcoords = roi.coords(subs(:,2), :);
        vcoords = roi.coords(voxel_indices(subs(:,2)), :);
        vcoords(:,2) = vcoords(:,2) + 0.1;
        angs = (0.5*vals * ones(1,3)) .* orient(subs(:,1), :);
%         angs = 0.5*orient(subs(:,1), :);

%         quiver3(vcoords(:,1),vcoords(:,2),vcoords(:,3),angs(:,1),angs(:,2),angs(:,3),0,'b');
        h = quiver3(ax,vcoords(:,1),vcoords(:,2),vcoords(:,3),angs(:,1),angs(:,2),angs(:,3),0,'r','LineWidth',1','MaxHeadSize',0.8);
%         h = quiver3(ax,vcoords(:,1),vcoords(:,2),vcoords(:,3),angs(:,1),angs(:,2),angs(:,3),0.4,'g','LineWidth',1','MaxHeadSize',0.8);
        
        hs.fs(f,1) = h;
        
        fgExpert = fe.fg.fibers{f}' + 1;
        b = plot3(fgExpert(:,1),fgExpert(:,2),fgExpert(:,3),'k','LineWidth',2);

        
%         endings = vcoords + angs;
%         scatter3(endings(:,1),endings(:,2),endings(:,3), '.');
    end
    fprintf('\n');
%     quiver3(vcoords(:,1),vcoords(:,2),vcoords(:,3),angs(:,1),angs(:,2),angs(:,3),0);
    %[subs, vals] = find(Phi);
    
end
