%load('../data/fe_struct_with_predicted_signal_from_Arcuate_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33.mat')
load '../data/fe_struct_with_predicted_signal_from_ARC_SLF_normFP_96dirs_b2000_1p5iso_PROB_lmax10_NUM01_L33'

load('../data/stage1/Phi_OMP_Large.mat');
Phi_sp1 = Phi;

load('../data/stage1/Phi_GD_Large.mat');
Phi_sp2 = Phi;

load('../data/newsubsets/compact_Phi.mat');

orient = fe.life.M.Atoms.orient;

orient_norm = orient.* orient;

orient_norm = sum(orient_norm);

orient_norm = sqrt(orient_norm);

norm_mat = orient_norm' * orient_norm;

inner_product = abs(orient' * orient);

inner_product = inner_product ./ norm_mat;

Na = size(Phi_sp1,1);
Nv = size(Phi_sp1,2);
Nf = size(Phi_sp1,3);

angs_diff = [];
for v = 1:Nv
    myV = full(spmatrix(Phi_sp1(:,v,:)));
    myVExp = full(spmatrix(Phi(:,v,:)));
    myVF = find(sum(myV)~=0);
    myVExpF = find(sum(myVExp)~=0);
    % the fascicles active for both expert and our phi
    vIS = intersect(myVF, myVExpF);
    for i = 1:length(vIS)
        myF = myV(:,vIS(i));
        myFExp = myVExp(:,vIS(i));
        myFindex = find(myF);
        myFExpIndex = find(myFExp);
        temp = inner_product(myFindex,myFExpIndex);
        res = max(temp,[],1);
        angs_diff = [angs_diff,res];
    end    
end
sum(real(acos(angs_diff))*180/pi)
histogram(real(acos(angs_diff))*180/pi,'Normalization','probability','DisplayStyle','stairs');
hold on

Na = size(Phi_sp2,1);
Nv = size(Phi_sp2,2);
Nf = size(Phi_sp2,3);

angs_diff = [];
for v = 1:Nv
    myV = full(spmatrix(Phi_sp2(:,v,:)));
    myVExp = full(spmatrix(Phi(:,v,:)));
    myVF = find(sum(myV)~=0);
    myVExpF = find(sum(myVExp)~=0);
    % the fascicles active for both expert and our phi
    vIS = intersect(myVF, myVExpF);
    for i = 1:length(vIS)
        myF = myV(:,vIS(i));
        myFExp = myVExp(:,vIS(i));
        myFindex = find(myF);
        myFExpIndex = find(myFExp);
        temp = inner_product(myFindex,myFExpIndex);
        res = max(temp,[],1);
        angs_diff = [angs_diff,res];
    end
end
sum(real(acos(angs_diff))*180/pi)
histogram(real(acos(angs_diff))*180/pi,'Normalization','probability','DisplayStyle','stairs');

xlabel('Angle Differences of Fascicle Directions');
ylabel('P(Angle Differences of Fascicle Directions)');
legend('OMP', 'Greedy');