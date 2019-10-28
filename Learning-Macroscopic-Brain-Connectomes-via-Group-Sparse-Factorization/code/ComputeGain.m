% called by group_reg_greedy_atom_selection.m
% calculate gained information by including new atoms
function [gS,V] = ComputeGain(BTB,C,S,cInv,b)
gS = [];
V = [];
for a = 1:size(BTB,1)
    cSa = BTB(S,a);
    c = cInv*cSa;
    v = 1/(C(a)-cSa'*c);
    V = [V,v];
    gS =[gS;v*(b(S)'*c-b(a))^2];
end