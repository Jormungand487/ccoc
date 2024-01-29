clear
clc
close all

%load the processed data in example 3 in the paper
source = 'ex14';
p = readmatrix(['../../data/',source,'/mS.csv']);
p = p(2:end, 2:end);
p = p';
q = readmatrix(['../../data/',source,'/mT.csv']);
q = q(:, 2:end);
q = q';
q0 = readmatrix(['../../data/',source,'/mU.csv']);
q0(1,:) = [];
q0(:,3) = 1;
q0 = spconvert(q0);
q0 = q0';

Cx_truth = readtable(['../../data/',source,'/lS.csv']);
Cx_truth = Cx_truth(:,2);
[~,~,Cx_truth] = unique(Cx_truth, 'stable');

Cy_truth = readtable(['../../data/',source,'/lT.csv']);
Cy_truth = Cy_truth(:,2);
[~,~,Cy_truth] = unique(Cy_truth, 'stable');



pp = parpool(64);
%%coupleCoC+
%setting the values of hyperparameters
nrowcluster1=6;nrowcluster2=6;ncolcluster=15;ncolcluster0=6;iter=20;
lambda=1.5;beta=1.0;gamma=1;nsub=6;
[Cy, Cz, Cz0, cluster_p, cluster_q, cluster_q0] = coupleCoC_plus(p,q,q0,Cx_truth,nrowcluster2,ncolcluster,ncolcluster0,iter,lambda,beta);

%%results
[TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx_truth, Cy);


[ind_X, ind_Y]=swap_label_plus(cluster_p,cluster_q,nsub);
delete(pp);

save(['../../data/',source,'Cy.mat'],'Cy');
save(['../../data/',source,'Cz.mat'],'Cz');
save(['../../data/',source,'Cz0.mat'],'Cz0');
save(['../../data/',source,'TAB_Y.mat'],"TAB_Y");
save(['../../data/',source,'Eval_tab.mat'],"Eval_tab");
save(['../../data/',source,'ind_X.mat'],"ind_X");
save(['../../data/',source,'ind_Y.mat'],"ind_Y");




