function [ind_X, ind_Y]=swap_label_plus(cluster_p,cluster_q,cluster_pt)
%normalize the cluster distribution
cluster_p0 = (cluster_p./sum(cluster_p,2))/sum(cluster_p,'all');
cluster_q0 = (cluster_q./sum(cluster_q,2))/sum(cluster_q,'all');
cluster_pt0 = (cluster_pt./sum(cluster_pt,2))/sum(cluster_pt,'all');



Dis = zeros(size(cluster_q0,1),size(cluster_p0,1));
for i = 1:size(cluster_q0,1)
    for j = 1:size(cluster_p0,1)
        Dis(i,j) = JSDiv_v3(cluster_q0(i,:),cluster_p0(j,:));
    end
end

[~,ind_Y] = min(Dis,[],2);

Dis = zeros(size(cluster_p0,1),size(cluster_pt0,1));
for i = 1:size(cluster_p0,1)
    for j = 1:size(cluster_pt0,1)
        Dis(i,j) = JSDiv_v3(cluster_p0(i,:),cluster_pt0(j,:));
    end
end

[~,ind_X] = min(Dis,[],2);

clearvars -except ind_X ind_Y dis