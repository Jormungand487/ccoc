import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score as Cal_ARI

def Cal_NMI(clusterlabels, truelabels):
    LC = set(clusterlabels)
    LT = set(truelabels)
    l = len(clusterlabels)
    I = 0
    Hc = 0
    Ht = 0
    for c in LC:
        Pc = clusterlabels.count(c)/l
        Hc -= Pc * np.log(Pc)
        for t in LT:
            Pct = sum([(a == c) and (b == t) for a, b in zip(clusterlabels, truelabels)])/l
            if Pct != 0:
                Pt = truelabels.count(t)/l
                I += Pct * np.log(Pct / (Pc * Pt))
    for t in LT:
        Pt = truelabels.count(t)/l
        Ht -= Pt * np.log(Pt)
    return I / np.sqrt(Hc * Ht)

def Cal_purity(cluster_labels, true_labels):
    n = len(true_labels)
    
    # 统计每个簇中各类别样本的数量
    cluster_count = {}
    for cluster_id, true_label in zip(cluster_labels, true_labels):
        if cluster_id not in cluster_count:
            cluster_count[cluster_id] = {}
        if true_label not in cluster_count[cluster_id]:
            cluster_count[cluster_id][true_label] = 0
        cluster_count[cluster_id][true_label] += 1
    
    # 计算每个簇中最多的类别样本数量，然后求和
    total = 0
    for cluster_id in set(cluster_labels):
        total += max(cluster_count[cluster_id].values())
    
    # 纯度计算
    purity = total / n
    return purity

def rand_index(cluster_labels, true_labels):
    n = len(true_labels)
    a, b = 0, 0
    
    for i in range(n):
        for j in range(i + 1, n):
            same_true = (true_labels[i] == true_labels[j])
            same_cluster = (cluster_labels[i] == cluster_labels[j])
            if same_true and same_cluster:
                a += 1
            elif not same_true and not same_cluster:
                b += 1
    
    # 计算RI
    ri = (a + b) / (n * (n - 1) / 2)
    return ri

a = pd.read_csv("label.csv",header = None)
a = a[0].tolist()
# a = a.to_numpy()
b = pd.read_csv("pred.csv",header = None)
# b = b.to_numpy()
b = b[0].tolist()
print(f"purity is {Cal_purity(a,b)}.")
# print(f"RI is {rand_index(a,b)}.")
print(f"ARI is {Cal_ARI(a,b)}.")
print(f"NMI is {Cal_NMI(a,b)}.")