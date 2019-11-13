

function [nmi,ACC,AR,f,p,r,RI] = clustering(S, cls_num, gt)

C              = SpectralClustering(S,cls_num);
[A,nmi,avgent] = compute_nmi(gt,C);
ACC            = Accuracy(C,double(gt));
[f,p,r]        = compute_f(gt,C);
[AR,RI,MI,HI]  = RandIndex(gt,C);



