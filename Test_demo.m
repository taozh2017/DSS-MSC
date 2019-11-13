
% -------------------------------------------------------------------------
% Code for paper entitled "Dual Shared-Specific Multiview Subspace Clustering"
% $ History
% Created by T. Zhou, Dec., 2017.
% Released by T. Zhou, Nov., 2019.
% 
% NOTE: you need to tune parameters when use different datasets.
% -------------------------------------------------------------------------

clc;clear;close all;
addpath(genpath('.'));

% --------------------------------------- load data
data_id = 'bbcsport-mtv';
[X,gt]  = data_load(data_id);


% --------------------------------------- paras setting
paras.lambda = 0.1;
paras.beta   = 0.001;
paras.Ns     = 30;
paras.Nc     = 80;
        
% --------------------------------------- clustering 
[nmi,ACC,AR,f,p,r,RI] = DSS_MVC_Learning(X,gt,paras);
disp(['ACC=',num2str(ACC),'  NMI=',num2str(nmi)])

 

