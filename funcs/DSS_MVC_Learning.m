
% -------------------------------------------------------------------------
% Demo for paper entitled "Dual Shared-Specific Multiview Subspace Clustering";
% $ History
% Created by T. Zhou, Sep., 2017.
% Revised
% -------------------------------------------------------------------------

function [nmi,ACC,AR,f,p,r,RI] = DSS_MVC_Learning(X,gt,paras)

% ---------------------------------- data pre-processing
V = size(X,2);      % number of views
N = size(X{1},2);   % number of samples
cls_num  = size(unique(gt),1);

lambda = paras.lambda;
beta   = paras.beta;
Ns     = paras.Ns;
Nc     = paras.Nc;


for i = 1:V
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
end

% ---------------------------------------------- initialize variables

for i = 1:V
   
    Ep{i} = zeros(Ns+Nc,N);
    Eh{i} = zeros(Ns+Nc,N);
    Yp{i} = zeros(Ns+Nc,N);
    Yh{i} = zeros(Ns+Nc,N);
    Ys{i} = zeros(N,N);
    Hv{i} = rand(Nc,N);
    Zv{i} = zeros(N,N);  
end
Hs   = rand(Ns,N);
Z    = zeros(N,N);
Yz   = zeros(N,N);


% ----------------------------------------------
IsConverge = 0;
mu         = 1e-4;
pho        = 1.5;  % 1.5
max_mu     = 1e6;
max_iter   = 100;
iter       = 1;
thresh     = 1e-6;

% -------------------------------------------------------------------------
while (IsConverge == 0&&iter<max_iter+1)
    
    sum_DtD = zeros(N,N);
    sum_DtC = zeros(N,N);
    
    for i = 1:V
        
        % ---------------------------------- update Pv
        P{i} = updatePP(Yp{i},mu,[Hs;Hv{i}] + Ep{i}, X{i});
        
        % ---------------------------------- update Hv
        ZZ    = (eye(N)-Z-Zv{i});
        Hv{i} = (P{i}(Ns+1:end,:)*X{i} - Ep{i}(Ns+1:end,:) + Yp{i}(Ns+1:end,:)/mu + Eh{i}(Ns+1:end,:)*ZZ' - Yh{i}(Ns+1:end,:)*ZZ'/mu) / (eye(N)+ZZ*ZZ'+eye(N)*1e-8 ); % revised 09/19
        
        % ---------------------------------- update Zv
        tepB  = [Hs;Hv{i}];
        tepA  = [Hs;Hv{i}] - [Hs;Hv{i}]*Z - Eh{i};
        Qv{i} = ComputeSaprseValues(Zv{i} + Ys{i}/mu, beta/mu);
        Zv{i} = (tepB'*tepB + eye(N) ) \ (tepB'*tepA + Qv{i} - diag(diag(Qv{i})) + (tepB'*Yh{i} - Ys{i})/mu);
        clear tepB tepA;
        
        % ---------------------------------- update Ev
        G = [P{i}*X{i} - [Hs;Hv{i}] + Yp{i}/mu; [Hs;Hv{i}] - [Hs;Hv{i}]*Z - [Hs;Hv{i}]*Zv{i} + Yh{i}/mu];
        E = solve_l1l2(G,lambda/mu);
        
        Ep{i} = E(1:Ns+Nc,:);
        Eh{i} = E(Ns+Nc+1:end,:);
        
        % ----------------------------------
        tepC    = [Hs;Hv{i}] - [Hs;Hv{i}]*Zv{i} - Eh{i};
        tepD    = [Hs;Hv{i}];
        sum_DtD = sum_DtD + tepD'*tepD;
        sum_DtC = sum_DtC + (tepD'*(tepC + Yh{i}/mu));
        
    end
    
    % ---------------------------------- update J
    J = softth(Z+Yz/mu,1/mu);  %  + eye(N)*1e-8
    
    % ---------------------------------- update Z
    Z  =  inv(eye(size(sum_DtD,2)) + sum_DtD + eye(N)*1e-8)* (J - Yz/mu + sum_DtC);
    
    % ---------------------------------- update H
    tepHs = zeros(Ns,N);
    tepZ  = zeros(N,N);
    for i = 1:V
        Zvv   = eye(N) - Z - Zv{i};
        tepHs = tepHs + (P{i}(1:Ns,:)*X{i} - Ep{i}(1:Ns,:) + Yp{i}(1:Ns,:)/mu + Eh{i}(1:Ns,:)*Zvv' - Yh{i}(1:Ns,:)*Zvv'/mu);
        tepZ  = tepZ + eye(N) + Zvv*Zvv';
    end
    Hs  = tepHs * inv(tepZ + eye(N)*1e-8);   % revise 0926 tepHs / (tepZ + eye(N)*1e-8);
    
    % ---------------------------------- updata multipliers
    for i = 1:V
        Yp{i} = Yp{i} + mu*(P{i}*X{i}  - [Hs;Hv{i}]   - Ep{i});
        Yh{i} = Yh{i} + mu*([Hs;Hv{i}] - [Hs;Hv{i}]*Z - [Hs;Hv{i}]*Zv{i} - Eh{i});
        Ys{i} = Ys{i} + mu*(Zv{i} - Qv{i} + diag(diag(Qv{i})));
    end
    Yz = Yz + mu*(Z-J);
    mu = min(pho*mu, max_mu);
    
    % ----------------------------------- convergence conditions
    min_err = 0;
    for i = 1:V
        errp(i) = norm(P{i}*X{i}  - [Hs;Hv{i}]   - Ep{i},inf);
        errh(i) = norm([Hs;Hv{i}] - [Hs;Hv{i}]*Z - [Hs;Hv{i}]*Zv{i} - Eh{i},inf);
        errs(i) = norm(Zv{i} - Qv{i} + diag(diag(Qv{i})) , inf);
        
    end
    max_err0 = max([errp(:);errh(:);errs(:)]);
    max_err  = max([max_err0,norm(Z-J,inf)]);
    
    % -----------------------------------
    if max_err < thresh
        IsConverge = 1;
    end
    cov_val(iter) = max_err; % norm(Z-J,inf);
    iter          = iter + 1;
    
end

% -------------------------------------------------------------------------
Z_all = zeros(size(Z));

for j = 1:V
    Zi = Zv{j};
   
    Zi    = Z + Zi;
    Z_all = Z_all + (abs(Zi) + abs(Zi'));
    
end

Z_all = Z_all / V;
[nmi,ACC,AR,f,p,r,RI] = clustering(Z_all, cls_num, gt);
