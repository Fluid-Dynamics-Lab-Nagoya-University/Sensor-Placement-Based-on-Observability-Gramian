function [zhat, L, ztilde, Utilde, data, iter] = F_sensor_gramianCRSNC_dlyap(C, k, MAXITER, A, dnm)
%{
Solves the problem
	maximize log det (sum_{i=1}^m z_i Wo_i + kappa sum_{i=1}^m(log(z_i)+ log(1-z_i))
	subject to sum(z) = k
			   0 <= z_i <= 1, i=1,..., m
variable z in R^m
problem parameters kappa (>0), a_1, ..., a_m in R^n

for original work for static model, see paper Sensor Selection via Convex
Optimization, Nov 2007 Siddharth Joshi & Stephen Boyd
www.stanford.edu/~boyd/papers/sensor_selection.html

%}
%% Newton's method parameters
%MAXITER  = 1000;
NT_TOL = 1e-3;
GAP = 1.005;
%% Backtracking line search parameters
alpha = 0.1;
beta = 0.5;

[m, r, dimC] = size(C);
z = ones(m,1)*(k/m); % initialize
g = zeros(m,1);
ones_m = ones(m,1);
kappa = log(GAP)*r/m; 
% guarantees GM of lengths of semi-axes of ellipsoid corresponding to 
% ztilde <= 1.01 from optimal

% fprintf('\nIter.  Step_size  Newton_decr.  Objective  log_det log_det_dis time \n');
CdC=zeros(r,r);
zmat=repmat(z,1,r);
for ll=1:dimC
    CdC=CdC+(C(:,:,ll))'*(zmat.*C(:,:,ll));
end
%%   initialize - gram
Wo = dlyap(A',CdC);
fz = -log(det(Wo)) - kappa*sum(log(z) + log(1-z));

%% ITERATION
% fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', -fz, log(det(C'*diag(z)*C)));
iflg=0;
l=round(m/dnm);
lh=round(l/2);
idata=0;
iflgconverge=max(floor(m/l),1);
tic
for iter=1:MAXITER
    %% sketch
    [~,indz]=sort(z);
    S = sparse(1:l,[indz(end-lh+1:end); indz(randperm(m-lh,l-lh))],true,l,m);

    ones_l=ones(l,1);
    Cs = S*C;
    zz = S*z;
    
    %%  Sketched FIM and Obs. Gram
    iWo = inv(Wo);   %<------- inv of obs. gram.
    %% Sketched gradient
    AQA = dlyap(A,iWo);
    
    gtemp = (Cs*AQA).*Cs;
    g = sum(gtemp,2);
    g = - g - kappa*(1./zz - 1./(1-zz));  %<-------- gradient
    %% Sketched Hessian
    H = zeros(l,l);
    for j = 1:l
        AQAccAQA = dlyap(A,iWo*dlyap(A',Cs(j,:)'*Cs(j,:))*iWo); %! O(r^3)
        Htemp = (Cs*AQAccAQA).*Cs;                              %! O(nr^2)
        H(:,j) = -sum(Htemp, 2);        
    end

    H = - H + kappa*sparse(1:l,1:l,(1./(zz.^2) + 1./((1-zz).^2)),l,l);   %<-------- Hessian
    
    %%
    R = chol(H);                                                %! O(n^3) 
    Hinvg = (R\(R'\g));                                         %! O(n^3)
    Hinv1 = (R\(R'\ones_l));                                    %! O(n^3)
    dS = -Hinvg + ((ones_l'*Hinvg) / (ones_l'*Hinv1))*Hinv1;
    %% Stretched Newton step
    dz= S'* dS;
    deczi = find(dz < 0);
    inczi = find(dz > 0);
    s = min([1; 0.99*[-z(deczi)./dz(deczi) ; (1-z(inczi))./dz(inczi)]]);

    while (1)
        zp = z + s*dz;
        zpmat = repmat(zp,1,r);
        CdC=zeros(r,r);
        for ll=1:dimC
           CdC=CdC+(C(:,:,ll))'*(zpmat.*C(:,:,ll));
        end
        Wo = dlyap(A',CdC);                                     %! O(r^3)
        fzp = -log(det(Wo)) - kappa*sum(log(zp) + log(1-zp));  %<-------- determinant Gram
        
        if (fzp <= fz + alpha*s*g'*dS)
            break;
        end
        s = beta*s;
    end
    z = zp; fz = fzp;% CdC=CdCp;
    zsort=sort(z); thres=zsort(m-k); zhat=(z>thres);
    fzhat= -log(det(C(zhat,:)'*C(zhat,:)));
    if mod(iter,10)==0 
    idata=idata+1;
    time10=toc;
    data(idata,:)=[iter -fz log(det(CdC)) -fzhat time10];
    end 
    if(-g'*S*dz/2 <= NT_TOL)
        iflg=iflg+1;
        %break;
    else
        iflg=0;
    end
    if(iflg == iflgconverge)
        disp([mfilename ': CONVERGE @ L=' num2str(l) 'in iter=' num2str(iter)])
        break
    end
end

zsort=sort(z); thres=zsort(m-k); zhat=(z>thres);
dz = sparse(find(zhat),find(zhat),1,m,m);
L = log(det(C'*dz*C));
ztilde = z; 
Utilde = log(det(C'*dz*C)) + 2*m*kappa;
if iter == MAXITER
    disp([mfilename ': ENDED    @ L=' num2str(l) ' in iter=' num2str(iter)])
end

time10=toc;
% fprintf('%4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n', iter, s, -g'*dS/2, -fz, log(det(CdC)),-fzhat,time10);
idata=idata+1;
data(idata,:)=[iter -fz log(det(CdC)) -fzhat time10];