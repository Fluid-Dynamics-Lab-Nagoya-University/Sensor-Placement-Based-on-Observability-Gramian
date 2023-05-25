 function [zhat, L, ztilde, Utilde,test] = F_sensor_gramianconvex_approxnt(C, k, maxiteration, A, z)
%{
<memo>
gradient
hessian


%}
    %% sens_sel_approxnt_vec
%
% see paper Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% Nov 2007 Siddharth Joshi & Stephen Boyd

% Newton's method parameters
% maxiteration = 200;
% n = 1000;
% r = 10;
% A = randn(n,r);
% k = 5;

MAXITER  = maxiteration;
NT_TOL = 1e-3;
GAP = 1.005;
% GAP = 1+10^-16;
% Backtracking line search parameters
alpha = 0.1;
beta = 0.5;
[m n l] = size(C);
%%  initialize
if isempty(z)
z = ones(m,1)*(k/m); 
end
%%
g = zeros(m,1);
ones_m = ones(m,1);
kappa = log(GAP)*n/m; 
% guarantees GM of lengths of semi-axes of ellipsoid corresponding to 
% ztilde <= 1.01 from optimal

%fprintf('\nIter.  Step_size  Newton_decr.  Objective  log_det\n');
CdC=zeros(n,n);
zmat=repmat(z,1,n);
for ll=1:l
%     UdU=UdU+((U(:,:,ltmp))'*diag(z)*U(:,:,ltmp));
%     dC=diag(z)*C(:,:,ll);
    dC=zmat.*C(:,:,ll);
%     dU=dU+diag(z)*U(:,:,Itemp);
    CdC=CdC+(C(:,:,ll))'*dC;
end
% fz = -log(det(AdA)) - kappa*sum(log(z) + log(1-z)); %<-------- determinant

%%   initialize - gram
[r,~]=size(A);
AAI1 = kron(A',A')-eye(r^2);  %<------- AIA1
iAAI1 = inv(AAI1);
vecCC1 = reshape(CdC,[],1);  %<------- CdC1
vecWo = - iAAI1 * vecCC1;
Wo = reshape(vecWo,[r,r]);
fz = -log(det(Wo)) - kappa*sum(log(z) + log(1-z));

%%   initializa - AccA
for i=1:m
    matCC(:,:,i) = C(i,:)'*C(i,:);
end
vecCC2 = reshape(matCC,r^2,[]);
vecAccA = - iAAI1 * vecCC2;
AccA = reshape(vecAccA,[r,r,m]);
%%
AAI2 = kron(A,A)-eye(r^2);
iAAI2 = inv(AAI2);
%%   computation
% %fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', -fz, log(det(AdA)));
% fprintf('\nIter.  Step_size  Newton_decr.  Objective  obs. gramian\n');
% fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', fz, det(Wo));

for iter=1:MAXITER
%     disp(['itr: ', num2str(iter)])
    CdC=zeros(n,n);
    for ll=1:l
       CdC=CdC+((C(:,:,ll))'*(zmat.*C(:,:,ll)));
    end
    %%
%     iCdC = inv(CdC);
%     V=zeros(m,m);
%     for ll=1:l
%        V = V + C(:,:,ll)*iCdC*(C(:,:,ll))';   
%     end
%     g = -diag(V)- kappa*(1./z - 1./(1-z));  %<-------- gradient
%     H = V.*V + kappa*diag(1./(z.^2) + 1./((1-z).^2));   %<-------- Hessian
    
    %% vector identity
    vecCC1 = reshape(CdC,[],1);
    vecWo = - iAAI1 * vecCC1;
    Wo = reshape(vecWo,[r,r]);
    iWo = inv(Wo);   %<------- inv of obs. gram.
    
    vecCC2 = reshape(iWo,[],1);
    vecAQA = - iAAI2 * vecCC2;
    AQA = reshape(vecAQA,[r,r]);
    
%     col_iWo=zeros(r,r,m);
%     for i=1:m
%         col_iWo(:,:,i)=iWo;
%     end    
%     col_iWo = repmat(iWo,[1,1,m]);
%     vecCC3_old = reshape(col_iWo .* AccA .* col_iWo,r^2,[]);
    vecCC3 = reshape(pagemtimes(iWo,pagemtimes(AccA,iWo)),r^2,[]);
    vecAQAccAQA = - iAAI2 *  vecCC3;
    AQAccAQA = reshape(vecAQAccAQA,[r,r,m]);
    H=zeros(m,m);
    for j=1:m
       Htemp = (C*AQAccAQA(:,:,j)).*C;
       H(:,j) = -sum(Htemp, 2);
%        H(:,i) = diag(C*AQAccAQA(:,:,i)*C');
    end
%     H = pagemtimes(C,AQAccAQA).*repmat(C,[1,1,m]);
       
    
    %% gradient & Hessian
%     g = -diag(C*AQA*C')- kappa*(1./z - 1./(1-z));  %<-------- gradient
    gtemp = (C*AQA).*C;
    g = sum(gtemp,2);
    
    g = - g - kappa*(1./z - 1./(1-z));  %<-------- gradient
%     H = H + kappa*diag(1./(z.^2) + 1./((1-z).^2));   %<-------- Hessian
    H = - H + kappa*sparse(1:m,1:m,(1./(z.^2) + 1./((1-z).^2)),m,m);   %<-------- Hessian

    %%
%     if det(H)>10^300
%         det(H);
%     end
    R = chol(H);    %<-------- R'*R=H inv(H)=inv(R)*inv(R)'
    Hinvg = (R\(R'\g)); %<-------- inv(Hessian)*gradient
    Hinv1 = (R\(R'\ones_m));    %<-------- inv(Hessian)*ones
    dz = -Hinvg + ((ones_m'*Hinvg) / (ones_m'*Hinv1))*Hinv1;    %<-------- search step

    deczi = find(dz < 0);
    inczi = find(dz > 0); 
    s = min([1; 0.99*[-z(deczi)./dz(deczi) ; (1-z(inczi))./dz(inczi)]]);    %<-------- detect overrun from defined range (0,1)
    
    while (1)
        %disp([zp, z, s*dz])
        zp = z + s*dz; % keep zi in defined range
        zpmat = repmat(zp,1,n);
        CdC=zeros(n,n);
        for ll=1:l
%            CdC=CdC+((C(:,:,ll))'*diag(zp)*C(:,:,ll));
           CdC=CdC+(C(:,:,ll))'*(zpmat.*C(:,:,ll));
        end
%%   gram 
        vecCC4 = reshape(CdC,[],1);
        vecWo = - iAAI1 * vecCC4;
        Wo = reshape(vecWo,[r,r]);

%%
%         fzp = -log(det(AdA)) - kappa*sum(log(zp) + log(1-zp));
%         %<-------- determinantCTC
        fzp = -log(det(Wo)) - kappa*sum(log(zp) + log(1-zp));  %<-------- determinant Gram
%         disp(['fzp:', num2str(fzp), ',fz:', num2str(fz), ',alpha*s*gpri*dz:', num2str(alpha*s*g'*dz)])
        if (fzp <= fz + alpha*s*g'*dz)  % break if update doesnot overrun???
            break;
        end
        s = beta*s; % s becomes smaller and try again
    end
    z = zp; fz = fzp;
    zmat = repmat(z,1,n); 
    CdC=zeros(n,n);
    for ll=1:l
       CdC=CdC+((C(:,:,ll))'*(zmat.*C(:,:,ll)));     
%        CdC=CdC+((C(:,:,ll))'*diag(z)*C(:,:,ll));        
    end
    
    
% %     fprintf('%4d %10.3f %10.3f %10.3f %10.3f\n', iter, s, -g'*dz/2, -fz, log(det(AdA)));
% %     fprintf('%4d %10.3f %10.3f %10.3f %10.3f\n', iter, s, -g'*dz/2, fz, det(Wo));
%     disp(['gdz: ', num2str(-g'*dz/2),'  maximum: ', num2str(max(z))])
    %%
    if(-g'*dz/2 <= 10*NT_TOL)   % break if the amount of change in obj. (given by gradient * step size) is smaller than threshold
        [mzp,mzpind]=max(zp);
%         disp([mzp,mzpind])
        [mzp,mzpind]=max(s*dz);
%         disp([mzp,mzpind])
        [mzp,mzpind]=max(g);
%         disp([mzp,mzpind])
%         disp(-g'*dz/2)
%         disp(max(zp),min(zp))
    end
    if(-g'*dz/2 <= NT_TOL)   % break if the amount of change in obj. (given by gradient * step size) is smaller than threshpld
        break;
    end
    test.obj(iter,1) = fzp;
    test.sens(iter,:) = maxk(z,k);
end
%     disp(['gdz: ', num2str(-g'*dz/2),'  maximum: ', num2str(max(z))])
if iter==MAXITER
%     disp(k)
%     disp('max iteration')
%     sprintf('n = %d: %s',size(C,1),'max iteration')
    disp(['LCR_n = ' num2str(size(C,1)) ' : max iteration'])
end
zsort=sort(z); 
thres=zsort(m-k); 
zhat=(z>thres);
CdC=zeros(n,n);
for ll=1:l
   CdC=CdC+((C(:,:,ll))'*(zmat.*C(:,:,ll)));      
%    CdC=CdC+((C(:,:,ll))'*diag(z)*C(:,:,ll));        
end
L = log(det(CdC));
ztilde = z; 
Utilde = log(det(CdC)) + 2*m*kappa;


