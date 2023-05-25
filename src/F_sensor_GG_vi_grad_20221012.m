function [sensors]=F_sensor_GG_vi_grad_20221012(Aorg,Corg,sensors,p,amp_reg)
%% derivative of gram matrix
[~,r] = size(Corg);
% [r,~] = size(Aorg);
[ps,~]= size(sensors);
% H=zeros(p,n);
if ps
    C = Corg(sensors,:);
    Q = dlyap(Aorg',C'*C);
    Qprior = zeros(size(C'*C));
else
    Qprior = amp_reg*eye(r,r);
    Q = Qprior;
end
%%
% AAI = kron(Aorg,Aorg)-eye(r^2);  %<------- AIA1
% iAAI = inv(AAI); 
% for nn =1: n
%     vecCC = reshape(Corg(nn,:)'*Corg(nn,:),[],1);  %<------- CdC1
%     vecWo = - iAAI1 * vecCC;
%     Gram = reshape(vecWo,[r,r]);
% %     sys = ss(Aorg,Borg,Corg(nn,:),Dorg(nn,:),1);
% %     test.gram = gram(sys,'o');
%     test.rank(nn) = rank(Gram);
% end
%%
for k=ps+1:p
%     Qinv=inv(Q+eps*eye(r,r));
    Qinv=inv(Q);% sometimes singular thus "inv" causes error
%     vecCC = reshape(Qinv,[],1);  %<------- CdC1
%     vecWo = - iAAI * vecCC;
%     Gram = reshape(vecWo,[r,r]);
    
    Gram = dlyap(Aorg,Qinv);
    
    %%
%     obj=zeros(n,1);
    obj = sum((Corg*Gram).*Corg,2);
%     for i=1:n
%         obj(i)=Corg(i,:)*Gram*(Corg(i,:))';
%     end
    for j=1:k-1
        obj(sensors(j,1))=0;
    end
    obj;
    [M(k),Im]=max(obj);
    H(k,Im)=1;
    sensors(k,1)=Im;
%     vecCC = reshape(Corg(sensors,:)'...
%         *Corg(sensors,:),[],1);  %<------- CdC1
%     vecWo = - iAAI * vecCC;
%     Q = reshape(vecWo,[r,r]) + Qprior;

    C = Corg(sensors,:);
    Q = dlyap(Aorg',C'*C) + Qprior;

%     if rank(Q) < r
%         disp([mfilename ': deficient'])
%     end
%     C=H*Corg;
%     D=H*Dorg;
%     sys=ss(Aorg,Borg,C,D,1);
%     Q=(gram(sys,'o'));
%     Det2(k,itr_data)=det(Q);
    rank2(k)=rank(Q);
end
end