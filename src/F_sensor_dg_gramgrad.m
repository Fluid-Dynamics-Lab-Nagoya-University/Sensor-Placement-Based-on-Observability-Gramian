function [sensors] = F_sensor_dg_gramgrad(Aorg,Corg,sensors,p,amp_reg)
%% derivative of gram matrix
[~,r] = size(Corg);
[ps,~]= size(sensors);
if ps
    C = Corg(sensors,:);
    Q = dlyap(Aorg',C'*C);
    Qprior = zeros(size(C'*C));
else
    Qprior = amp_reg*eye(r,r);
    Q = Qprior;
end
%%
for k=ps+1:p
    Qinv=inv(Q); %! sometimes singular thus "inv" causes error    
    Gram = dlyap(Aorg,Qinv);
    
    %%
    obj = sum((Corg*Gram).*Corg,2);
    for j=1:k-1
        obj(sensors(j,1))=0;
    end
    [~,Im]=max(obj);
    sensors(k,1)=Im;

    C = Corg(sensors,:);
    Q = dlyap(Aorg',C'*C) + Qprior;
end
end