function [sensors,H,Det,pfull]=F_sensor_GG_vi_20221002(A,Corg,sensors,p)
[ps,~]=size(sensors);
[n,r]=size(Corg);
[H]=F_calc_sensormatrix(1, ps, n, sensors);
Cprev = H*Corg;
AAI1 = kron(A',A')-eye(r^2);
iAAI1 = inv(AAI1);
rankflag = 0;
pfull = p;
for k=(ps+1):p
    %%   
    if ~rankflag
        maxrank=0;
        for i=1:n
            if ismember(i,sensors)
                G(:,:,i)=zeros(r,r);
            else
                C = [Cprev;Corg(i,:)];   % Greedy
                %             D = [Dprev;Dorg(i)];    % Greedy
                %         [ iR_temp ]= F_calc_Rinv_scalar(Un,Sn,[sensors; i],iR);

%                 vecCC = reshape(C'*C,[],1);  %<------- CdC1
%                 vecWo = - iAAI1 * vecCC;
%                 G(:,:,i) = reshape(vecWo,[r,r]);
                G(:,:,i) = dlyap(A',C'*C);
                if isnan(G(:,:,i))
                    disp([i k])
                end
                maxrank=max(rank(G(:,:,i)),maxrank);
            end
        end
        if maxrank == r
            rankflag = 1;
            pfull = k;
        end
    else
        for i=1:n
            if ismember(i,sensors)
                G(:,:,i)=zeros(r,r);
            else
                C = [Cprev;Corg(i,:)];   % Greedy
                %             D = [Dprev;Dorg(i)];    % Greedy
                %         [ iR_temp ]= F_calc_Rinv_scalar(Un,Sn,[sensors; i],iR);

                vecCC = reshape(C'*C,[],1);  %<------- CdC1
                vecWo = - iAAI1 * vecCC;
                G(:,:,i) = reshape(vecWo,[r,r]);
                if isnan(G(:,:,i))
                    disp([i k])
                end
%                 maxrank=max(rank(G(:,:,i)),maxrank);
            end
        end
        maxrank = r;
    end
    %%   
    obj=zeros(n,1);
    for i=1:n
        Gram=G(:,:,i);
        
        [Ug,Sg,Vg]=svd(Gram);
        Ugr=Ug(:        ,1:maxrank);
        Sgr=Sg(1:maxrank,1:maxrank);
        Vgr=Vg(:        ,1:maxrank);
        obj(i)=(det(Vgr'*(Ugr*Sgr)));
    end
    obj;
    [M(k),Im]=max(obj);
   sensors(k,1)=Im;
%     iR = F_calc_Rinv_vec(Unoi,Snoi,sensors,1,iR);
%     [ iR ]= F_calc_Rinv_scalar(Un,Sn,sensors,[]);
    Gram=G(:,:,Im);
    M(k);
    Det(k)=det(Gram);
    rank_gram(k)=rank(Gram);
    H(k,Im)=1;
    
    Cprev = [Cprev;Corg(Im,:)]; 
end

end