function [sensors,pfull] = F_sensor_GG(A,Corg,sensors,p,F_obj)
%{
    Greedy selection with any onjective function provided. Note that F_obj
    is to be maximized, so "-" to minimizing optimization s.t. -tr(inv(W)).
%}

[ps,~]      = size(sensors);
[n,r]       = size(Corg);

Cprev       = Corg(sensors,:);
rankflag    = false;
pfull       = p;

for k = (ps+1):p
    %%   
    if ~rankflag %! for rank deficient case
        maxrank = 0;
        for i = 1:n
            if ismember(i,sensors)
                G(:,:,i) = NaN(r,r);
            else
                C = [Cprev;Corg(i,:)];
                G(:,:,i) = dlyap(A',C'*C);
                if isnan(G(:,:,i))
                    disp([i k])
                end
                maxrank=max(rank(G(:,:,i)),maxrank);
            end
        end
        if maxrank == r
            rankflag    = 1;
            pfull       = k;
        end
    else %! for rank sufficient case
        for i = 1:n
            if ismember(i,sensors)
                G(:,:,i) = NaN(r,r);
            else
                C = [Cprev;Corg(i,:)];
                G(:,:,i) = dlyap(A',C'*C);
                if isnan(G(:,:,i))
                    disp([i k])
                end
            end
        end
        maxrank = r;
    end
    %%   
    obj  = F_obj(G,maxrank);
    [M(k),Im]   = max(obj);
    sensors(k,1) = Im;    
    Cprev       = [Cprev;Corg(Im,:)]; 
end

end