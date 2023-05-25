function [Wo]=F_calc_Gram_dlyap(A,U,sensors,iR)
[ps,~]  = size(sensors);
[n,r]   =size(U);
C       = U(sensors,:);
if isempty(iR)
    iR = eye(ps);
end
Wo      = dlyap(A',C'*iR*C);
end