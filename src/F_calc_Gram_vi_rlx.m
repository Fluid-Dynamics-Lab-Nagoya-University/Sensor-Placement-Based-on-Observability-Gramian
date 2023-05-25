function [Wo] = F_calc_Gram_vi_rlx(A,U,z)
% [ps,~]=size(sensors);
[n,r]=size(U);
% [H]=F_calc_sensormatrix(1, ps, n, sensors);
AAI = kron(A',A')-eye(r^2);
iAAI = inv(AAI);
CC = U'*(repmat(z,1,r).*U);
vecCC = reshape(CC,[],1);  %<------- CdC1
vecWo = - iAAI * vecCC;
Wo = reshape(vecWo,[r,r]);
end