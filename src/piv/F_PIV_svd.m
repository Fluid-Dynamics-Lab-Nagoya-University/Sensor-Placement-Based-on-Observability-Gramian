function [u, loc] = F_PIV_svd(snap_ind, PIVpath, case_num)
PIVdata = matfile(PIVpath); size(PIVdata,'ufield_allt');
fname = fieldnames(PIVdata);

u.case      = case_num;
% u.field     = PIVdata.(char(fname(end)))(:,:,snap_ind);
u.field     = PIVdata.ufield_allt(:,:,snap_ind);
[sy,sx,sz]  = size(u.field);
u.vector    = reshape(u.field, sy*sx,[]);
u.mean      = mean(u.vector,2);

loc.NaN     = isnan(u.mean);
loc.Inf     = isinf(u.mean);
loc.normal  = ~(loc.NaN|loc.Inf);

u.fluc      = u.vector-repmat(u.mean, [1 sz]);
[u.Uorg, u.Sorg, u.Vorg] = svd(u.fluc(loc.normal,:), 'econ');

end