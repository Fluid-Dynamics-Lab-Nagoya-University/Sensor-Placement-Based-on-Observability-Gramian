function [z, cvx_cputime, cvx_optbnd, cvx_optval,...
    cvx_slvitr, cvx_slvtol, cvx_status] = F_cvx_gramian_sdp_lmi(A, C, p)
solver_now = cvx_solver;
[r,~] = size(A);
[n,~] = size(C);
cvx_solver('Mosek')
% cvx_solver_settings('MSK_IPAR_NUM_THREADS',1)
cvx_begin sdp
    variable z(n) nonnegative
    variable Z(n,n) symmetric
    variable X(r,r) symmetric semidefinite
    maximize( det_rootn(X) )
    subject to 
        A'*X*A-X >= -C'*(repmat(z,1,r).*C);
        trace(Z) <= p;
        diag(Z) == z;
        [Z,z;z',1] >= 0;
cvx_end
cvx_solver(solver_now);
end