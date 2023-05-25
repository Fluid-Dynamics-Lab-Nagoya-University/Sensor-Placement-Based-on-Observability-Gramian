function f = F_calc_det(W,maxrank)

%{
    *** NOTE ***
    W must be a symmetric matrix, such like FIM and Gramians.
    if maxrank is used, calculation considers the 
%}

% if ~issymmetric(W(:,:,1))
%     error('W is not symmetric matrix.... aborted')
% end

if isempty(maxrank)
    maxrank = min(size(W,1),size(W,2));
end

if maxrank < size(W,1)
    for j = 1:size(W,3)
        if ~isnan(W(:,:,j)) 
            [~,S] = eigs(W(:,:,j),maxrank);
            f(j) = det(S);
        else
            f(j) = NaN;
        end
    end
else
    for j = 1:size(W,3)
        if ~isnan(W(:,:,j))
            %         [~,S] = eigs(W(:,:,j),maxrank);
            f(j) = det(W(:,:,j));
        else
            f(j) = NaN;
        end
    end
end


end