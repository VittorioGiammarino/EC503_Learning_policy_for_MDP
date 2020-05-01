function y = DataEncoding(Xtr,P)
% DataEncoding encodes the training data before computing the features in
% the Apprenticeship learning

[n,m]=size(Xtr);
for i=1:n
        [~,I]=max(P(:,Xtr(i,1),Xtr(i,2)));
        Y(i) = I;
end

y = Y(find(Y~=0))';

end