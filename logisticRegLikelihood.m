function [fval, gradval] = logisticRegLikelihood(weight, X, C, m, lambda)

% Compute cost function and gradient

[number_samples, number_features] = size(X);
number_features = number_features/m;

gradval = zeros(number_features,1);

W = kron(eye(m),weight);
XWc = X*W; %inner product Training set weights

fpt2 = sum(XWc(sub2ind([number_samples,m],[1:number_samples]',C)));
XWc = exp(XWc);
sXWc = sum(XWc,2);
fpt1 = sum(log(sXWc));
fval = (fpt1 - fpt2)/number_samples + lambda * sum(abs(weight)); % Cost function

XWc = XWc./repmat(sXWc,1,m);

for i = 1:m
  gradval = gradval + ...
    ((XWc(:,i) - (C==i))'*X(:,(i-1)*number_features+1:i*number_features))';    
end

gradval = gradval/number_samples + lambda * sign(weight); %gradient
end
