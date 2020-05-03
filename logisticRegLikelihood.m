function [fval, gradval] = logisticRegLikelihood(weight, X, C, m, lambda)

% Compute cost function and subgradient

[number_samples, number_features] = size(X);
number_features = number_features/m;

gradval = zeros(number_features,1);

%NLL(theta)=sum_samples(log(sum_actions(exp(theta'phi(x_i,u_j))))) ...
% - sum_sample(theta'phi(x_i,u_i))

W = kron(eye(m),weight); %theta
XWc = X*W; %theta'phi(x_i,u_i)

fpt2 = sum(XWc(sub2ind([number_samples,m],[1:number_samples]',C))); %sum_sample(theta'phi(x_i,u_i))
XWc = exp(XWc); %exp(theta'phi(x_i,u_i)) --> policy
sXWc = sum(XWc,2); %sum_actions(exp(theta'phi(x_i,u_j)))
fpt1 = sum(log(sXWc)); %sum_samples(log(sum_actions(exp(theta'phi(x_i,u_j)))))
fval = (fpt1 - fpt2)/number_samples + lambda * sum(abs(weight)); % Cost function = NLL(theta)

XWc = XWc./repmat(sXWc,1,m); %normalize policy

%Grad NLL = -sum_samples(phi - sum_actions(XWc*phi))
for i = 1:m
  gradval = gradval + ...
    ((XWc(:,i) - (C==i))'*X(:,(i-1)*number_features+1:i*number_features))'; % sum_samples[sum_actions(XWc(action)*phi(x_i,u_j)] - sum_samples(phi(x_i,u_i))
end

gradval = gradval/number_samples + lambda * sign(weight); %subgradient required since L1 is not differentiable
end
