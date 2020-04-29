clear all 
close all

trans = [0.95,0.05; 0.10,0.90];
emis = [1/6 1/6 1/6 1/6 1/6 1/6;
   1/10 1/10 1/10 1/10 1/10 1/2];
[m,n]=size(emis);

% Condition 1:
if rank(trans)==length(trans) && rank(trans)==rank(emis)
    fprintf("all set\n")
else
    fprintf("Condition 1 not met\n")
end

obs=1000000;

[seq,states] = hmmgenerate(obs,trans,emis);
[estimateTR,estimateE] = hmmestimate(seq,states);

nobs=max(length(emis));
P1=zeros(nobs,1);
for i=1:obs
    P1(seq(i))=P1(seq(i))+1;
end
P1=P1/obs;

Pair21 = zeros(nobs,nobs);
for i=1:obs-1
    Pair21(seq(i+1),seq(i))=Pair21(seq(i+1),seq(i))+1;
end
Pair21(:,:)=Pair21(:,:)/(obs-1);
Pair3x1 = zeros(nobs,nobs,nobs);
for i=1:obs-2
    Pair3x1(seq(i+2),seq(i),seq(i+1))=Pair3x1(seq(i+2),seq(i),seq(i+1))+1;
end
Pair3x1(:,:,:)=Pair3x1(:,:,:)/(obs-2);

[U,S,V]=svd(Pair21);

figure()
semilogy(S, 'xb', 'LineWidth', 2, 'MarkerSize', 12);
grid on
title('Singular values model 00');

U1=U(1:n,1:m);
rank(U1'*emis');
B=zeros(m,m,nobs);

b1 = U1'*P1;
b_inf = pinv(Pair21'*U1)*P1;
for i=1:nobs
    B(:,:,i) = U1'*Pair3x1(:,:,i)*pinv(U1'*Pair21);
end

Pair31 = zeros(nobs,nobs);
for i=1:obs-2
    Pair31(seq(i+2),seq(i))=Pair31(seq(i+2),seq(i))+1;
end
Pair31(:,:)=Pair31(:,:)/(obs-2);

O=zeros(m,nobs);
for i=1:nobs
    O(:,i) = eig(U1'*Pair3x1(:,:,i)*pinv(U1'*Pair31));
end

pi = (pinv(O')*P1)';
T = (pinv(O')*Pair21*pinv(O')'/diag(pi))';



