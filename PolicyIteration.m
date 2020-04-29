function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Policy iteration
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (k x k x L)-matrix containing the transition probabilities
%           between all states in the state space for all control inputs.
%           The entry P(i, j, l) represents the transition probability
%           from state i to state j if control input l is applied.
%
%       G:
%           A (k x L)-matrix containing the stage costs of all states in
%           the state space for all control inputs. The entry G(i, l)
%           represents the cost if we are in state i and apply control
%           input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (k x 1)-matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (k x 1)-matrix containing the index of the optimal control
%       	input for each element of the state space. Mapping of the
%       	terminal state is arbitrary (for example: HOVER).

global K HOVER

%% Handle terminal state

global TERMINAL_STATE_INDEX

L=5;

P(TERMINAL_STATE_INDEX,:,:)=[];
P(:,TERMINAL_STATE_INDEX,:)=[];
G(TERMINAL_STATE_INDEX,:)=[];

%%% INITIALIZATION POLICY
mu=zeros(K-1,2);
mu(:,1)=HOVER;

%%% COMPUTATION OF G FOR THE INITIAL POLICY
GG=zeros(K-1,1);
for k=1:K-1
    GG(k)=G(k,mu(k,1));
end

%%% COMPUTATION OF P FOR INITIAL POLICY
PP=zeros(K-1);
for k=1:K-1
    for j=1:K-1
        PP(k,j)=P(k,j,mu(k,1));
    end
end

%%% POLICY ITERATION ALGORITHM


JJ=zeros(K-1,2);
UU=zeros(K-1,1);
J=zeros(K-1,L);
I=eye(K-1);




JJ(:,1)=linsolve(I-PP,GG);
%Err=1;
%Err1=1;

%checklow=0;
%for n=1:N % Iterations
a=1;
n=0;
while a==1
    
    n=n+1;
    for k=1:K-1
        for u=1:L % Computation of cost to go for each input
            CTG=zeros(L,1);
            if n>1
                JJ(:,1)=JJ(:,2);
                mu(:,1)=mu(:,2);
            end
            for j=1:K-1
                CTG(u)=CTG(u)+P(k,j,u)*JJ(j,1); % Expected cost to go
            end
            
            J(k,u)=G(k,u)+CTG(u); % Total cost for each input: E[step_cost] + E[cost to go]
            
        end
        
        [UU(k),mu(k,2)]=min(J(k,:));  % Optimization (minimization) over control inputs
        for j=1:K-1
            PP(k,j)=P(k,j,mu(k,2));
        end
        GG(k)=G(k,mu(k,2));
    end
    
    JJ(:,2)=linsolve(I-PP,GG);
    checklow=1;
    for k=1:K-1
        Err=JJ(k,2)-JJ(k,1);
        if Err>0 
           checklow=0; 
        end
    end
    if checklow==1
        for k=1:K-1
            Err1=JJ(k,2)-JJ(k,1);
            if abs(Err1)>0
                checklow=0;
            end
        end
    end
     if checklow==1 %&& Err1<tol
         a=0;
     end
%      if n==100
%          a=0;
%      end
end
J_opt(1:TERMINAL_STATE_INDEX-1,1)=JJ(1:TERMINAL_STATE_INDEX-1,2);
J_opt(TERMINAL_STATE_INDEX,1)=0;
J_opt(TERMINAL_STATE_INDEX+1:K,1)=JJ(TERMINAL_STATE_INDEX:K-1,2);

u_opt_ind(1:TERMINAL_STATE_INDEX-1,1)=mu(1:TERMINAL_STATE_INDEX-1,2);
u_opt_ind(TERMINAL_STATE_INDEX,1)=HOVER;
u_opt_ind(TERMINAL_STATE_INDEX+1:K,1)=mu(TERMINAL_STATE_INDEX:K-1,2);
end


