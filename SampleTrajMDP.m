function traj = SampleTrajMDP(P,u_opt,T)
% Simulates MDP

global TERMINAL_STATE_INDEX K
traj = {};
switch nargin
    case 2
for t=1:K
x(t,1)=t;
j=1;
while(x(t,j)~=TERMINAL_STATE_INDEX)
x_k_possible=find(P(x(t,j),:,u_opt(x(t,j)))~=0);
while(isempty(x_k_possible))
    disp('not doable policy, random action is picked')
    u_opt(x(t,j))=randi(5);
    x_k_possible=find(P(x(t,j),:,u_opt(x(t,j)))~=0);
end
prob=P(x(t,j),x_k_possible,u_opt(x(t,j)))/min(P(x(t,j),x_k_possible,u_opt(x(t,j))));
for i=2:length(prob)
    prob(i)=prob(i)+prob(i-1);
end
draw=rand/min(P(x(t,j),x_k_possible,u_opt(x(t,j))));
index_x_plus1=min(find(draw<prob));
x(t,j+1)=x_k_possible(index_x_plus1);
j=j+1;
end
traj{1,t}=x(t,find(x(t,:)~=0));
traj{2,t} = u_opt(x(t,find(x(t,:)~=0)))';
end
            
    case 3
for t=1:T
x(t,1)=randi([1,K]);
j=1;
while(x(t,j)~=TERMINAL_STATE_INDEX)
x_k_possible=find(P(x(t,j),:,u_opt(x(t,j)))~=0);
while(isempty(x_k_possible))
    disp('not doable policy, random action is picked')
    u_opt(x(t,j))=randi(5);
    x_k_possible=find(P(x(t,j),:,u_opt(x(t,j)))~=0);
end
prob=P(x(t,j),x_k_possible,u_opt(x(t,j)))/min(P(x(t,j),x_k_possible,u_opt(x(t,j))));
for i=2:length(prob)
    prob(i)=prob(i)+prob(i-1);
end
draw=rand/min(P(x(t,j),x_k_possible,u_opt(x(t,j))));
index_x_plus1=min(find(draw<prob));
x(t,j+1)=x_k_possible(index_x_plus1);
j=j+1;
end
traj{1,t}=x(t,find(x(t,:)~=0));
traj{2,t} = u_opt(x(t,find(x(t,:)~=0)))';
end
end


end