function traj = SampleTrajMDP_stochastic(P,u_opt,T,initial_state)
% Simulate MDP for Randomized Stationary Policy

global TERMINAL_STATE_INDEX K
traj = {};
switch nargin
    case 2
for t=1:K
x(t,1)=t;
j=1;
while(x(t,j)~=TERMINAL_STATE_INDEX)
    
%select action
prob_u=u_opt(x(t,j),:)/min(u_opt(x(t,j),:));
for i=2:length(prob_u)
    prob_u(i)=prob_u(i)+prob_u(i-1);
end
draw_u=rand/min(u_opt(x(t,j),:));
u=min(find(draw_u<prob_u));
u_final(t,j)=u;

%next state
x_k_possible=find(P(x(t,j),:,u)~=0);
while(isempty(x_k_possible))
    disp('not doable policy, random action is picked')
    u=randi(5);
    x_k_possible=find(P(x(t,j),:,u)~=0);
end
prob=P(x(t,j),x_k_possible,u)/min(P(x(t,j),x_k_possible,u));
for i=2:length(prob)
    prob(i)=prob(i)+prob(i-1);
end
draw=rand/min(P(x(t,j),x_k_possible,u));
index_x_plus1=min(find(draw<prob));
x(t,j+1)=x_k_possible(index_x_plus1);
j=j+1;
end
traj{1,t}=x(t,find(x(t,:)~=0));
traj{2,t} = u_final;
end
            
    case 3
for t=1:T
x(t,1)=randi([1,K]);
j=1;

while(x(t,j)~=TERMINAL_STATE_INDEX)
    
%select action
prob_u=u_opt(x(t,j),:)/min(u_opt(x(t,j),:));
for i=2:length(prob_u)
    prob_u(i)=prob_u(i)+prob_u(i-1);
end
draw_u=rand/min(u_opt(x(t,j),:));
u=min(find(draw_u<prob_u));
u_final(t,j)=u;
%next state
x_k_possible=find(P(x(t,j),:,u)~=0);
while(isempty(x_k_possible))
    disp('not doable policy, random action is picked')
    u=randi(5);
    x_k_possible=find(P(x(t,j),:,u)~=0);
end
prob=P(x(t,j),x_k_possible,u)/min(P(x(t,j),x_k_possible,u));
for i=2:length(prob)
    prob(i)=prob(i)+prob(i-1);
end
draw=rand/min(P(x(t,j),x_k_possible,u));
index_x_plus1=min(find(draw<prob));
x(t,j+1)=x_k_possible(index_x_plus1);

j=j+1;
end
traj{1,t}=x(t,find(x(t,:)~=0));
traj{2,t} = u_final;
end

    case 4
t=1;
x(t,1)=initial_state;
j=1;

while(x(t,j)~=TERMINAL_STATE_INDEX)
    
%select action
prob_u=u_opt(x(t,j),:)/min(u_opt(x(t,j),:));
for i=2:length(prob_u)
    prob_u(i)=prob_u(i)+prob_u(i-1);
end
draw_u=rand/min(u_opt(x(t,j),:));
u=min(find(draw_u<prob_u));
u_final(t,j)=u;
%next state
x_k_possible=find(P(x(t,j),:,u)~=0);
while(isempty(x_k_possible))
    disp('not doable policy, random action is picked')
    u=randi(5);
    x_k_possible=find(P(x(t,j),:,u)~=0);
end
prob=P(x(t,j),x_k_possible,u)/min(P(x(t,j),x_k_possible,u));
for i=2:length(prob)
    prob(i)=prob(i)+prob(i-1);
end
draw=rand/min(P(x(t,j),x_k_possible,u));
index_x_plus1=min(find(draw<prob));
x(t,j+1)=x_k_possible(index_x_plus1);

j=j+1;
end
traj{1,t}=x(t,find(x(t,:)~=0));
traj{2,t} = u_final;
end


end