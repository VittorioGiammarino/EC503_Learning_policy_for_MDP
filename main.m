%% Vittorio Gaimmarino, Boston University, 2020
clear all
close all
clc
%% Options For the world
% [M, N]
mapsize = [15, 20];

%% Global problem parameters
% Do not add or remove any global parameter 
global GAMMA R Nc P_WIND
GAMMA  = 0.2; % Shooter gamma factor
R = 2; % Shooter range
Nc = 10; % Time steps required to bring drone to base when it crashes
P_WIND = 0.1; % Gust of wind probability

% IDs of elements in the map matrix
global FREE TREE SHOOTER PICK_UP DROP_OFF BASE 
FREE = 0;
TREE = 1;
SHOOTER = 2;
PICK_UP = 3;
DROP_OFF = 4;
BASE = 5;

% Index of each action in the P and G matrices. Use this ordering
global NORTH SOUTH EAST WEST HOVER
NORTH  = 1;
SOUTH = 2;
EAST = 3;
WEST = 4;
HOVER = 5;

%% Generate map
% map(m,n) represents the cell at indices (n,m) according to the axes
% specified in the PDF.
[map] = Generate_world(mapsize,20,3);

%% Generate state space
disp('Generate state space');
% Generate a (K x 3)-matrix 'stateSpace', where each accessible cell is
% represented by two rows (with and without carrying a package).
stateSpace = [];
for m = 1 : size(map, 1)
    for n = 1 : size(map, 2)
        if map(m, n) ~= TREE
            stateSpace = [stateSpace;
                          m, n, 0;
                          m, n, 1];
        end
    end
end
% State space size
global K
K=size(stateSpace,1);
global TERMINAL_STATE_INDEX
TERMINAL_STATE_INDEX = ComputeTerminalStateIndex(stateSpace, map); % Compute terminal state
P = ComputeTransitionProbabilities(stateSpace, map); % Compute Transition Probability
O = ComputeObservationProbabilities(stateSpace, map); % Compute observation Probability if in HMM
G = ComputeStageCosts(stateSpace, map); % Compute cost for each cell in the map

disp('Solve stochastic shortest path problem with Policy Iteration');    
[ J_opt_pi, u_opt_ind_pi ] = PolicyIteration(P, G);   % Solve the problem using Policy iteration
if size(J_opt_pi,1)~=K || size(u_opt_ind_pi,1)~=K
   disp('[ERROR] the size of J and u must be K')
end
disp('Done'); 
%% Plot Optimal solution
plotOptimalSolution(map,stateSpace,u_opt_ind_pi) % Optimal solution for the expert


%% HMM 
% See generate HMM samples

% obs=100;
% 
% [seq,states] = hmmgenerate(obs,P,O);
% [estimateTR,estimateE] = hmmestimate(seq,states);
% % 
% % %%
% simulation(map,states,stateSpace);
% % 
% % %%
% figure()
% plot(stateSpace(states(:),2)-0.5,stateSpace(states(:),1)-0.5,'rx','LineWidth',10)

%% Sampling from optimal solition to learn Sarsa
T=100; %Number of trajectories

traj = SampleTrajMDP(P,u_opt_ind_pi,T);

%simulation(map,traj{1,6},stateSpace); % Uncomment to see simulation of a
%trajectory
Xtr=[];
for t=1:length(traj)
    Xtr=[Xtr; traj{1,t}' traj{2,t}'];
end

%% SARSA
actions=5;
%initQ=randi([1 actions],[K actions]);
initQ=ones(K, actions);
epsilon=0.1;
gamma=0.75;
alpha=0.25;
T=10000;
steps=1000;
disp('running SARSA')
[Q,reward_Sarsa] = SARSA(map,stateSpace,P,initQ,epsilon,gamma,alpha,T,steps);
temp_Q=Q';
[~,u] = (max(temp_Q));
u_Sarsa=u';
disp('done')
figure(5)
plot(1:1:T,reward_Sarsa)
ylabel('reward')
xlabel('iter')
ylim([-30 110])
xlim([-1 T+10])
title('cumulative reward SARSA')

close 2 3
plotOptimalSolution(map,stateSpace,u_Sarsa);


%% SARSA from experts
initQ=ones(K, actions);

for i=1:length(Xtr)
    initQ(Xtr(i,1),Xtr(i,2))=100;
end
T=1000;
epsilon=0.01;
steps=500;
disp('running SARSA')
[Q,reward_Sarsa_exp] = SARSA(map,stateSpace,P,initQ,epsilon,gamma,alpha,T,steps);
temp_Q=Q';
[~,u] = (max(temp_Q));
u_Sarsa_exp=u';
disp('done')
figure(6)
plot(1:1:T,reward_Sarsa_exp)
ylabel('reward')
xlabel('iter')
ylim([-30 110])
xlim([-1 T+10])
title('cumulative reward SARSA from expert')

close 2 3
plotOptimalSolution(map,stateSpace,u_Sarsa_exp);

%% Sampling from optimal solition to learn inverse RL
T=100; %Number of trajectories

traj = SampleTrajMDP(P,u_opt_ind_pi);

%simulation(map,traj{1,6},stateSpace); % Uncomment to see simulation of a
%trajectory
Xtr=[];
for t=1:length(traj)
    Xtr=[Xtr; traj{1,t}' traj{2,t}'];
end

%% Learning Policies for MDP logistic regression
% Define Waypoints needed for the features

%define waypoints
waypoints=[];
waypoints(1,1)=ComputePickUpStateIndex(stateSpace, map);
waypoints(2,1)=TERMINAL_STATE_INDEX;
for i=2:2:mapsize(2)
    for j=2:2:mapsize(1)
        if map(j,i)~=TREE
            for k=1:K
                if (stateSpace(k,1)==j && stateSpace(k,2)==i)
                    waypoints = [waypoints; k]; 
                end
            end
        end
    end
end

%% Compute Feature(samples,actions) 
% Feature(samples,actions*waypoints) is needed in the Logistic regression

c=0.0001; %tune this parameter
y_hat=waypoints; % encoding for waypoints
y = DataEncoding(Xtr,P); % encoding for the rest of the points in the training set
Gaussian_kernel = @(x,y) exp(-c*norm(x-y)^2); %Kernel used as feature, other can be tried

Data_P_actions=[];
for t=1:length(Xtr)
    actions=[];
    actions=[Xtr(t,1) NORTH; Xtr(t,1) SOUTH; Xtr(t,1) EAST; Xtr(t,1) WEST; Xtr(t,1) HOVER];
    Data_P_actions(:,t) = DataEncoding(actions,P); % computes encoding for each combination sample action
end
actions = 5;
Features_sample_action=[];
for t=1:length(Xtr)
    for i=1:length(y_hat)
        for a=1:actions
               Features_sample_action(t,a +(i-1)*actions)=Gaussian_kernel(y_hat(i,1),Data_P_actions(a,t)); %compute feature(samples,actions*waypoints)
        end
    end
end


%% Logistic Regression

number_samples=size(Xtr,1);

options = optimoptions(@fminunc,'Algorithm','quasi-newton',... %Try different algorithms
  'HessUpdate','bfgs',...
  'MaxFunEvals',1e14,...
  'TolX',1e-20,...
  'TolFun',1e-20,...
  'MaxIter',1e14,...
  'Display','off',...
  'GradObj','on',...
  'Hessian','off',...
  'DerivativeCheck','off');

r0 = zeros(length(y_hat),1); % weights initialization
n = 40;
lambda = logspace(-8, -4, n);

weights = zeros(length(y_hat), n);
log_likelihood = zeros(1, n);
validation = zeros(1, n);

ratio = 0.7; %ratio samples testing/validation
num_train = ceil(ratio*number_samples);
RandID = randperm(number_samples);
samples_control=Xtr(:,2);

Train_samples_feature_selected = Features_sample_action(RandID(1:num_train),:);
Train_samples_control = samples_control(RandID(1:num_train));
Test_samples_feature_selected = Features_sample_action(RandID(num_train+1:end),:);
Test_samples_control = samples_control(RandID(num_train+1:end));


parfor i = 1:n
  disp(i)
  [coeffs,fval] = fminunc(...
    @(r)logisticRegLikelihood(r, Train_samples_feature_selected, Train_samples_control, actions, lambda(i))...
    ,r0,options);
  weights(:, i) = coeffs;
  log_likelihood(i) = fval;
  validation(i) = logisticRegLikelihood(coeffs, Test_samples_feature_selected, Test_samples_control, actions, 0);
end

[~,bestid] = min(validation);

figure()
plot(lambda,validation)
ylabel('logisticRegLikelihood')
xlabel('\lambda')
title('validation')
%% Build policy

Final_weights=weights(:,bestid);

Final_policy_actions=[];
for t=1:K
    actions=[];
    actions=[t NORTH; t SOUTH; t EAST; t WEST; t HOVER];
    Final_policy_actions(:,t) = DataEncoding(actions,P);
end
actions = 5;
Final_features_sample_action=[];
for t=1:K
    for i=1:length(y_hat)
        for a=1:actions
               Final_features_sample_action(t,a +(i-1)*actions)=Gaussian_kernel(y_hat(i,1),Final_policy_actions(a,t));
        end
    end
end

W=kron(eye(actions),Final_weights);
XWc=Final_features_sample_action*W;
XWc = exp(XWc);
for i=1:K
u_stoch(i,:) = XWc(i,:)/sum(XWc(i,:)); %Stochastic policy
[~,I]=max(u_stoch(i,:));
u_learnt(i,1)=I; %Deterministic policy
end

error=numel(find(u_opt_ind_pi~=u_learnt))/numel(u_opt_ind_pi);

%% Test policy
close 2 3
plotOptimalSolution(map,stateSpace,u_learnt)
%% Code for video simulation
% T=1; %Number of trajectories
% traj_u = SampleTrajMDP_stochastic(P,u_stoch,T,41);
% simulation(map,traj_u{1,1},traj_u{2,1},stateSpace);




