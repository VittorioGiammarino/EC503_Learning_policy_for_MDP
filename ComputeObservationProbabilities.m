function O = ComputeObservationProbabilities(stateSpace, map)

% Function used to simulate the HMM, given a state I have a noisy sensor
% and it might measure something else

global BASE PICK_UP TREE TERMINAL_STATE_INDEX K

O=zeros(K,K);
[M,N]=size(map);
[i_base,j_base]=find(map==BASE);        % (i,j) for the base cell
[i_pickup,j_pickup]=find(map==PICK_UP); % (i,j) for pickup cell

for psi=0:1
    for i =1:M
        for j=1:N
            k=find(ismember(stateSpace, [i j psi], 'rows'));
            if map(i,j)~=TREE
                if k == TERMINAL_STATE_INDEX
                    O(k,k) = 1;
                else
                    switch j
                        case N
                          switch i
                                case M
                                    O(k,k)=0.4;
                                    if map(i-1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j-1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i j-1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i-1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                case 1
                                    O(k,k)=0.4;
                                    if map(i+1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j-1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i j-1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i+1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                              otherwise
                                    O(k,k)=0.5;
                                    if map(i-1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i+1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i+1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                           end
                        case 1
                           switch i
                                case M
                                    O(k,k)=0.4;
                                    if map(i-1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i j+1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i-1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j+1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                case 1
                                    O(k,k)=0.4;
                                    if map(i+1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j+1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i j+1 psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                                    if map(i+1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j psi], 'rows'));
                                    O(k,adj)=0.2;
                                    end
                               otherwise
                                    O(k,k)=0.5;
                                    if map(i+1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i+1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                           end
                        otherwise
                            switch i
                                case M
                                    O(k,k)=0.5;
                                    if map(i,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                case 1
                                    O(k,k)=0.5;
                                    if map(i+1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i+1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i+1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                otherwise
                                    O(k,k)=0.2;
                                    if map(i+1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j+1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j+1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i-1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i-1 j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i+1,j-1)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j-1 psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                                    if map(i+1,j)~=TREE
                                    adj = find(ismember(stateSpace, [i+1 j psi], 'rows'));
                                    O(k,adj)=0.1;
                                    end
                            end         
                    end
                end
            end
        end
    end
    k=find(ismember(stateSpace, [i_base j_base psi], 'rows'));
    O(k,k)=1;
    k=find(ismember(stateSpace, [i_pickup j_pickup psi], 'rows'));
    O(k,k)=1;
end
end