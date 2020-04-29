function simulation(map,states,actions,stateSpace)

global TREE SHOOTER PICK_UP DROP_OFF BASE

mapsize=size(map);
[i_base,j_base]=find(map==BASE);        % (i,j) for the base cell
[i_pickup,j_pickup]=find(map==PICK_UP); % (i,j) for pickup cell
[i_DO,j_DO]=find(map==DROP_OFF);        % (i,j) for the base cell
[i_trees(:),j_trees(:)]=find(map==TREE); % (i,j) for pickup cell
[i_S(:),j_S(:)]=find(map==SHOOTER);        % (i,j) for the base cell


base = [j_base i_base];
pick_up = [j_pickup i_pickup];
trees = [j_trees' i_trees'];
shooters = [j_S' i_S'];
drop_off = [j_DO i_DO];
    
ntrees=max(length(trees));
nshooters=max(length(shooters));


if ntrees~=0
j=1;
for i=1:ntrees
    if j==nshooters+1
        j=1;
    end
figure(4)
    plot([0 mapsize(2) mapsize(2) 0 0],[0 0 mapsize(1) mapsize(1) 0],'k-',...
         [base(1)-1 base(1)-1 base(1) base(1) base(1)-1],...
         [base(2)-1 base(2) base(2) base(2)-1 base(2)-1],'k-',...
         [pick_up(1)-1 pick_up(1)-1 pick_up(1) pick_up(1) pick_up(1)-1],...
         [pick_up(2)-1 pick_up(2) pick_up(2) pick_up(2)-1 pick_up(2)-1],'k-',...
         [drop_off(1)-1 drop_off(1)-1 drop_off(1) drop_off(1) drop_off(1)-1],...
         [drop_off(2)-1 drop_off(2) drop_off(2) drop_off(2)-1 drop_off(2)-1],'k-',...
         [shooters(j,1)-1 shooters(j,1)-1 shooters(j,1) shooters(j,1) shooters(j,1)-1],...
         [shooters(j,2)-1 shooters(j,2) shooters(j,2) shooters(j,2)-1 shooters(j,2)-1],'k-',...
         [trees(i,1)-1 trees(i,1)-1 trees(i,1) trees(i,1) trees(i,1)-1],...
         [trees(i,2)-1 trees(i,2) trees(i,2) trees(i,2)-1 trees(i,2)-1],'k-')
    hold on
    fill([base(1)-1 base(1)-1 base(1) base(1) base(1)-1],...
         [base(2)-1 base(2) base(2) base(2)-1 base(2)-1],'r', ...
         [pick_up(1)-1 pick_up(1)-1 pick_up(1) pick_up(1) pick_up(1)-1],...
         [pick_up(2)-1 pick_up(2) pick_up(2) pick_up(2)-1 pick_up(2)-1],'y',...
         [drop_off(1)-1 drop_off(1)-1 drop_off(1) drop_off(1) drop_off(1)-1],...
         [drop_off(2)-1 drop_off(2) drop_off(2) drop_off(2)-1 drop_off(2)-1],'b',...
         [shooters(j,1)-1 shooters(j,1)-1 shooters(j,1) shooters(j,1) shooters(j,1)-1],...
         [shooters(j,2)-1 shooters(j,2) shooters(j,2) shooters(j,2)-1 shooters(j,2)-1],'c',...
         [trees(i,1)-1 trees(i,1)-1 trees(i,1) trees(i,1) trees(i,1)-1],...
         [trees(i,2)-1 trees(i,2) trees(i,2) trees(i,2)-1 trees(i,2)-1],'g')
     hold on
         text(base(1)-0.5, base(2)-0.5, 'B')
     hold on
         text(pick_up(1)-0.5, pick_up(2)-0.5,'P')
     hold on
         text(drop_off(1)-0.5, drop_off(2)-0.5, 'D') 
     hold on
         text(shooters(j,1)-0.5, shooters(j,2)-0.5, 'S')
     hold on
    grid on
    xlim([-2 mapsize(2)+2])
    ylim([-2 mapsize(1)+2])
    j=j+1;
end   
else
    figure(4)
    plot([0 mapsize(2) mapsize(2) 0 0],[0 0 mapsize(1) mapsize(1) 0],'k-',...
         [base(1)-1 base(1)-1 base(1) base(1) base(1)-1],...
         [base(2)-1 base(2) base(2) base(2)-1 base(2)-1],'k-',...
         [pick_up(1)-1 pick_up(1)-1 pick_up(1) pick_up(1) pick_up(1)-1],...
         [pick_up(2)-1 pick_up(2) pick_up(2) pick_up(2)-1 pick_up(2)-1],'k-',...
         [drop_off(1)-1 drop_off(1)-1 drop_off(1) drop_off(1) drop_off(1)-1],...
         [drop_off(2)-1 drop_off(2) drop_off(2) drop_off(2)-1 drop_off(2)-1],'k-',...
         [observer(1)-1 observer(1)-1 observer(1) observer(1) observer(1)-1],...
         [observer(2)-1 observer(2) observer(2) observer(2)-1 observer(2)-1],'k-')
    hold on
    fill([base(1)-1 base(1)-1 base(1) base(1) base(1)-1],...
         [base(2)-1 base(2) base(2) base(2)-1 base(2)-1],'r', ...
         [pick_up(1)-1 pick_up(1)-1 pick_up(1) pick_up(1) pick_up(1)-1],...
         [pick_up(2)-1 pick_up(2) pick_up(2) pick_up(2)-1 pick_up(2)-1],'y',...
         [drop_off(1)-1 drop_off(1)-1 drop_off(1) drop_off(1) drop_off(1)-1],...
         [drop_off(2)-1 drop_off(2) drop_off(2) drop_off(2)-1 drop_off(2)-1],'b',...
         [observer(1)-1 observer(1)-1 observer(1) observer(1) observer(1)-1],...
         [observer(2)-1 observer(2) observer(2) observer(2)-1 observer(2)-1],'m-')
     hold on
         text(base(1)-0.5, base(2)-0.5, 'B')
     hold on
         text(pick_up(1)-0.5, pick_up(2)-0.5,'P')
     hold on
         text(drop_off(1)-0.5, drop_off(2)-0.5, 'D') 
     hold on
    grid on
    xlim([-2 mapsize(2)+2])
    ylim([-2 mapsize(1)+2])
end

figure(4)
text(stateSpace(states(1),2)-0.5,stateSpace(states(1),1)-0.5,'start')
hold on
text(stateSpace(states(end),2)-0.5,stateSpace(states(end),1)-0.5,'Drop-off')

for s=1:length(states)-1
    figure(4)
    hold on
    switch actions(s)
    case 1
        txt='\uparrow';
    case 2
        txt='\downarrow';
    case 3
        txt='\rightarrow';
    case 4 
        txt='\leftarrow';
    case 5 
        txt='\times';
    end
    hold on
    text(stateSpace(states(s),2)-0.5,stateSpace(states(s),1)-0.5,txt)
    pause(1)
    
    drawnow
    F(s) = getframe(gcf);
end

video = VideoWriter('MDP_sim.avi', 'Motion JPEG AVI');
video.FrameRate = 0.7;
open(video)
writeVideo(video, F);
close(video)

end