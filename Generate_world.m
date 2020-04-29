function map = Generate_world(mapsize,ntrees,nshooters)
    if length(mapsize)~=2
        error("Mapsize wrong size, please specify the size of your 2D world")
    end
    map = zeros(mapsize(1),mapsize(2));
    trees=zeros(ntrees,2);
    for i=1:ntrees
        trees(i,1)=randi([1 mapsize(2)],1,1);
        trees(i,2)=randi([1 mapsize(1)],1,1);
        map(trees(i,2),trees(i,1))=1;
    end
    shooters=zeros(nshooters,2);
    for i=1:nshooters
        shooters(i,1)=randi([1 mapsize(2)],1,1);
        shooters(i,2)=randi([1 mapsize(1)],1,1);
        map(shooters(i,2),shooters(i,1))=2;
    end
    pick_up = [randi([1 mapsize(2)],1,1) randi([1 mapsize(1)],1,1)];
    drop_off = [randi([1 mapsize(2)],1,1) randi([1 mapsize(1)],1,1)];
    while(drop_off(1) == pick_up(1) && drop_off(2) == pick_up(2))
        drop_off = [randi([1 mapsize(2)],1,1) randi([1 mapsize(1)],1,1)];
    end
    base = [randi([1 mapsize(2)],1,1) randi([1 mapsize(1)],1,1)];
    while((base(1) == pick_up(1) && base(2) == pick_up(2)) || ...
            (base(1) == drop_off(1) && base(2) == drop_off(2)))
        base = [randi([1 mapsize(2)],1,1) randi([1 mapsize(1)],1,1)];
    end
    map(pick_up(2),pick_up(1))=3;
    map(drop_off(2),drop_off(1))=4;
    map(base(2),base(1))=5;
    
%     observer = [randi([1 mapsize(2)],1,1) randi([1 mapsize(1)],1,1)];
%     while((observer(1) == base(1) && observer(2) == base(2)))
%         observer = [randi([1 mapsize(2)],1,1) randi([1 mapsize(1)],1,1)];
%     end
if ntrees~=0
j=1;
for i=1:ntrees
    if j==nshooters+1
        j=1;
    end
figure(1)
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
    grid on
    xlim([-2 mapsize(2)+2])
    ylim([-2 mapsize(1)+2])
    j=j+1;
end   
else
    figure(1)
    plot([0 mapsize(2) mapsize(2) 0 0],[0 0 mapsize(1) mapsize(1) 0],'k-',...
         [base(1)-1 base(1)-1 base(1) base(1) base(1)-1],...
         [base(2)-1 base(2) base(2) base(2)-1 base(2)-1],'k-',...
         [pick_up(1)-1 pick_up(1)-1 pick_up(1) pick_up(1) pick_up(1)-1],...
         [pick_up(2)-1 pick_up(2) pick_up(2) pick_up(2)-1 pick_up(2)-1],'k-',...
         [drop_off(1)-1 drop_off(1)-1 drop_off(1) drop_off(1) drop_off(1)-1],...
         [drop_off(2)-1 drop_off(2) drop_off(2) drop_off(2)-1 drop_off(2)-1],'k-')
    hold on
    fill([base(1)-1 base(1)-1 base(1) base(1) base(1)-1],...
         [base(2)-1 base(2) base(2) base(2)-1 base(2)-1],'r', ...
         [pick_up(1)-1 pick_up(1)-1 pick_up(1) pick_up(1) pick_up(1)-1],...
         [pick_up(2)-1 pick_up(2) pick_up(2) pick_up(2)-1 pick_up(2)-1],'y',...
         [drop_off(1)-1 drop_off(1)-1 drop_off(1) drop_off(1) drop_off(1)-1],...
         [drop_off(2)-1 drop_off(2) drop_off(2) drop_off(2)-1 drop_off(2)-1],'b')
     hold on
         text(base(1)-0.5, base(2)-0.5, 'B')
     hold on
         text(pick_up(1)-0.5, pick_up(2)-0.5,'P')
     hold on
         text(drop_off(1)-0.5, drop_off(2)-0.5, 'D') 
    grid on
    xlim([-2 mapsize(2)+2])
    ylim([-2 mapsize(1)+2])
end

    
    
end