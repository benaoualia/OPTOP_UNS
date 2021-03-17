function [ ] = Display3D(x,coord,connectiv,C)
    face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
    set(gcf,'Name','ISO display','NumberTitle','off');
    for i=1:C
         if (x(i) >= 0.5) 
             vert=[coord(connectiv(i,:),:)];
             vert(:,[2 3]) = vert(:,[3 2]);
             vert(:,2,:) = -vert(:,2,:);
             vert([3 4],:)=vert([4 3],:);
             vert([7 8],:)=vert([8 7],:);
             patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-x(i)),0.2+0.8*(1-x(i)),0.2+0.8*(1-x(i))]);
             hold on;
         end
    end
    axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end

