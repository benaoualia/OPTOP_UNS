function [ ] = Display2D(x,coord,connectiv,C)
    face = [1 2 3 4];
    set(gcf,'NumberTitle','off');
    for i=1:C
         vert=[coord(connectiv(i,:),:)];
         vert(:,2,:) = -vert(:,2,:);
         patch('Faces',face,'Vertices',vert,'FaceColor',[(1-x(i)),(1-x(i)),(1-x(i))]);
         hold on;
    end
    axis equal; axis tight; axis off; pause(1e-6);
end

