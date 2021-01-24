function [voisins,centers]=calcul_voisins(coord,connectiv,rmin,N)
    C=length(connectiv);    
    voisins=cell(C,1);
    centers=zeros(C,2);
    for i=1:C
        x2=max([coord(connectiv(i,1),1),coord(connectiv(i,2),1),coord(connectiv(i,3),1),coord(connectiv(i,4),1)]);
        x1=min([coord(connectiv(i,1),1),coord(connectiv(i,2),1),coord(connectiv(i,3),1),coord(connectiv(i,4),1)]);
        y2=max([coord(connectiv(i,1),2),coord(connectiv(i,2),2),coord(connectiv(i,3),2),coord(connectiv(i,4),2)]);
        y1=min([coord(connectiv(i,1),2),coord(connectiv(i,2),2),coord(connectiv(i,3),2),coord(connectiv(i,4),2)]);
        centers(i,:)=[(x2-x1)/2+x1 (y2-y1)/2+y1];
    end
    for i=1:C
        for j=1:C
            dist=sqrt((centers(i,1)-centers(j,1))^2+((centers(i,2)-centers(j,2))^2));
            if (dist<=rmin)
                voisins{i}=[voisins{i} j];
            end
        end
    end
end