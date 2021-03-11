function [voisins,centers,distances]=calcul_voisins(coord,connectiv,rmin)
    C=length(connectiv);    
    voisins=cell(C,1);
    distances=cell(C,1);
    centers=zeros(C,2);
    for i=1:C
        centers(i,:)=[(coord(connectiv(i,1),1)+coord(connectiv(i,3),1))/2 (coord(connectiv(i,1),2)+coord(connectiv(i,3),2))/2];
    end
    for i=1:C
        for j=1:C
            dist=sqrt((centers(i,1)-centers(j,1))^2+((centers(i,2)-centers(j,2))^2));
            if (dist<=rmin)
                voisins{i}=[voisins{i} j];
                distances{i}=[distances{i} dist];
            end
        end
    end
end