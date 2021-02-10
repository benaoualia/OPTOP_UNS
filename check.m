%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(C,rmin,x,dc,voisins,distances)
dcn=zeros(C,1);
for i = 1:C
    sum=0.0; 
    voisin=voisins{i};
    distance=distances{i};
    m=length(voisin);
    for j=1:m
        fac = rmin-distance(j);
        sum = sum+fac;
        dcn(i) = dcn(i) +fac*x(voisin(j)).*dc(voisin(j));
    end
    dcn(i) = dcn(i)/(x(i)*sum);
end
end
