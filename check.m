%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(C,rmin,x,dc,voisins,centers)
dcn=zeros(C,1);
for i = 1:C
    sum=0.0; 
    voisin=voisins{i};
    m=length(voisin);
    for j=1:m
        fac = rmin-sqrt((centers(i,1)-centers(voisin(j),1))^2+((centers(i,2)-centers(voisin(j),2))^2));
        sum = sum+fac;
        dcn(i) = dcn(i) +fac*x(j).*dc(j);
      end
    end
    dcn(i) = dcn(i)/(x(i)*sum);
end
