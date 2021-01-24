clear all ; clc
rmin=1.2;nelx=60;nely=60;volfrac=0.5;penal=3;
l=1;
for i=0:nelx
    for j=0:nely
        coord(l,1)=i;
        coord(l,2)=j;               % Cordinates of nodes
        l=l+1;
    end
end
l=1;
for i=1:(nelx+1)*(nely+1)
    if(not(mod(i,nely+1)==0))
        connectiv(l,1)=i;
        connectiv(l,2)=i+1;            % 4 connective of elements
        connectiv(l,3)=i+nely+1;
        connectiv(l,4)=i+nely+2;
        l=l+1;
    end
    if (l>nelx*nely) 
       break
    end
end
N=length(coord);      %nombre des noeuds
C=length(connectiv); %nombre des éléments
[voisins,centers]=calcul_voisins(coord,connectiv,rmin,N);
    F(2,1) = -1;
    fixeddofs   = union([1:2:2*(nely+1)],[2*N]);
top(coord,connectiv,volfrac,penal,rmin,nelx,nely,N,C,voisins,centers, fixeddofs,F);