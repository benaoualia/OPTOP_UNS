clear ; clc

rmin=1.2;nelx=10;nely=10;volfrac=0.5;penal=3;

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
        connectiv(l,2)=i+nely+1;            % 4 connective of elements
        connectiv(l,3)=i+nely+2;
        connectiv(l,4)=i+1;
        l=l+1;
    end
    if (l>nelx*nely) 
       break
    end
end


N=length(coord);      %nombre des noeuds
N2=2*N;


[voisins,centers,distances]=calcul_voisins(coord,connectiv,rmin);

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(N2,1);
F(2,1) = -1;
fixeddofs   = union([1:2:2*(nely+1)],[N2]);


top(coord,connectiv,volfrac,penal,rmin,N,voisins,F,fixeddofs,distances);
