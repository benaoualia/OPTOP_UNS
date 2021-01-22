clc;close all;clear all;
%-- INITIALIZE -------------------------------------------------
rmin=1.2;nelx=60;nely=60;volfrac=0.5;penal=3;
x(1:nely*nelx) = volfrac;
nodes=zeros((nelx+1)*(nely+1),2);
DDL=zeros((nelx+1)*(nely+1),2);
connective=zeros((nelx*nely),4);
coord=zeros((nelx*nely),2);
voisin=cell(nelx*nely,1);
loop = 0; 
change = 1.;
%-- MATRIX FILLING ---------------------------------------------
l=1;
for i=0:nelx
    for j=0:nely
        nodes(l,1)=i;
        nodes(l,2)=j;               % Cordinates of nodes
        l=l+1;
    end
end
for i=1:(nelx+1)*(nely+1)
    DDL(i,1)=2*i-1;                 % DDL of nodes
    DDL(i,2)=2*i;
end
l=1;
for i=1:(nelx+1)*(nely+1)
    if(not(mod(i,nely+1)==0))
        connective(l,1)=i;
        connective(l,2)=i+1;            % 4 connective of elements
        connective(l,3)=i+nely+1;
        connective(l,4)=i+nely+2;
        l=l+1;
    end
    if (l>nelx*nely) 
       break
    end
end
l=1;c=1;
for i=1:nelx*nely
    coord(i,1)=l;
    coord(i,2)=c;
    c=c+1;                              % cordinates of each element
    if (c>nely)
        l=l+1;
        c=1;
    end
end
for i = 1:nelx*nely
    for k = max(coord(i,1)-floor(rmin),1):min(coord(i,1)+floor(rmin),nelx)
        for l = max(coord(i,2)-floor(rmin),1):min(coord(i,2)+floor(rmin),nely)
            n=l+(k-1)*nely;
            if (ne(n,i)==1)                                                     % voisins
                voisin{i}=[voisin{i} n];
            end
        end
    end
end
% START ITERATION
while change > 0.01  
loop = loop + 1;
xold=x;
%-- FE analysis --------------------------------------------------
[U]=FE(nelx,nely,x,penal,connective,coord,DDL);
%-- objective function and sensitivity analysis ------------------
[KE] = lk;
  c = 0.;
for i = 1:nely*nelx
      edof=sort([DDL(connective(i,:),1) ; DDL(connective(i,:),2)]);
      Ue = U(edof,1);
      c = c + x(i)^penal*Ue'*KE*Ue;
      dc(i,1) = -penal*x(i)^(penal-1)*Ue'*KE*Ue;
end
%-- FILTERING OF SENSITIVITIES --------------------------------------------
  [dc]=check(nelx,nely,rmin,x,dc,voisin,coord);
%-- DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD -----------------------
x_oc=zeros(nely,nelx);
dc_oc=zeros(nely,nelx);
for i=1:nelx*nely
    x_oc(coord(i,1),coord(i,2))=x(i);
    dc_oc(coord(i,1),coord(i,2))=dc(i);
end
  [x]= OC(nelx,nely,x_oc,volfrac,dc_oc);
xx=x;
x=zeros(1,nelx*nely);
for i=1:nelx*nely
    x(i)=xx(coord(i,1),coord(i,2));
end
%-- PRINT RESULTS ---------------------------------------------------------
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
%-- PLOT DENSITIES --------------------------------------------------------
colormap(gray); imagesc(-xx); axis equal; axis tight; axis off;pause(1e-6);
end




















