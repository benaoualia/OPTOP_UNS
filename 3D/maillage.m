clear all;clc;
nely=20;nelx=60;nelz=4;
l=0;o=0;t=0;            %Initialization
for k=0:nelz
    coord([1+o:(nelx+1)*(nely+1)+o],3)=k;
    o=o+(nelx+1)*(nely+1);
    for j=0:(nelx)
       coord((1+l*(nely+1):(nely+1)+l*(nely+1)),1)=j;
       l=l+1;
       f=0;
        for i=0:nely
            coord(1+t,2)=nely-f;
            t=t+1;
            f=f+1;
        end
    end
end
j=1;k=1;
for i=1:nelx*nely*nelz
    connectiv(i,1)=j;
    connectiv(i,2)=j+1;
    connectiv(i,3)=j+nely+1;
    connectiv(i,4)=j+nely+2;
    c=(nelx+1)*(nely+1);
    connectiv(i,5)=j+c;
    connectiv(i,6)=j+1+c;
    connectiv(i,7)=j+nely+1+c;
    connectiv(i,8)=j+nely+2+c;
    j=j+1;
    if(mod(j,nely+1)==0)
        j=j+1;
    end
    if (mod(connectiv(i,4),(nelx+1)*(nely+1))==0)
        j=1+k*(nelx+1)*(nely+1);
        k=k+1;
    end
end
rmin = 1.5; %voisinnage
penal=3;%pénalisation
volfrac=0.4;
E0 = 1;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio

C=size(connectiv,1);  %nombre d'éléments
N=length(coord); %Nombre de noeuds
N3=N*3; %Nombre de ddl
alldofs     = [1:N3];
F = sparse(N3,1);
F((nelx+1)*(nely+1)*3-1,1) = -1;F((nelx+1)*(nely+1)*3*3-1,1) = -1;F((nelx+1)*(nely+1)*5*3-1,1) = -1;F((nelx+1)*(nely+1)*2*3-1,1) = -1;F((nelx+1)*(nely+1)*4*3-1,1) = -1;
fixeddofs=1;
for m=0:nelz
    fixeddofs =union(fixeddofs,[3*(1+(nelx+1)*(nely+1)*m)-2:3*((nelx+1)*(nely+1)*m+21)]);
end
freedofs    = setdiff(alldofs,fixeddofs);
x(1:C) = volfrac; 
loop = 0; 
change = 1;
 [voisins,centers,distances]=calcul_voisins(coord,connectiv,rmin,C);

 while change > 0.01  
      loop = loop + 1;
      xold = x;
      
    % FE-ANALYSIS
    [U]=FE(N3,C,connectiv,x,penal,F,fixeddofs,freedofs,nu,Emin,E0);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
     [KE] = lk_H8(nu);
      c = 0.;
      for i= 1:C
            S=connectiv(i,:);
            Ue = U(  [3*S(2)-2;3*S(2)-1; 3*S(2);
                            3*S(4)-2;3*S(4)-1; 3*S(4); 
                            3*S(3)-2;3*S(3)-1; 3*S(3);
                            3*S(1)-2;3*S(1)-1; 3*S(1); 
                            3*S(6)-2;3*S(6)-1; 3*S(6); 
                            3*S(8)-2;3*S(8)-1; 3*S(8); 
                            3*S(7)-2;3*S(7)-1; 3*S(7); 
                            3*S(5)-2;3*S(5)-1; 3*S(5)],1);
            ce=Ue'*KE*Ue;
            c = c + ( Emin +((x(i)^penal)*(E0-Emin)))*ce;
            dc(i) =-penal*(E0-Emin)*x(i)^(penal-1)*ce;
      end 
   % FILTERING OF SENSITIVITIES
      [dc]   = check(C,rmin,x,dc,voisins,distances);    
      
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
      [x]    = OC(x,volfrac,dc,C);    
      change = max(max(abs(x-xold)));
 end
 Display3D(x,coord,connectiv,C)