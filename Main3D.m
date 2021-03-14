clear all ; clc; 
rmin = 1.5; %voisinnage
penal=3;%pénalisation
volfrac=0.4;
E0 = 1;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio

load('coord.mat')
load('connectiv.mat')

C=size(connectiv,1);  %nombre d'éléments
N=length(coord); %Nombre de noeuds
N3=N*3; %Nombre de ddl
alldofs     = [1:N3];
F = sparse(N3,1);
F(29,1) = -1;F(59,1) = -1;F(89,1) = -1;
fixeddofs   =[1,2,3,4,5,6,31,32,33,34,35,36,61,62,63,64,65,66];
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
