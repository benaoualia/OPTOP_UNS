%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(N,C,coord,connectiv,x,penal,nelx,nely)
    [KE] = lk; 
    N2=2*N;
    K = sparse(N2, N2);
    F = sparse(N2,1);
    U = zeros(N2,1);
    for i=1:C
        e=sort(connectiv(i,:));
        edof = [2*e(1)-1; 2*e(1); 2*e(3)-1; 2*e(3);2*e(4)-1; 2*e(4); 2*e(2)-1; 2*e(2)];
        K(edof,edof) = K(edof,edof) + x(i)^penal*KE;
    end
    
    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F(2,1) = -1;
    fixeddofs   = union([1:2:2*(nely+1)],[N2]);
    alldofs     = [1:N2];
    freedofs    = setdiff(alldofs,fixeddofs);
    
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
    U(fixeddofs,:)= 0;
end
