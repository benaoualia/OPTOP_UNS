function [U]=FE(N3,C,connectiv,x,penal,F,fixeddofs,freedofs,nu,Emin,E0)
    [KE] = lk_H8(nu);
    K = sparse(N3, N3);
    U = zeros(N3,1);
    for i=1:C
        e=connectiv(i,:);
        edof = [3*e(2)-2;3*e(2)-1; 3*e(2);
                     3*e(4)-2;3*e(4)-1; 3*e(4); 
                     3*e(3)-2;3*e(3)-1; 3*e(3);
                     3*e(1)-2;3*e(1)-1; 3*e(1); 
                     3*e(6)-2;3*e(6)-1; 3*e(6); 
                     3*e(8)-2;3*e(8)-1; 3*e(8); 
                     3*e(7)-2;3*e(7)-1; 3*e(7); 
                     3*e(5)-2;3*e(5)-1; 3*e(5)];
        K(edof,edof) = K(edof,edof) +( Emin +((x(i)^penal)*(E0-Emin)))*KE;
    end
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
    U(fixeddofs,:)= 0;
end