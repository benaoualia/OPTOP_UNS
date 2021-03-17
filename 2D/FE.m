%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(N2,C,connectiv,x,penal,F,fixeddofs,freedofs)
    [KE] = lk;
    K = sparse(N2, N2);
    U = zeros(N2,1);
    for i=1:C
        e=connectiv(i,:);
        edof = [2*e(1)-1; 2*e(1);2*e(2)-1; 2*e(2); 2*e(3)-1; 2*e(3); 2*e(4)-1; 2*e(4)];
        K(edof,edof) = K(edof,edof) + x(i)^penal*KE;
    end
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
    U(fixeddofs,:)= 0;
end