function top(coord,connectiv,volfrac,penal,rmin,nelx,nely,N,C,voisins,centers)
% INITIALIZE

x(1:C) = volfrac; 
loop = 0; 
change = 1.;

% START ITERATION
while change > 0.01  
      loop = loop + 1;
      xold = x;

    % FE-ANALYSIS
      [U]=FE(N,C,coord,connectiv,x,penal,nelx,nely);  
      
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      [KE] = lk;
      c = 0.;
      for i= 1:C
            S=connectiv(i,:);
            Ue = U( [2*S(1)-1; 2*S(1); 2*S(2)-1; 2*S(2); 2*S(3)-1; 2*S(3);2*S(4)-1; 2*S(4)],1);
            c = c + x(i)^penal*Ue'*KE*Ue;
            dc(i) = -penal*x(i)^(penal-1)*Ue'*KE*Ue;
      end
      
    % FILTERING OF SENSITIVITIES
      [dc]   = check(C,rmin,x,dc,voisins,centers);    
      
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
      [x]    = OC(N,x,volfrac,dc); 
      
%     % PRINT RESULTS
      change = max(max(abs(x-xold)));
      disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
           ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(N*volfrac)) ...
            ' ch.: ' sprintf('%6.3f',change )])
        
    % PLOT DENSITIES  
    xx=reshape(x,[nelx,nely]);
      colormap(gray); imagesc(-xx); axis equal; axis tight; axis off;pause(1e-6);
end 
