%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc,voisin,coord)
dcn=zeros(nely*nelx,1);
for i=1:nelx*nely
    sum=0.0;
    for j=1:length(voisin{i,1})
        fac = rmin-sqrt((coord(i,1)-coord(voisin{i,1}(j),1))^2+(coord(i,2)-coord(voisin{i,1}(j),2))^2);
        sum = sum+max(0,fac);
        dcn(i) = dcn(i) + max(0,fac)*x(i)*dc(i);
    end
    dcn(i) = dcn(i)/(x(i)*sum);
end

% for i = 1:nelx
%   for j = 1:nely
%     sum=0.0; 
%     for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
%       for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
%         fac = rmin-sqrt((i-k)^2+(j-l)^2);
%         sum = sum+max(0,fac);
%         dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
%       end
%     end
%     dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
%   end
% end