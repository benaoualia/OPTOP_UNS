clear;clc
i=-1;
j=21;
k=0;
c=1;
for n=1:(21*21*11)
    if (c==1)
        i=i+1;
    end
    if (mod(n,2)==1|| i==0)
        coord(n,:)=[i j k];
       c=1;
    else 
        coord(n,:)=[i j-1 k];
        c=c+1;
    end
    if (j==1)
        j=21;
        k=k+2;
    end
    if (i==20&&c==1)
        i=-1;
        j=j-2; 
    end
end