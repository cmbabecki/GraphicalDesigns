%generates adjacency matrix of expander G_p = (Z_p,E) where x is connected
%to x+1, x-1, and - x^{-1}

function  A =getGraphGp(p)
A = zeros(p,p);
A(1,1) =1;
A(1,2)=1;
A(1,p)=1;
for ii=1:p-1
    [ ~, ~,C] = gcd(p,ii);
    ModMultInv =  -1*mod(C,p);
    A(ii+1, mod(ii+1,p)+1) =1;
    A(ii+1, mod(ii-1,p)+1) =1;
    A(ii+1, mod(ModMultInv,p)+1)=  A(ii+1, mod(ModMultInv,p)+1)+1;
end
