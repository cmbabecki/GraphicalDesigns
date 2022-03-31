%creates adj matrix of line graph of G (d-regular, no loops or multiedges), given Adj matrix of G.  
function [N,A]=makeLineGraph(Adj)
n= size(Adj,1);
d = sum(Adj(1,:));
N=[];
for ii =1:n
    for jj = (ii+1): n
        if Adj(ii,jj) == 1
            eij = zeros(1, n);
            eij(ii) = 1;
            eij(jj) = 1;
            N=[N; eij];
        end
    end
end
N=N';
A= N'*N - 2*eye(d*n/2);
        
