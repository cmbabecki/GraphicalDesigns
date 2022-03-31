%generates adjacency matrix of the graph G_i defined in my paper for the
%Johnson scheme

function A = getGraphJG(j,n,k)
ind = nchoosek([1:n], k);
indicatorvecs = zeros(size(ind,1),n);
for ii = 1:size(ind,1)
    indicatorvecs(ii,ind(ii,:))= 1;
end
A= zeros(nchoosek(n,k),nchoosek(n,k));
for ii =1:nchoosek(n,k)
    for jj =1:nchoosek(n,k)
        if indicatorvecs(ii,:)* indicatorvecs(jj,:)'==k-j
            A(ii,jj) = 1;
        end
    end
end
end