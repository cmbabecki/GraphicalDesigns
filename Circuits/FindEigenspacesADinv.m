%% Input the adjacency matrix of an undirected graph. This program computes 
%% the eigenspaces of L= D-A and separates them in the 'eigenspaces' variable. 
%% The eigenspaces are ordered from low to high frequency (outside in on the spectrum
%% [0,2maxDeg]).

function [eigenspaces, uniqueEigvals] = FindEigenspacesADinv(A)

%% computing eigenspaces
n = size(A,1);
D = zeros(n,n);
for i = 1:n
    D(i,i) = sum(A(i,:));
end
maxDeg = max(diag(D));

[U, eigvals] = eig(sym(A*inv(D)));
U = real(double(U));
eigvals = round(real(double(diag(eigvals))), 4);


%% ordering the eigenspaces

[~, I] = sort(abs(eigvals),'descend');
eigvals = eigvals(I);
U = U(:, I);

%% separating the eigenspaces
uniqueEigvals = uniquetol(eigvals);
[~, I] = sort(abs(uniqueEigvals), 'descend');
uniqueEigvals = uniqueEigvals(I);
numEigenspaces = size(uniqueEigvals, 1);
multiplicities = zeros(1, numEigenspaces);
for i = 1:size(uniqueEigvals, 1)
    multiplicities(i) = sum(eigvals == uniqueEigvals(i));
end

eigenspaces = NaN(max(multiplicities), n, numEigenspaces);
counter = 1;
for i = 1:numEigenspaces
    for j = 1: multiplicities(i)
        eigenspaces(j, :, i) = U(:, counter)';
        counter = counter + 1;
    end
end


