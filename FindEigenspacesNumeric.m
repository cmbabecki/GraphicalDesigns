%% Input a regular graph. This program computes its eigenspaces
%% and separates them in the 'eigenspaces' variable. The eigenspaces 
%% are ordered from low to high frequency.

function [eigenspaces, uniqueEigvals, multiplicities] = FindEigenspacesNumeric(A)

%% computing eigenspaces
n = size(A,1);
degree = sum(A(:,1));
[U,eigvals]=eig((A));
U = real(double(U));
eigvals= round(real(double(diag(eigvals))), 4) + .00001; %+.00001 breaks a tie by ordering the positive eigenvalue first

%% ordering the eigenspaces

[~,I]=sort(abs(eigvals),'descend');
eigvals= eigvals(I);
U= U(:,I);

%% separating the eigenspaces
uniqueEigvals=uniquetol(eigvals);
[~,I]=sort(abs(uniqueEigvals),'descend');
uniqueEigvals = uniqueEigvals(I);
numEigenspaces= size(uniqueEigvals,1);
multiplicities = zeros(1,size(uniqueEigvals,1));
for i = 1:size(uniqueEigvals,1)
    multiplicities(i) = sum(eigvals==uniqueEigvals(i));
end

eigenspaces = NaN(max(multiplicities),n,numEigenspaces);
counter =1;
for i = 1:numEigenspaces
    for j = 1: multiplicities(i)
        eigenspaces(j,:,i) = U(:,counter)';
        counter = counter +1;
    end
end
uniqueEigvals = uniqueEigvals/degree;
