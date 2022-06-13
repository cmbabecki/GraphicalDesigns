%% Finds graphical designs through looking at circuits of the eigenvector matrix.  
%% Only use on small graphs. User selects a graph by uncommenting the desired
%% A = getGraph. User also selects operator, eigenspaces, and circuit size as prompted by program

addpath(genpath('C:\Users\cmbab\Documents\MATLAB\GraphData'));

clear
clc 
tol =  1e-05;

% A = getGraphGraphAntilog;
% A = getGraphArayaWeiner70;
% A = getGraphBlanusaSnark1;
% A = getGraphBlanusaSnark2;
% A = getGraphBrouwerHaemers;
% A = getGraphCelminsSwartSnark1;
% A = getGraphCelminsSwartSnark2;
% A = getGraphCn(7);
% A = getGraphNCube(7); 
% A = getGraphCubeVertexTransitive66;
% A = getGraphDoubleStarSnark;
% A = getGraphFlowerSnarkJ5;  
% A = getGraphFlowerSnarkJ7;
% A = getGraphFrucht;
% A = getGraphGosset;
% A = getGraphHadamard121;
% A = getGraphHarriesWong; 
% A = getGraphHoffmanSingleton;
% A = getGraphIcosahedral;
% A = getGraphKneser(5,2); %(5,2) = petersen graph
% A = getGraphKnm(n,m); %only regular if n=m
% A = getGraphKn(n);
% A = getGraphLoupekinesSnark1;
% A = getGraphLoupekinesSnark2;
% A = getGraphOctohedron;
% A = getGraphPaley(prime);
% A = getGraphPaley9; %several other p^t paley graphs available
% A = getGraphPeterson(5);
% A = getGraphRobertsonWegner;
% A = getGraphSTDPetersen;
% A = getGraphStarfish;
% A = getGraphSuetake;
% A = getGraphSuzuki;
% A = getGraphSzekeresSnark; 
% A = getGraphTruncTetraDual;
% A = getGraphTruncTetra;
% A = getGraphTruncCube;
% A = getGraphTruncatedOcta;
% A = getGraphTruncatedCubeOcta;
% A = getGraphTutte12Cage; 
% A = getGraphTutte8Cage;
% A = getGraphVanLintSchrijver;
% A = getGraphWatkinsSnark;
% A = getGraphZara;
A = getGraphSmallIrregular;

%[~,A] = makeLineGraph(A);
% A = createRandRegGraph(12,3);

n = size(A,1);

Laplacian = input("What matrix do you want to use? Press 1 for AD^{-1}, 2 for D-A:   ");
if Laplacian == 1
    [eigenspaces, eigenvalues] = FindEigenspacesADinv(A);
elseif Laplacian == 2
    [eigenspaces, eigenvalues] = FindEigenspacesDminA(A); 
end

   
U = [];
fprintf('m = # of eigenspaces = %d \n', size(eigenspaces,3))
eigenspaceIndices = input("Select eigenspaces from {2,...,m} . Enter indices as a vector []:   ");  
for i = eigenspaceIndices
   U = [U; rmmissing(eigenspaces(:,:,i))];  
end
circuitSize = input("What size circuit?   ");


%% the rows of circuitIndex records the vertices in every circuit of size circuitSize
circuitIndex = [];
ind = nchoosek(1:n,circuitSize);
for i = 1: nchoosek(n,circuitSize)
   if rank(U(:,ind(i,:))) < circuitSize
       circuitIndex =  [circuitIndex; i];
   end
end
circuitIndex =  ind(circuitIndex,:); 
numCircuits = size(rmmissing(circuitIndex),1);



circuitMatrices = zeros(size(U,1),circuitSize,numCircuits);
for i = 1: numCircuits
    circuitMatrices(:,:,i) = U(:, circuitIndex(i,:));
end

%% circuitWeights records the weights of each circuit of the chosen size
circuitWeights = NaN(circuitSize,numCircuits);
counter = 0;
for i = 1: size(circuitIndex,1)
    if size(null(circuitMatrices(:,:,i)),2) == 0
        counter = counter + 1;
    elseif size(null(circuitMatrices(:,:,i)),2) > 1
        nullSpace = null(circuitMatrices(:,:,i));
        circuitWeights(:, i) = nullSpace(:,1);
        for j = 1:circuitSize
            if abs(circuitWeights(j, i))< tol
                circuitWeights(j, i)=0;
            end
        end
    else  
        circuitWeights(:, i) = null(circuitMatrices(:,:,i));
        for j = 1:circuitSize
            if abs(circuitWeights(j, i))< tol
                circuitWeights(j, i)=0;
            end
        end
    end
end


%% posCircuitIndex records the vertices of each positive circuit
%% numPosCircuits records the number of positive circuits
%% posCircuitWeights records the weights of each positive circuit.
posCircuitWeights = NaN(circuitSize,length(circuitIndex));
nonDesigns = 0;
combCircuitWeights =NaN(circuitSize,length(circuitIndex));
for i = 1: numCircuits
    if sum(circuitWeights(:, i) > 0) == circuitSize
        posCircuitWeights(:,i) = circuitWeights(:, i);
        if posCircuitWeights(1,i) > 0
            if abs( posCircuitWeights(:,i) / posCircuitWeights(1,i) - ones(circuitSize,1)) < tol 
                combCircuitWeights(:,i) = circuitWeights(:, i);
            end
        end
    elseif sum(circuitWeights(:, i) < 0) == circuitSize
        posCircuitWeights(:,i) = -circuitWeights(:, i);
        if posCircuitWeights(1,i) > 0
            if abs( posCircuitWeights(:,i) / posCircuitWeights(1,i) - ones(circuitSize,1)) < tol 
                combCircuitWeights(:,i) = circuitWeights(:, i);
            end
        end
    elseif abs(ones(1,circuitSize) * circuitWeights(:, i))< tol
        nonDesigns = nonDesigns + 1;
    end
end
posCircuitIndex = find(~isnan(sum(posCircuitWeights,1)));
posCircuitVertices = circuitIndex(posCircuitIndex,:);
numPosCircuits = size(posCircuitIndex,2);

combCircuitIndex = find(~isnan(sum(combCircuitWeights,1)));
combCircuitVertices = circuitIndex(combCircuitIndex,:);
numCombCircuits =  sum(~isnan(sum(combCircuitWeights)));

if counter == length(circuitIndex)
    fprintf('  no circuits of this size\n')
else
    fprintf('  There are %d circuits.\n', numCircuits)
    fprintf('  There are %d positive circuits.\n', numPosCircuits)
    if numPosCircuits > 0
        fprintf('  Of the %d positive circuits, %d have equal weights.\n', numPosCircuits, numCombCircuits)
    end
    if numCircuits - numPosCircuits > 0
        fprintf('  Of the %d nonpositive circuits, %d are orthogonal to the all-ones vector. \n', numCircuits -numPosCircuits, nonDesigns)
    end
end
