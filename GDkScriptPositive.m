% A script to find a positively weighted k-design for each possible choice
% of k. Uncomment the graph adjacency matrix you wish to compute. Make sure
% you have working paths to Gurobi, YALMIP, and the adjacency matrix you are interested in. 
% The optimization solves a 1-norm relaxation of what is truly a 0-norm
% problem, so there is no guarantee of optimality.

clc
clear
addpath(genpath('\YALMIP-master'));
addpath(genpath('C:\gurobi911\win64\matlab'));
addpath(genpath('C:\Users\cmbab\Documents\MATLAB\GraphData')); %add your own paths
rng(1)

%% parameters + select graph
tol = 1e-5;
bf_threshold = 5e6;

% A = getGraphGraphAntilog;
% A = getGraphArayaWeiner70; 
% A = getGraphBinaryRootedTree(d); %d = depth of tree
% A = getGraphBlanusaSnark1;
% A = getGraphBlanusaSnark2;
% A = getGraphBrouwerHaemers;
% A = getGraphCelminsSwartSnark1;
% A = getGraphCelminsSwartSnark2;
% A = getGraphCn(n); %n = length of cycle
% A = getGraphQd(d); %d = dimension of cube
% A = getGraphCubeVertexTransitive66;
% A = getGraphDoubleStarSnark;
% A = getGraphFlowerSnarkJ5;  
% A = getGraphFlowerSnarkJ7;
% A = getGraphFrucht;
% A = getGraphGosset;
% A = getGraphGp(p);
% A = getGraphHadamard121;
% A = getGraphHoffmanSingleton;
% A = getGraphKneser(n,k); %Kneser graph KG(n,k)
% A = getGraphKnm(n,m); %complete bipartite graph
% A = getGraphKn(n); % complete graph
% A = getGraphLoupekinesSnark1;
% A = getGraphLoupekinesSnark2;
% A = getGraphOctohedron;
% A = getGraphPaley(p); %only works for p a prime
% A = getGraphPaley9; %several other p^t paley graphs available
% A = getGraphPetersen;
% A = getGraphRobertsonWegner;
% A = getGraphStarfish;
% A = getGraphSuetake;
% A = getGraphSuzuki;
% A = getGraphSzekeresSnark; 
% A = getGraphTruncTetraDual;
% A = getGraphTruncTetra;
% A = getGraphTruncatedCubeOcta;
% A = getGraphTutte12Cage;
% A = getGraphTutte8Cage;
% A = getGraphVanLintSchrijver;
% A = getGraphWatkinsSnark;
% A = getGraphZara;

% A = createRandRegGraph(vertices,regularity);

% [~,A] = makeLineGraph(A);  %line graph
n = size(A,1);
% A = ones(n,n) - eye(n) - A;  % complement graph

%% setup
assert(all(all(A==A')),'A is non symetrical')
assert(var(sum(A))==0,'graph is not regular')

%% get eigenspaces, eigenvalues, and their multiplicities for the graph

[eigenspaces, uniqueEigvals, multiplicities] = FindEigenspacesNumeric(A);
m = length(multiplicities); % # of distinct eigenspaces

%% initialize vars
gdk_min_sub_size = nan(1,m);
U_gdk_opt = nan(n,m);
gdk_opt = nan(n,m);
totalruntime = 0;

%% run GD k for 1 <= k < m
for k= 1:m-1 %can never integrate all eigenspaces with a proper subset
    k
        [gdk_min_sub_size(k), gdk_opt(:,k), U_gdk_opt(:,k), runtime] = findGDkPositive(eigenspaces,k);
        totalruntime = totalruntime + runtime;
        if gdk_min_sub_size(k) == n
            break
        end
end


%% Calculate efficacy of each design

efficacy = zeros(1,size(gdk_opt,2));

for j=1:size(gdk_opt,2)
    efficacy(j)= gdk_min_sub_size(:,j)/sum(multiplicities(1:j));
end


%% find an optimal design

gd_opt_index = find(efficacy == min(efficacy));
gd_opt = find (gdk_opt(:, gd_opt_index(1)))';  %records vertices of the optimal design
gd_opt_vec = gdk_opt(:,gd_opt_index); %weight vector of optimal design 
 

% there is always an extremal design with positive weights, so we do not
% compute maximal designs

%% check if 1-neighborhood of design is the whole graph

%initialize vars
gdk_min_sub_size = gdk_min_sub_size(~isnan(gdk_min_sub_size));
N=size(gdk_min_sub_size,2);
nbhd= NaN(n,N);
nbhdIsGraph = NaN(1,N);

%compute neighborhood
for i = 1:N
	nbhd(:,i) = gdk_opt(:,i);
    
    for j=1:n
        nbhd(:,i) = nbhd(:,i)+A(:,j)*gdk_opt(j,i);
    end
    %check if neighborhood is the whole graph
	if nbhd(:,i) >= ones(n,1)
        nbhdIsGraph(i) = 1;  % 1 means yes, NaN means no 
   	end
end



%% plot k by N
figure;
subplot(3,1,1)
bar([gdk_min_sub_size]')
xlabel('# eigenspaces integrated');
ylabel('minimum subset size');

