% integer program for finding a minimum cardinality combinatorial k-graphical design
% inputs: 
%   eigenspaces = eigenspaces of a graph ordered by frequency (from FindEigenspacesNumeric) 
%   k = number of eigenspaces to average 
% outputs:
%   min_sub_size = size of a minimum k-design (combinatorial)
%   x_opt =  indicator vector of this minimum k-design
%   Ux_opt = is 0 on rows corresponding to eigenvectors which x_opt averages
%   runtime = tracks run time of program

function [min_sub_size, x_opt, Ux_opt, runtime] = findGDk(eigenspaces,k)

opts = sdpsettings;
opts.verbose = 1;
tol = 1e-8;
n = size(eigenspaces,2);

%constructs matrix of eigenspaces 2-k
Uk = [];
for i=2:k
    Uk = [Uk; rmmissing(eigenspaces(:,:,i))];
end

%constructs full eigenbasis matrix
U = [ones(1,n);Uk];
for i= (k+1): size(eigenspaces,3)
    U = [U; rmmissing(eigenspaces(:,:,i))];
end

%setting up integer program
t_start = tic;
x = binvar(n,1);

A= [Uk; -Uk; -ones(1,n)];
b = [zeros(2*(size(Uk,1)),1); -1 ];

optimize(A*x <=b , ones(1,n)*x,opts)

% computes outputs
runtime = toc(t_start);

x_opt = round(value(x));
Ux_opt = U*(x_opt/(ones(1,n)*x_opt)-1/n);
min_sub_size = ones(1,n)*x_opt;
