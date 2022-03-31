function A= getGraphBinaryRootedTree(depth)
vertices = 1:(2^depth-1);
A = zeros((2^depth-1),(2^depth-1));
for ii=1: 2^(depth-1)
    A(ii, 2*ii)=1;
    A(ii, 2*ii+1)=1;
    A(2*ii, ii)=1;
    A(2*ii+1, ii)=1;
end
