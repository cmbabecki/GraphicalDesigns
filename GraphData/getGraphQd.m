%% from stack exchange 
function A = getGraphQd(d)
A=[0];
for k=0:d-1
    I=eye(2^k);
    A=[A,I;I,A];
end