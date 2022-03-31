function A= getGraphPaley(q)
% Create a Paley graph with q vertices
% Only enter q a prime power equiv 1 mod 4.

squares = unique(mod([1:q].^2,q));
squares = squares(2:end);

A= zeros(q,q);
for ii=1:q
    for jj=(ii+1):q
        if ismember(ii-jj, squares) == 1 || ismember(jj-ii, squares) ==1
            A(ii,jj)=1;
            A(jj,ii)=1;
        end
    end
end
            


