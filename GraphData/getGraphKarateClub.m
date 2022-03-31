function A = getGraphKarateClub

load karate2

A = sparse(edges(:,1),edges(:,2),1,34,34);
A = A+A';
