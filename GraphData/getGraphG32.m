function A = getGraphG32
%Distance graph on 3 cube where xy is an edge if d_H(x,y) \leq 2. 
A = ones(8,8)-eye(8)-rot90(eye(8)); 
end

