function A = getGraphCn(n)
c = zeros(1,n);
c(1,2)= 1;
c(1,n) =1;
A = zeros(n,n);
for i=0:(n-1)
    A(i+1,:) = circshift(c,i);
end

    

%easy to make any circulant matrix this way