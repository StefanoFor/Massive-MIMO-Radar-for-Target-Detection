function [MGLRT,Lambda] = MGLRT(n2,v,x,N)

Lambda = abs(v'*x)^2/n2/(x'*x);
          
MGLRT = 2*N*log(1/(1-Lambda));

end

