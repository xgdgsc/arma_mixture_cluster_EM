
function lse = logSumExp ( X )

% function lse = logSumExp ( X )
% logSumExp - Log ( Sum ( Exponential ) )
%
% usage
% lse = logSumExp ( X )
%
% input
% 	X : (M x N)-matrix of column vectors ( with log values )
%
% output
% 	lse : The log sum exponential over all columns.
%
% description
% 	Computes the log sum exponential for some vector of values. It is used
% 	to avoid numerical underflow when calculcating the log of a sum of small 
%		numbers like probabilities. 
%
% 	In a mixture model the probability of an event x is 
% 	P(x) = pi1*P1(x) + pi2*P2(x)... 
% 	The problem is usually that  P1, P2,... are small at x which makes 
% 	underflow happen. To fix underflow one usually operates in the log 
% 	domain i.e. log(P(x)) = log(pi1*P1(x) + pi2*P2(x)...)
% 	The problem with this is that the log cannot decompose sums and we 
% 	still get underflow. To fix this(somewhat) we can write:
% 	log(P(x)) = log( exp(log(pi1) + log(P1(x)) ) + 
%		log( exp(log(pi2) + log(P2(x)) ) + ...). Now by finding the max value
%		of pi1*P1(x),pi2*P2(x),... and deducting it we can remove most of the
%		value in the equation and get it out. It is simple if one looks at the 
%		following calculations
%
%		log(p) = log(p1 + p2) = log(exp(log(p1))+exp(log(p2)))
%		pMax = max([log(p1),log(p2)])
%		log(p) = log( exp(pMax) * ( exp(log(p1)-pMax)+exp(log(p2-pMax)) )
%		log(p) = pMax + log( exp(log(p1)-pMax) + exp(log(p2)-pMax) )
% 	
%		Now if we for example assume log(p1)>log(p2) then 
%		log(p) = pMax + log( exp(0) + exp(log(p2)-pMax) ) =
%		pMax + log( 1 + exp(log(p2)-pMax) )
%		
% 	This means that we gotten out most of the probability mass from the 
%		sum and we have avoided summing several small numbers. Hopefully 
%		the exp(log(p2)-pMax) will be nice as well. 
% 	
% 	
%
% author
%     Martin Hjelm, mar.hjelm@gmail.com


%%%%%%%%%% CHECK INPUT ETC. %%%%%%%%%%%%

% Check erroneous input
  if nargin < 1
      error('PCA.m: Too few input arguments. For help type help PCA.\n');    
  end  
  
 	[M,~] = size(X);
 	% Get max in all columns
 	xMax = max ( X );
 	% Subtract max in all columns from all column values and 
 	% calculate the final log sum exp
 	lse = xMax + log ( sum ( exp( X - ones(M,1) * xMax ) ) );
    %上式中，将X的每一列都减去这一列中的最大值，然后对每个数取指数，然后对每列取和，再取对数，然后再加上最大值
    %在仿真1中，lse最终的大小为1×30

 	% Check for infinity and if so take xMax as the value
  lseInf = ~isfinite(lse);
  lse(lseInf) = xMax(lseInf); 

end