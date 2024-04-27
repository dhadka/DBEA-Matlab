function [f,g,x] = wfg1(objnum,x)
persistent L K;
% to keep total number of variables =24
L=[20, 20, 20, 20, 20, 14, 18, 17, 16, 15, 14, 13, 12, 11, 10 ];
K=[4, 4, 4, 4, 4, 10, 6, 7, 8, 9, 10, 11, 12, 13, 14];

if nargin == 1
	prob.nf = objnum;
	prob.ng = 0;   
    prob.nx =K(objnum)+L(objnum);
	for i = 1:prob.nx
		prob.range(i,:) = [0,2*i];
	end
	f = prob;
else
	[f,g] = wfg1_true(x,objnum,K(objnum),L(objnum));
end
return


function [f,g] = wfg1_true(x,objnum,k,l)
tp=1;
f=wfg(x, objnum, k, l, tp);
g = [];
return