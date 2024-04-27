function [f,g,x] = zdt3(x)
if nargin == 0
	prob.nx = 10;
	prob.nf = 2;
	prob.ng = 0;
	prob.range = cell(1, prob.nx);
	for i = 1:prob.nx
		prob.range{i} = Range('range', [0.0, 1.0]);
	end
	f = prob;
else
	[f,g] = zdt3_true(x);
end
return


function [f,g] = zdt3_true(x)
N = length(x);
f(1) = x(1);
g = 1 + 9/(N-1) * sum(x(2:N));
h = 1 - sqrt(f(1)/g) - (f(1)/g)*sin(10*pi*f(1));
f(2) = g*h+1;

g = [];
return;
