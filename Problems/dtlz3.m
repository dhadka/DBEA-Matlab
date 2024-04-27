function [f,g,x] = dtlz3(objnum,x)

if nargin == 1
	prob.nf = objnum;
	prob.ng = 0;
    if objnum <= 3
        prob.k=5;
    else
        prob.k = 10;
    end
	prob.nx = prob.nf+prob.k-1;
	for i = 1:prob.nx
		prob.range(i,:) = [0,1];
	end
	f = prob;
else
	[f,g] = dtlz3_true(x,objnum);
end
return
function [f,g] = dtlz3_true(x,objnum)

M = objnum;
if objnum <= 3
    k=5;
else
    k = 10;
end
temp = sum((x(:,M:M+k-1)-0.5).^2-cos(20*pi*(x(:,M:M+k-1)-0.5)),2);
G = 100*(k+temp);
temp1 = x(:,1:M-1)*pi/2;
f(:,1) = (1+G).*prod(cos(temp1),2);
for j = 2:M-1
    f(:,j)= (1+G).*prod(cos(temp1(:,1:M-j)),2).*sin(temp1(:,M-j+1));
end
f(:,M) = (1+G).*sin(x(:,1)*pi/2);
g = [];
 
 return