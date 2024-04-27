%% nd_sort() - Find non-dominated points from a given set
% dir = 1, minimize all values
% dir = 2, maximimze all values
function [ndset,idx] = nd_rank1(set1,dir)
if nargin == 1
	dir = 1;
end
[N,M] = size(set1);
switch(dir)
	case 1
		dom = nd_sort_min(set1, M, N);
	case 2
		dom = nd_sort_max(set1, M, N);
	otherwise
		error('wrong value of dir');
end

idx = [];
for i = 1:N
	if dom(i) == 0, idx = [idx i]; end
end
ndset = set1(idx(:),:);
return


%% nd_sort() for minimization
function [dom] = nd_sort_min(set1, M, N)
dom = zeros(1, N);
for i = 1:N
	f = repmat(set1(i,:), N, 1);
	dom_less = sum(f <= set1, 2);
	dom_more = sum(f >= set1, 2);
	for j = 1:N
		if dom_less(j) == M && dom_more(j) < M
			dom(j) = dom(j) + 1;
		end
	end
end
return


%% nd_sort() for maximization
function [dom] = nd_sort_max(set1, M, N)
dom = zeros(1, N);
for i = 1:N
	f = repmat(set1(i,:), N, 1);
	dom_less = sum(f <= set1, 2);
	dom_more = sum(f >= set1, 2);
	for j = 1:N
		if dom_more(j) == M && dom_less(j) < M
			dom(j) = dom(j) + 1;
		end
	end
end
