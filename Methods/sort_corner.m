%% Corner sort
function [rank] = sort_corner(f)
nf=size(f,2);
[tmp, f_unique] = unique(f, 'rows');
dup=setdiff(1:size(f,1),f_unique)';
for i = 1:nf
    [tmp, ID{i}] = sort(f(f_unique,i));
end
f_square = f(f_unique,:) .^ 2;
for i = 1:nf
    idx = (1:nf);
    idx(i) = [];
    [tmp, ID{nf+i}] = sort(sum(f_square(:,idx),2));
end
r1 = [];
cur_id = 1;
cur_f = 1;
while length(r1) < length(f_unique)
    r = ID{cur_f}(cur_id);
    if ~ismember(r,r1)
        r1 = [r1 r];
    end
    cur_f = cur_f + 1;
    if cur_f > 2*nf
        cur_f = 1;
        cur_id = cur_id + 1;
    end
end
rank =[f_unique(r1);dup];
end

