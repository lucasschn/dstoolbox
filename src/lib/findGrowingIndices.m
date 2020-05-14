function i_grow = findGrowingIndices(vect)
% Returns i_grow, a cell array whose cells contains sequences of
% indices at which vect grows continuously.
%
% Example 1 : i_grow = findGrowingIndices([1,2,0,0,1,2])
% i_grow = 
% 
% {[1,2]} {[4,5,6]}
%
% Example 2 : i_grow = findGrowingIndices([1,2,0,0,1,0])
% i_grow = 
% 
% {[1,2]} {[4,5]}


dv = diff(vect);
i_grow = {};
l = 1;
tmp = [];
for k=1:length(dv)
    if dv(k)>0
        tmp(end+1) = k;
        if k == length(dv)
            tmp(end+1) = k+1; 
            i_grow{l} = tmp;
        end
    elseif dv(k)<=0 && k>1 && dv(k-1)>0
        tmp(end+1) = k;
        i_grow{l} = tmp;
        tmp = [];
        l = l+1;
    end
end
end