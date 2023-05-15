function [s_new, w_new] = systematic_resample(s,w,n)
    s_new=zeros(3,n);
    w_new=zeros(1,n);
    w=w';
    v = cumsum(w);
    first = rand;
    values = sort(mod((first+(1:1:size(w,1))'/size(w,1)),1));
    idx_ = arrayfun(@(i) find(values(i)<v,1,'first'), 1:size(w,1),'uniformoutput',false);
    idx = cell2mat(idx_);
    
    s_new(1,:) = s(1,idx);
    s_new(2,:) = s(2,idx);
    s_new(3,:) = s(3,idx);
    w_new = w(idx);
    w_new=w_new/sum(w_new);
    
    w_new=w_new';

end