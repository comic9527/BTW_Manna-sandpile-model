function [univ,his] = whisto(A)
univ = unique(A);                       
his0 = histcounts(A,[univ, max(univ)+1]);
his = his0./sum(his0(:));
plot(univ,his,'linewidth',2);
end
