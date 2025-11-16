% ===============================================
% Algorithme de Thomas pour m systèmes tridiago-
% -naux utilisant le calcul vectoriel.
% ===============================================
function [X] = Thomas_C2D(aa,ab,ac,fi)
  [m,n]=size(fi);
  beta = zeros(m, n);
  gamma = zeros(m, n);
  x = zeros(m, n);
  beta(:,1)=aa(:,1);
  gamma(:,1)=fi(:,1)/ab(:,1);


end




