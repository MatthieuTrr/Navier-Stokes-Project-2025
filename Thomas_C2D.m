% ===============================================
% Algorithme de Thomas pour m systèmes tridiago-
% -naux utilisant le calcul vectoriel.
% ===============================================
function [x] = Thomas_C2D(aa,ab,ac,fi)
% On utilise le calcul vectoriel pour réaliser
% exactement les mêmes calculs que dans l'algo-
% -rithme de Thomas pour m sytèmes en même temps.
  [m,n]=size(fi);
  beta = zeros(m, n);
  gamma = zeros(m, n);
  x = zeros(m, n);
  beta(:,1)=ab(:,1);
  gamma(:,1)=fi(:,1)/ab(:,1);
  for i=2:n
    beta(:,i)=ab(:,i)-ac(:,i-1)./beta(:,i-1).*aa(:,i);
    gamma(:,i)=(fi(:,i)-aa(:,i).*gamma(:,i-1))./beta(:,i);
  end
  x(:,n)=gamma(:,n);
  for j=n-1:-1:1
      x(:,j) = gamma(:,j) - ac(:,j) ./ beta(:,j) .* x(:,j+1);
  end



end





