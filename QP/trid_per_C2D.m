% ===============================================
% fonction utilisant la transformation de matrices
% périodiques et l'algorithme de Thomas pour ré-
% -soudre m systèmes linéaires à matrices périodi-
% -ques.
% ===============================================
function [fi]=trid_per_C2D(aa,ab,ac,fi)
 % On souhaite résoudre m systèmes linéaires périodiques.
 % L'idée est de "stocker" nos systèmes tridiagonaux dans
 % 3 matrices: aa, ab, ac; contenant respectivement les 3
 % diagonales de tous les systèmes en colonnes. aa est donc
 % une matrice de taille m*n dont les colonnes sont les sous-
 % diagonales de chaque système.
[n, m] = size(ab);
v=zeros(n,m);
v(1,:)=aa(1,:);
v(n,:)=ac(n,:);
aa_star=aa(2:n,:);
ac_star=ac(1:n-1,:);
ab_star = ab;
ab_star(1,:)=ab(1,:)-aa(1,:);
ab_star(n,:)=ab(n,:)-ac(n,:);
X1=Thomas_C2D(aa_star,ab_star,ac_star,fi);
X2=Thomas_C2D(aa_star,ab_star,ac_star,v);
X_star=(X1(1,:)+X1(n,:))/(1+X2(1,:)+X2(n,:));
sol=X1-X2.*X_star;
fi=sol;
end
