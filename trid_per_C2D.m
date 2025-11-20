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
v = zeros(n, m);
sol = zeros (m);
for j = 1:m
    Xk = zeros(n);
    a = aa(:, j);
    b = ab(:, j);
    c = ac(:, j);
    b(1) = b(1) - a(1);   
    b(n) = b(n) - c(n);   
    v(1, j) = a(1);
    v(n, j) = c(n);
    X1 = Thomas_C2D(a,b,c,fi(j));
    X2 = Thomas_C2D(a,b,c,v);
    Xstar = (X1(1) + X1(n))/(1+X2(1)+X2(n));
    Xk = X1-(X2.*Xstar);
    sol(:, j) = Xk; 
end
fi = sol; 
end