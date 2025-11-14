function [X] = ThomasAlgorithm(Tridiag, f)
    n = length(f);
    X = zeros(n, 1); % Initialisation du vecteur solution
    
    % Extraction des diagonales
    A = diag(Tridiag, -1); % Sous-diagonale (longueur n-1)
    B = diag(Tridiag);     % Diagonale principale (longueur n)
    C = diag(Tridiag, 1);  % Sur-diagonale (longueur n-1)
    
    % Initialisation des vecteurs temporaires
    Beta = zeros(n, 1);
    Gamma = zeros(n, 1);
    
    % Phase descendante
    Beta(1) = B(1);
    Gamma(1) = f(1) / Beta(1);
    
    for i = 2:n
        Beta(i) = B(i) - A(i-1) * C(i-1) / Beta(i-1);
        Gamma(i) = (f(i) - A(i-1) * Gamma(i-1)) / Beta(i);
    end
    
    % Phase remontante
    X(n) = Gamma(n);
    for k = n-1:-1:1
        X(k) = Gamma(k) - C(k) * X(k+1) / Beta(k);
    end
end
