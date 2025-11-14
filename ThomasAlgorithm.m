function [X] = ThomasAlgorithm(Tridiag, f)
    n = length(f);
    X = zeros(n, 1); % Initialisation du vecteur solution
    
    % Extraction des diagonales
    A = diag(Tridiag, -1); % Extraction de la sous-diagonale
    B = diag(Tridiag);     % Extraction de la diagonale principale 
    C = diag(Tridiag, 1);  % Extraction de la sur-diagonale
    
    % Initialisation des vecteurs temporaires
    Beta = zeros(n, 1);
    Gamma = zeros(n, 1);
    
    % Calcul des premiers Beta_k et Gamma_k
    Beta(1) = B(1);
    Gamma(1) = f(1) / Beta(1);
    
    for i = 2:n % Récurrence descendante
        Beta(i) = B(i) - A(i-1) * C(i-1) / Beta(i-1);
        Gamma(i) = (f(i) - A(i-1) * Gamma(i-1)) / Beta(i);
    end
    
    % Récurrence remontante
    X(n) = Gamma(n);
    for k = n-1:-1:1
        X(k) = Gamma(k) - C(k) * X(k+1) / Beta(k);
    end
end

