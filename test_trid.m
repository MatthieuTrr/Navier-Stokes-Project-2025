function test_trid()
    tic
    n = 30;       
    m = 10;               
    aa = rand(n, m);        
    ab = rand(n, m) + 2;     
    ac = rand(n, m);         
    fi = rand(n, m);
    fi_sol = trid_per_C2D(aa, ab, ac, fi);
    err = zeros(1, m);

    for k = 1:m
        A = diag(ab(:,k)) ...
            + diag(aa(2:end,k), -1) ...
            + diag(ac(1:end-1,k), 1);
  
        A(1,end) = aa(1,k);     
        A(end,1) = ac(end,k);   
        A(1,1) = ab(1,k);
        A(end,end) = ab(end,k);

        x_matlab = A \ fi(:,k);

        err(k) = norm( fi_sol(:,k) - x_matlab );
    end
    temps = toc;
    fprintf('----------------------------------------\n');
    fprintf('Test de la fonction trid_per_C2D\n');
    fprintf('Taille du système n = %d, nombre m = %d\n', n, m);
    fprintf('Erreur (norme) par système :\n');

    disp(err);

    fprintf('Erreur max  : %.3e\n', max(err));
    fprintf('Erreur moy : %.3e\n', mean(err));
    fprintf('----------------------------------------\n');
    fprintf('Temps écoulé : %.6f s\n', temps);
end
