#ifndef __TEST_TRID_H__
#define __TEST_TRID_H__

#include "trid_per_C2D_0.h"
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream> // Pour écrire le fichier de données
using namespace Eigen;
using namespace std;

/**
 * @brief Teste la précision et la performance du solveur tridiagonal périodique.
 * * @param n Taille du système (nombre de lignes)
 * @param m Nombre de systèmes résolus simultanément (nombre de colonnes)
 */
void test_trid(int m, int n) {
    // 1. Génération de matrices de test
    // On ajoute 5.0 à la diagonale pour garantir une diagonale dominante
    // Cela assure que le système possède une solution unique et stable.
    MatrixXd aa = MatrixXd::Random(n, m);
    MatrixXd ac = MatrixXd::Random(n, m);
    MatrixXd ab = MatrixXd::Random(n, m).array() + 5.0; 
    MatrixXd fi = MatrixXd::Random(n, m);

    // 2. Mesure du temps pour VOTRE méthode (résolution de m systèmes)
    auto start_custom = chrono::high_resolution_clock::now();
    
    MatrixXd sol_custom = trid_per_C2D(aa, ab, ac, fi);
    
    auto end_custom = chrono::high_resolution_clock::now();
    double duration_custom = chrono::duration<double>(end_custom - start_custom).count();

    // 3. Mesure du temps et précision pour la MÉTHODE DE RÉFÉRENCE (Eigen LU)
    double duration_ref_total = 0;
    double total_error = 0;

    cout << "--- Lancement du test (n=" << n << ", m=" << m << ") ---" << endl;

    for (int j = 0; j < m; j++) {
        // Reconstruction de la matrice complète périodique pour la colonne j
        MatrixXd M = MatrixXd::Zero(n, n);
        for (int i = 0; i < n; i++) {
            M(i, i) = ab(i, j);                       // Diagonale
            M(i, (i + 1) % n) = ac(i, j);             // Sur-diagonale (périodique)
            M(i, (i - 1 + n) % n) = aa(i, j);         // Sous-diagonale (périodique)
        }

        // Résolution via le solveur LU d'Eigen (Référence stable)
        auto start_ref = chrono::high_resolution_clock::now();
        
        VectorXd sol_eigen = M.fullPivLu().solve(fi.col(j));
        
        auto end_ref = chrono::high_resolution_clock::now();
        duration_ref_total += chrono::duration<double>(end_ref - start_ref).count();

        // Calcul de l'erreur relative sur cette colonne
        double err = (sol_custom.col(j) - sol_eigen).norm() / sol_eigen.norm();
        total_error += err;
    }

    // 4. Affichage des résultats
    cout << "Performance :" << endl;
    cout << "  > Temps Custom (total)    : " << duration_custom << " s" << endl;
    cout << "  > Temps Reference (total) : " << duration_ref_total << " s" << endl;
    
    if (duration_custom < duration_ref_total) {
        cout << "  > Gain de rapidite       : x" << duration_ref_total / duration_custom << endl;
    }

    cout << "\nPrecision :" << endl;
    cout << "  > Erreur relative moyenne : " << total_error / m << endl;
    cout << "------------------------------------------" << endl;
}

void test_vitesse_graphe() {
    // Configuration des paramètres
    std::vector<int> valeurs_m = {5, 10, 15, 20, 30};
    int n_min = 3;
    int n_max = 100;
    
    // On crée un fichier de données par valeur de m pour faciliter le tracé
    cout << "Generation des donnees de performance pour trid_per_C2D..." << endl;

    for (int m : valeurs_m) {
        std::string filename = "perf_m" + std::to_string(m) + ".dat";
        std::ofstream data_file(filename);
        data_file << "# n  temps_execution_secondes" << std::endl;

        for (int n = n_min; n <= n_max; ++n) {
            // --- Préparation des données ---
            MatrixXd aa = MatrixXd::Random(n, m);
            MatrixXd ac = MatrixXd::Random(n, m);
            MatrixXd ab = MatrixXd::Random(n, m).array() + 5.0; 
            MatrixXd fi = MatrixXd::Random(n, m);

            // --- Mesure du temps Custom uniquement ---
            auto start = chrono::high_resolution_clock::now();
            
            MatrixXd sol = trid_per_C2D(aa, ab, ac, fi);
            
            auto end = chrono::high_resolution_clock::now();
            double duration = chrono::duration<double>(end - start).count();

            data_file << n << " " << duration << std::endl;
        }
        data_file.close();
        cout << "Fichier pour m = " << m << " termine." << endl;
    }

    // --- Génération du script Gnuplot multi-courbes ---
    std::ofstream gnu_script("plot_perf.gp");
    gnu_script << "set title 'Temps d execution de trid_per_C2D selon n et m'\n";
    gnu_script << "set xlabel 'Taille du systeme (n)'\n";
    gnu_script << "set ylabel 'Temps (secondes)'\n";
    gnu_script << "set grid\n";
    gnu_script << "set term png size 800,600\n";
    gnu_script << "set output 'graphe_performance_m.png'\n";
    
    // Construction de la commande de plot dynamique
    gnu_script << "plot ";
    for (size_t i = 0; i < valeurs_m.size(); ++i) {
        int m = valeurs_m[i];
        gnu_script << "'perf_m" << m << ".dat' with lines lw 2 title 'm = " << m << "'";
        if (i < valeurs_m.size() - 1) gnu_script << ", ";
    }
    gnu_script << "\n";
    gnu_script.close();

    cout << "\nTermine ! Executez 'gnuplot plot_perf.gp' pour voir le resultat." << endl;
}

#endif