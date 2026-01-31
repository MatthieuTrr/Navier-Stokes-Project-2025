#ifndef __TEST_TRID_H__
#define __TEST_TRID_H__

#include "trid_per_C2D_0.h"
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath> // Pour fabs et calculs d'erreur
using namespace Eigen;
using namespace std;

/**
 * @brief Teste la précision et la performance du solveur tridiagonal périodique.
 * 
 * @param m Nombre de systèmes résolus simultanément (nombre de colonnes)
 * @param n Taille du système (nombre de lignes)
 * @param verbose Si true, affiche les détails du test
 * @return pair<double, double> Erreur moyenne et temps d'exécution de votre méthode
 */
pair<double, double> test_trid(int m, int n, bool verbose = true) {
    // 1. Génération de matrices de test
    // On ajoute 5.0 à la diagonale pour garantir une diagonale dominante
    MatrixXd aa = MatrixXd::Random(n, m);
    MatrixXd ac = MatrixXd::Random(n, m);
    MatrixXd ab = MatrixXd::Random(n, m).array() + 5.0; 
    MatrixXd fi = MatrixXd::Random(n, m);

    // 2. Mesure du temps pour VOTRE méthode (résolution de m systèmes)
    auto start_custom = chrono::high_resolution_clock::now();
    
    MatrixXd sol_custom = trid_per_C2D(aa, ab, ac, fi);
    
    auto end_custom = chrono::high_resolution_clock::now();
    double duration_custom = chrono::duration<double>(end_custom - start_custom).count();

    // 3. Calcul de l'erreur moyenne par rapport à la référence (Eigen LU)
    double total_error = 0.0;
    int error_count = 0;
    
    // Méthode plus robuste pour calculer l'erreur
    for (int j = 0; j < m; j++) {
        // Reconstruction de la matrice complète périodique pour la colonne j
        MatrixXd M = MatrixXd::Zero(n, n);
        for (int i = 0; i < n; i++) {
            M(i, i) = ab(i, j);                       // Diagonale
            M(i, (i + 1) % n) = ac(i, j);             // Sur-diagonale (périodique)
            M(i, (i - 1 + n) % n) = aa(i, j);         // Sous-diagonale (périodique)
        }

        // Solution de référence via LU
        VectorXd sol_ref = M.fullPivLu().solve(fi.col(j));
        
        // Calcul de l'erreur (plus robuste que la norme relative)
        double max_val_ref = sol_ref.array().abs().maxCoeff();
        double max_diff = (sol_custom.col(j) - sol_ref).array().abs().maxCoeff();
        
        if (max_val_ref > 1e-12) { // Évite division par zéro
            total_error += max_diff / max_val_ref;
            error_count++;
        } else if (max_diff > 1e-12) {
            // Solution de référence proche de zéro mais erreur significative
            total_error += max_diff;
            error_count++;
        }
        // Sinon, les deux solutions sont proches de zéro - pas d'erreur à ajouter
    }
    
    double avg_error = (error_count > 0) ? total_error / error_count : 0.0;

    // 4. Affichage si demandé
    if (verbose) {
        cout << "--- Test (n=" << n << ", m=" << m << ") ---" << endl;
        cout << "Temps execution : " << duration_custom << " s" << endl;
        cout << "Erreur relative moyenne : " << avg_error << endl;
        
        // Interprétation de l'erreur
        if (avg_error < 1e-10) {
            cout << "Precision : Excellente (erreur < 1e-10)" << endl;
        } else if (avg_error < 1e-6) {
            cout << "Precision : Tres bonne (erreur < 1e-6)" << endl;
        } else if (avg_error < 1e-3) {
            cout << "Precision : Correcte (erreur < 1e-3)" << endl;
        } else {
            cout << "Precision : Probleme (erreur > 1e-3)" << endl;
        }
        cout << "------------------------" << endl;
    }

    return make_pair(avg_error, duration_custom);
}

/**
 * @brief Teste le solveur sur une gamme de tailles et retourne les erreurs moyennes
 * 
 * @param n_min Taille minimale du système
 * @param n_max Taille maximale du système
 * @param step Pas d'incrémentation pour n
 * @param m Nombre de systèmes à résoudre simultanément
 */
void test_range(int n_min, int n_max, int step, int m = 10) {
    cout << "\n=== Test sur une gamme de tailles (m=" << m << ") ===" << endl;
    cout << "n\tErreur moyenne\tTemps (s)" << endl;
    cout << "--------------------------------" << endl;
    
    double max_error = 0.0;
    int worst_n = n_min;
    
    for (int n = n_min; n <= n_max; n += step) {
        auto result = test_trid(m, n, false);
        double avg_error = result.first;
        double duration = result.second;
        
        printf("%3d\t%.2e\t%.2e\n", n, avg_error, duration);
        
        if (avg_error > max_error) {
            max_error = avg_error;
            worst_n = n;
        }
    }
    
    cout << "\nSynthèse:" << endl;
    cout << "Erreur maximale: " << max_error << " pour n = " << worst_n << endl;
}

/**
 * @brief Teste la robustesse avec différentes valeurs de m
 * 
 * @param n Taille fixe du système
 * @param m_values Différentes valeurs de m à tester
 */
void test_different_m(int n, const vector<int>& m_values) {
    cout << "\n=== Test avec n=" << n << " fixe ===" << endl;
    cout << "m\tErreur moyenne\tTemps (s)\tTemps par système" << endl;
    cout << "------------------------------------------------" << endl;
    
    for (int m : m_values) {
        auto result = test_trid(m, n, false);
        double avg_error = result.first;
        double duration = result.second;
        double time_per_system = duration / m;
        
        printf("%3d\t%.2e\t%.2e\t%.2e\n", m, avg_error, duration, time_per_system);
    }
}

/**
 * @brief Test de validation avec des cas spécifiques
 */
void validation_tests() {
    cout << "\n=== Tests de validation ===" << endl;
    
    // Test 1: Petit système
    cout << "\n1. Petit système (n=5, m=3):" << endl;
    test_trid(3, 5, true);
    
    // Test 2: Système moyen
    cout << "\n2. Système moyen (n=50, m=10):" << endl;
    test_trid(10, 50, true);
    
    // Test 3: Grand système
    cout << "\n3. Grand système (n=200, m=20):" << endl;
    test_trid(20, 200, true);
    
    // Test 4: Beaucoup de petits systèmes
    cout << "\n4. Beaucoup de systèmes (n=10, m=100):" << endl;
    test_trid(100, 10, true);
}

#endif
