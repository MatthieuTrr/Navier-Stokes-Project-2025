#ifndef __TRID_PER_C2D_0__
#define __TRID_PER_C2D_0__
#include <iostream>
using namespace std;
#include <Eigen/Dense>
using namespace Eigen;

// Fonction pour tester la compilation sans tester le code
void parler(){
    cout<<"Le code compile"<<endl;
    MatrixXd aa = MatrixXd::Random(3, 3);
    cout<<"("<<aa.rows()<<","<<aa.cols()<<")"<<endl;
}


MatrixXd trid_per_C2D(MatrixXd aa, MatrixXd ab, MatrixXd ac, MatrixXd fi){
    //======================================================================================
    // Modification de la matrice périodique
    //======================================================================================
    int n=aa.rows();
    int m=aa.cols();
    MatrixXd v=MatrixXd::Zero(n,m);
    v.row(0)=aa.row(0);
    v.row(n-1)=ac.row(n-1);
    MatrixXd aa_star=aa.block(1,0,n-1,aa.cols());
    MatrixXd ac_star=ac.block(0,0,n-2,ac.cols());
    MatrixXd ab_star=ab;
    ab_star.row(0)=ab.row(0)-aa.row(0);
    ab_star.row(n-1)=ab.row(n-1)-ac.row(n-1);
    //======================================================================================
    // Résolution simultanée des deux systèmes tridiagonaux via l'algorithme de Thomas
    //======================================================================================
    // Calcul des betas et gammas
    MatrixXd beta = MatrixXd::Zero(n,m), 
         gamma1 = MatrixXd::Zero(n,m),
         gamma2 = MatrixXd::Zero(n,m),
         X1 = MatrixXd::Zero(n,m),
         X2 = MatrixXd::Zero(n,m),
         sol = MatrixXd::Zero(n,m);
    // beta(1,:)=ab(1,:);
    beta.row(0)=ab.row(0);
    // gamma(1,:)=fi(1,:)./ab(1,:);
    gamma1.row(0)=fi.row(0).cwiseQuotient(ab.row(0));
    gamma2.row(0)=v.row(0).cwiseQuotient(ab.row(0));
    for (int i=1; i<n; i++){
        // beta(i,:)=ab(i,:)-ac(i-1,:).*aa(i-1,:)./beta(i-1,:);
        beta.row(i)=ab.row(i)-ac.row(i-1).cwiseProduct(aa.row(i-1)).cwiseQuotient(beta.row(i-1));
        // gamma(i,:)=(fi(i,:)-aa(i-1,:).*gamma(i-1,:))./beta(i,:);
        gamma1.row(i)=(fi.row(i)-aa.row(i-1).cwiseProduct(gamma1.row(i-1))).cwiseQuotient(beta.row(i));
        gamma2.row(i)=(v.row(i)-aa.row(i-1).cwiseProduct(gamma2.row(i-1))).cwiseQuotient(beta.row(i));
    }

    //Construction des X1 et X2 par récurrence inverse
    X1.row(n-1)=gamma1.row(n-1);
    X2.row(n-1)=gamma2.row(n-1);
    for (int j=n-2; j>=0; j--){
        // x(j,:) = gamma(j,:) - ac(j,:).* x(j+1,:) ./ beta(j,:) ;
        X1.row(j)=gamma1.row(j)-ac.row(j).cwiseProduct(X1.row(j+1)).cwiseQuotient(beta.row(j));
        X2.row(j)=gamma2.row(j)-ac.row(j).cwiseProduct(X2.row(j+1)).cwiseQuotient(beta.row(j));
    }
    //======================================================================================
    // construction de la solution finale à partir de X1 et X2
    //======================================================================================
    RowVectorXd X_star(m);
    // X_star=(X1(1,:)+X1(n,:))/(1+X2(1,:)+X2(n,:));
    X_star=(X1.row(0)+X1.row(n-1)).cwiseQuotient(X2.row(0)+X2.row(n-1)+RowVectorXd::Constant(m,1.0));
    sol=X1-X2*X_star.asDiagonal();
    return sol;
}

#endif