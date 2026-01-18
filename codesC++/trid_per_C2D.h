#ifndef __TRID_PER_C2D__
#define __TRID_PER_C2D__
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
    int m=aa.rows();
    int n=aa.cols();
    //======================================================================================
    // Modification de la matrice périodique
    //======================================================================================
    ab.col(0)=ab.col(0)-aa.col(0);
    ab.col(n-1)=ab.col(n-1)-ac.col(n-1);
    MatrixXd xs2=MatrixXd::Zero(m,n);
    xs2.col(0)=aa.col(0);
    xs2.col(n-1)=ac.col(n-1);
    //======================================================================================
    // Résolution simultanée des deux systèmes tridiagonaux via l'algorithme de Thomas
    //======================================================================================
    MatrixXd alph=MatrixXd::Zero(m,n);
    alph.col(0)=RowVectorXd::Constant(m,1.0).cwiseQuotient(ab.col(0));
    for (int i=1;i<=n;i++){
        alph.col(i)=RowVectorXd::Constant(m,1.0).cwiseQuotient(ab.col(i)-aa.col(i).cwiseProduct(ac.col(i-1)).cwiseProduct(alph.col(i-1)));  
        aa.col(i)=aa.col(i).cwiseProduct(alph.col(i-1));
        ac.col(i-1)=ac.col(i-1).cwiseProduct(alph.col(i-1));
    }
    ab.col(0)=xs2.col(0);
    for (int i=1;i<=n;i++){
        ab.col(i)=xs2.col(i)-aa.col(i).cwiseProduct(ab.col(i-1));
    }
    xs2.col(n-1)=ab.col(n-1).cwiseProduct(alph.col(n-1));
    for(int i=n-2;i>=0;i--){
        xs2.col(i)=ab.col(i).cwiseProduct(alph.col(i))-ac.col(i).cwiseProduct(xs2.col(i+1));
    }
    for (int i=1;i<=n-1;i++){
        fi.col(i)=fi.col(i)-aa.col(i).cwiseProduct(fi.col(i-1));
    }
    fi.col(n-1)=fi.col(n-1).cwiseProduct(alph.col(n-1));
    for (int i=n-2;i>=0;i--){
        fi.col(i)=fi.col(i).cwiseProduct(alph.col(i))-ac.col(i).cwiseProduct(fi.col(i+1));
    }
    RowVectorXd xs1=RowVectorXd::Constant(m,0.0);
    xs1=(fi.col(0)+fi.col(n-1)).cwiseQuotient(RowVectorXd::Constant(m,1.0)+xs2.col(0)+xs2.col(n-1));

    for (int i=0; i<n-1;i++){
        fi.col(i)=fi.col(i)-xs1.cwiseProduct(xs2.col(i));
    }
    return fi;
}

#endif