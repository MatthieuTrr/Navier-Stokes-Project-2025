#ifndef __CALC_LAP_H__
#define __CALC_LAP_H__
#include <Eigen/Dense>
using namespace Eigen;

float Lx=1;
float Ly=2;
int nx=21;
int ny=51;
int cfl=1;
float dx=Lx/(nx-1.0);
float dy=Ly/(ny-1.0);

// Vecteurs servant à la discrétisation des opérateurs différentiels
VectorXi ip= VectorXi::Linspaced(nx-1,1,nx-1); 
ip(nx-2)=0;
VectorXi jp= VectorXi::Linspaced(ny-1,1,ny-1); 
jp(ny-2)=0;
VectorXi im= VectorXi::Linspaced(nx-1,-1,nx-2); 
im(0)=nx-2;
VectorXi jm= VectorXi::Linspaced(ny-1,-1,ny-2); 
jm(0)=ny-2;
MatrixXd calc_Lap(MatrixXd u){ // u est un tableau (= matrice) de taille nx-1 * ny-1
    // Pour calculer une valeur approchée du Laplacien, on utilise la formule énoncée dans le sujet:
    // hc = ( u(ip,jc) - 2 * u + u(im,jc) ) / (dx * dx) + ( u(ic,jp) - 2 * u + u(ic,jm) ) / (dy * dy) 
    // Pour ordonner les lignes et colonnes de u selon ip,im,jp,jm; les coordonnées des vecteurs doivent être entre 0 et n-1 (pour une matrice n*n)
    // On doit donc faire ip-1,im-1,... pour pouvoir modifier u correctement.
    PermutationMatrix<Dynamic,Dynamic> permip(ip);
    PermutationMatrix<Dynamic,Dynamic> permjp(jp);
    PermutationMatrix<Dynamic,Dynamic> permim(im);
    PermutationMatrix<Dynamic,Dynamic> permjm(jm);
    return (u*permip-2*u+u*permim)/(dx*dx)+(permjp*u-2*u+permjm*u)/(dy*dy);
}
#endif
