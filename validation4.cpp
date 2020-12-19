#include "Tenseur.h"
#include <iostream>
using namespace std;



int main(int argc, char const *argv[]) {
  //1.
  int d = 3;
  int dims[] = {2, 2, 2};
  Tenseur T(d, dims);
  T.affiche();
  //2.
  Tenseur U(d, dims, ones(8));
  U.affiche();
  //3.
  Tenseur V = U + T;
  V.affiche();

  Tenseur W = U - T;
  W.affiche();

  //4.
  int indices[] = {1, 1, 1};
  U[phi(d, dims, indices)] = -1;
  U.affiche();
  V.affiche();
  //5.

  T.mode(1).affiche();

  //6.
  float c1_t[] = {1, 4};
  float c2_t[] = {3, 1./3};
  float c3_t[] = {0, 1.5};
  float c4_t[] = {-1, 2};
  Vecteur c1(2, c1_t);
  Vecteur c2(2, c2_t);
  Vecteur c3(2, c3_t);
  Vecteur c4(2, c4_t);
  Vecteur m[] = {c1, c2, c3, c4};
  Matrice M(4, m);
  //M.affiche();
  Tenseur G(d, dims, 1, M);
  T = G;
  //G.affiche();
  cout << "6. tenseur T defini par mod2: \n";
  T.affiche();
  T.mode(1).affiche();
  //7.
  float a1_t[] = {3, 0, 0};
  float a2_t[] = {-1, 6, -3};
  Vecteur a1(3, a1_t);
  Vecteur a2(3, a2_t);
  Vecteur a[] = {a1, a2};
  Matrice A(2, a);

  Tenseur S = T.pmod(A, 2);
  S.affiche(); // Je pense il y a une erreur dans le slide de la validation
  S.mode(2).affiche();

  Tenseur R = S + S;
  R.affiche();



  return 0;
}
