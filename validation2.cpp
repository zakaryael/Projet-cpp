#include "Vecteur.h"
#include "Matrice.h"
#include <cassert>
#include <iostream>
using namespace std;


int main(int argc, char const *argv[]) {
  //1.creation et affichage de A et B:
  //____A:
  float c1_t[] = {1, 1, 0};
  float c2_t[] = {-0.5, 2, -1};
  float c3_t[] = {0, -1, 1};
  Vecteur c1(3, c1_t);
  Vecteur c2(3, c2_t);
  Vecteur c3(3, c3_t);
  Vecteur m[] = {c1, c2, c3};
  Matrice A(3, m);
  cout << "affichage de la matrice A: \n";
  A.affiche();
  //____B:
  float bc1_t[] = {-2, 0};
  float bc2_t[] = {3, 1};
  Vecteur bc1(2, bc1_t);
  Vecteur bc2(2, bc2_t);
  Vecteur bm[] = {bc1, bc2};
  Matrice B(2, bm);
  cout << "affichage de la matrice B: \n";
  B.affiche();

  //2.copie de B dans C, modification de B, et affichage:

  Matrice C(B);
  B[1][0] = 0;
  cout << "affichage des matrice B et C: \n";

  B.affiche();
  C.affiche();

  //3.Construction et affichage de D:

  Matrice D = A.submat(0,2,0,1);
  cout << "affichage de la matrice D: \n";
  D.affiche();

  //4.Construction et affichage de E:
  float v_t[3] = {3, 2, 1};
  Vecteur v(3, v_t);
  Matrice E(v);
  cout << "affichage de la matrice E: \n";
  E.affiche();

  //5.Calcul et affichage de B+C, C-D, D*C
  cout << "affichage de la matrice B+C: \n";
  (B+C).affiche();
  cout << "affichage de la matrice C-B: \n";
  (C-B).affiche();
  (D*C).affiche();

  //6.Calcul de la norme de C:
  float norm_de_C = C.norm();

  //7....
  cout << "affichage de la matrice (B+B')/2: \n";
  (0.5 * (B + B.transpose())).affiche();

  //outer(v, v).affiche(); // test de la fct outer

  return 0;
}
