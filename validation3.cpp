#include <iostream>
#include <cmath>
using namespace std;
#include "svd.h"

int main(int argc, char const *argv[]) {

  //1.
  float c,s;
  givens(1, 2, c, s);
  cout << c << endl << s << endl;
  //2.
  float beta = 2;
  Vecteur x(2); x[0] = -1;
  Vecteur y(2); y[0] = 1 / sqrt(2); y[1] = y[0];
  Vecteur z(1); z[0] = -4;

  householder(x, beta).affiche();
  cout << beta << endl;

  householder(y, beta).affiche();
  cout << beta << endl;

  householder(z, beta).affiche();
  cout << beta << endl;
  ones(5).affiche();
  cout << identity(5).shape()[1] << endl;

  Matrice D = 10 * identity(2);
  D[0][1] = -6;
  D[1][0] = -6;
  D.affiche();
  reductridiag(D);

  Matrice I = identity(2);
  qrsym(D, I);
  I.affiche();
  D.affiche();
  float t1[3] = {1, 0.5, 1.0/3};
  Vecteur c1(3, t1);
  float t2[3] = {0.5, 1.0/3, 0.25};
  Vecteur c2(3, t2);
  float t3[] = {1.0/3, 0.25, 0.2};
  Vecteur c3(3, t3);
  Vecteur m[] = {c1, c2, c3};
  Matrice M(3, m);

  Matrice Q = identity(3);
  Matrice Pi = qrpivot(M, Q);

  return 0;
}
