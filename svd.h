#ifndef SVD_H
#define SVD_H

#include "Matrice.h"
#include "Vecteur.h"

void givens(float x, float z, float& c, float& s);
Vecteur householder(Vecteur x, float& beta);
int signe(float);
Matrice reductridiag(Matrice &);
void qrsym(Matrice &A, Matrice &Q);
Matrice qrpivot(Matrice &A, Matrice &Q);
#endif
