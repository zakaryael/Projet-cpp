#ifndef MATRICE_H
#define MATRICE_H
#include "Vecteur.h"

class Matrice
{
public:
  Matrice(int, int);
  Matrice(Vecteur v);
  Matrice(int m, Vecteur* v);
  Matrice();
  ~Matrice();
  Matrice(const Matrice &);
  Matrice & operator=(const Matrice &M);
  Vecteur & operator[](int j);
  Matrice operator+(Matrice);
  Matrice operator-(Matrice);
  Matrice operator*(Matrice);
  friend Matrice operator*(float alpha, Matrice M);
  void affiche();
  Vecteur l(int); //abandon
  Vecteur mvprod(Vecteur v);
  Matrice transpose();
  Matrice submat(int i_l,int j_l, int i_c, int j_c);
  float norm();
  int * shape(){return dims;}
  void submod(int i_l, int j_l, int i_c, int j_c, Matrice M);
  bool isdiag();//test de diagonalité
  bool isnull();//test de nullité

private:
  int dims[2];
  Vecteur* mat;
  friend class Tenseur;
};

Matrice outer(Vecteur u, Vecteur v);
Matrice identity(int);

#endif
