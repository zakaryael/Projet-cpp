#ifndef TENSEUR_H
#define TENSEUR_H

#include "Matrice.h"

class Tenseur
{
public:
  Tenseur();
  Tenseur(int d, int *dims);
  Tenseur(int d, int *dims, Vecteur v);
  ~Tenseur();
  void affiche();
  Tenseur(const Tenseur &);
  Tenseur & operator=(const Tenseur &);
  Tenseur(int d, int *dims, int k, Matrice A);
  float& operator[](int i);
  Tenseur operator+(Tenseur T);
  Tenseur operator-(Tenseur T);
  Matrice mode(int k);
  Tenseur pmod(Matrice M, int k);

private:
  int order;
  int * dims;
  int nbelts;
  Vecteur elts;
};

int comp_nbelts(int d, int * dimensions);
int phi(int d, int * dims, int * indices);
int * prod_dims(int d, int * dimensions);
int * phi_inv(int d, int * dims, int i);
void pop_k(int , int *tab, int *tab2, int k);
void push_k(int *tab, int k, int ik);
void print_tab(int size, int *);
#endif
