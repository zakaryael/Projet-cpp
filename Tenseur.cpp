#include "Tenseur.h"
#include <iostream>
using namespace std;
Tenseur::Tenseur()
{
  order = 0;
  dims = nullptr;
  nbelts = 0;
  Vecteur v;
  elts = v;
}

int comp_nbelts(int d, int * dimensions) // fonction qui calcule le nb d'elts d'un tenseur..
{
  return prod_dims(d, dimensions)[0];
}

Tenseur::Tenseur(int d, int * dimensions)
:order(d), nbelts(comp_nbelts(d, dimensions))
{
  dims = new int[d];
  std::copy(dimensions, dimensions + d, dims);
  elts = Vecteur(nbelts);
}

void Tenseur::affiche()
{
  std::cout << "========Tenseur========\nordre: " << order << std::endl;
  cout << "les dimensions: ";
  print_tab(order, dims);
  std::cout << "le vecteur contenant les " << nbelts <<  " elts: \n";
  elts.affiche();
}

Tenseur::Tenseur(int d, int *dimensions, Vecteur v)
:order(d), nbelts(comp_nbelts(d, dimensions))
{
  dims = new int[d];
  std::copy(dimensions,dimensions + d, dims);
  elts = v;//elts(v);
}

Tenseur::Tenseur(const Tenseur &T)
{
  order = T.order;
  nbelts = T.nbelts;
  dims = new int[order];
  std::copy(T.dims, T.dims + order, dims);
  elts = T.elts;
}

Tenseur & Tenseur::operator=(const Tenseur &T)
{
  if(this != &T)
  {
    delete[] dims;
    order = T.order;
    nbelts = T.nbelts;
    dims = new int[order];
    std::copy(T.dims, T.dims + order, dims);
    elts = T.elts;
  }
  return *this;
}

Tenseur::~Tenseur()
{
  delete[] dims;
}

int phi(int d, int * dims, int * indices)
{
  //la fonction phi présenté dans le projet
  //les indices(en entré et en sortie) commencent de zéro
  int i = 0;
  int *f = prod_dims(d, dims);
  for(int j = 0; j < d; j++)
  {
    i += indices[j] * f[j+1];
  }
  return i;
}

int * prod_dims(int d, int * dims)
{
  //fonction auxiliaire pour le calcul du nbelts, phi et phi_inv
  //renvoie un tableau contenant n_d*...*n_k, k=1,...,d
  int *f = new int[d+1];
  f[d] = 1;
  for(int k = d; k > 0; k--)
  {
    f[k-1] = dims[k-1] * f[k]; // f[k-1] should = nd*...*nk
    //std::cout<< k << " : " <<f[k-1] << endl; //<< "  fk-1: " <<f[k-1] << "  dimk: " <<dims[k] << std::endl;
  }
  return f;
}

int * phi_inv(int d, int * dims, int i)
{
  //la reciproque de phi
  int *f = prod_dims(d, dims);
  int * indices = new int[d];
  for(int k = 0; k < d; k++)
  {
    indices[k] = i / f[k+1] % dims[k];//formule explicite, pas besoin de proceder par recurrence.
  }
  return indices;
}

Tenseur::Tenseur(int d, int *dimensions, int k, Matrice A)
:order(d), nbelts(comp_nbelts(d, dimensions))
{
  dims = new int[d];
  copy(dimensions,dimensions + d, dims);

  Vecteur v(nbelts);
  for(int i = 0; i < nbelts; i++)
  {
    int * ind = phi_inv(d, dims, i);
    int ik = ind[k-1];
    copy(ind + k, ind + d, ind + k - 1); //enlever la k-eme valeure pour appliquer phi
    copy(dims + k, dims + d, dims + k - 1);
    int jk = phi(d-1, dims, ind);
    v[i] = A[jk][ik];
    //cout << jk;
  }
  //v.affiche();
  //A.affiche();
  elts = v;
}

float & Tenseur::operator[](int i)
{
  return elts[i];
}

Tenseur Tenseur::operator+(Tenseur T)
{
  Tenseur sum(order, dims);
  sum.elts = elts + T.elts;
  return sum;
}

Tenseur Tenseur::operator-(Tenseur T)
{
  Tenseur diff(order, dims);
  diff.elts = elts - T.elts;
  return diff;
}

Matrice Tenseur::mode(int k)
{
  int n = dims[k-1]; int m = nbelts / n;
  int *dims_k = new int[order-1];// dimensions sans la k-eme copmposante
  copy(dims, dims + k - 1, dims_k);
  copy(dims + k, dims + order, dims_k + k - 1);
  int * indices = new int[order];
  Matrice A(n, m);
  for(int j = 0; j < m; j++)
  {
    for(int i = 0; i < n; i++)
    {
      int *ind_k = phi_inv(order-1, dims_k, j);//i1,,,ik-1, ik+1,,,id
      copy(ind_k, ind_k + k - 1 , indices);
      copy(ind_k + k - 1, ind_k + order-1, indices + k);
      indices[k-1] = i;
      A[j][i] = elts[phi(order, dims, indices)];
    }
  }
  return A;
}

void pop_k(int size, int *tab1, int *tab2, int k)
{
  //fonction auxiliaire à utiliser dans mode et le 3eme constructeur pour ameliorer la lisibilité du code
  //met le tableau tab1 ou la k-eme valeure supprimée dans tab2
  copy(tab1, tab1 + k - 1, tab2);
  copy(tab1 + k - 1, tab1 + size - 1, tab2 + k - 1);
}

void print_tab(int size, int *tab)
{
  //affiche les elts d'un tableau d'entiers
  cout << "[";
  for(int i = 0; i < size; i++) cout << tab[i] << "  ";
  cout << "]" << endl;
}

Tenseur Tenseur::pmod(Matrice M, int k)
{
  int *dims_prod = new int[order];
  copy(dims, dims + order, dims_prod);
  int mk = M.shape()[0];
  dims_prod[k-1] = mk;
  Tenseur prod(order, dims_prod);
  int ik;
  for(int i = 0; i < prod.nbelts; i++)
  {
    int *indices = phi_inv(order, dims_prod, i); //correspondants au indices du tenseur produit
    ik = indices[k-1]; //garder le k-eme indice
    for(int j = 0; j < dims[k-1]; j++)
    {
      indices[k-1] = j; //indices maintenant correspondent à ceux du tenseur initial
      prod.elts[i] += M[j][ik] * elts[phi(order, dims, indices)];
    }
  }
  return prod;
}
