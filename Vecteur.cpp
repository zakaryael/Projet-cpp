#include "Vecteur.h"
#include <iostream>
#include <cmath>
using namespace std;

void Vecteur::affiche()
{
  cout<<"______________ \ndim: " << dim << endl;
  cout<<"vecteur: [";
  for(int i = 0; i < dim; i++)
  {
    cout << *(tab + i); if(i < dim-1) cout << "; ";
  }
  cout << "]" << endl << "______________ \n\n";

}

Vecteur::Vecteur(int d, float* table)
: dim(d)
{
  if(d < 0) {throw runtime_error("la dimension doit etre un entier positif ");}
  tab = new float[d];
  copy(table,table + d, tab);
}

Vecteur::Vecteur()
:dim(0), tab(nullptr)
{
  //tab =  nullptr;
}
Vecteur::Vecteur(int d): dim(d)
{
if(d < 0) {throw runtime_error("la dimension doit etre un entier positif ");}
  tab = new float[d];
  for(int i = 0; i < d; i++)
  {
    tab[i] = 0;
  }
}

float& Vecteur::operator[](int i)
{
  if(i < 0 || i >= dim) throw runtime_error("out of range");
  return tab[i];
}

Vecteur Vecteur::operator+(Vecteur v)
{

  if(v.dim != dim) {throw runtime_error("dimensions incompatibles");}
  Vecteur sum(dim);
  for(int i = 0; i < dim; i++)
  {
    sum.tab[i] = v.tab[i] + tab[i];
  }
  return sum;
}

Vecteur Vecteur::operator-(Vecteur v)
{
  if(v.dim != dim) {throw runtime_error("dimensions incompatibles");}
  Vecteur diff(dim);
  for(int i = 0; i < dim; i++)
  {
    diff.tab[i] = tab[i] - v.tab[i];
  }
  return diff;
}

Vecteur Vecteur::operator*(float alpha)
{
  //check the two vectors are of the same dimension
  Vecteur scal(dim);
  for(int i = 0; i < dim; i++)
  {
    scal.tab[i] = alpha * tab[i];
  }
  return scal;
}
Vecteur operator*(float alpha, Vecteur v)
{
  Vecteur scal(v.dim);
  for(int i = 0; i < v.dim; i++)
  {
    scal.tab[i] = alpha * v.tab[i];
  }
  return scal;
}
float Vecteur::dot(Vecteur v)
{
  if(v.dim != dim) {throw runtime_error("dimensions incompatibles");}
  float prod = 0;
  for(int i = 0; i < dim; i++)
  {
    prod += v.tab[i] * tab[i];
  }
  return prod;
}

float Vecteur::norm()
{
  return sqrt((*this).dot(*this));
}


Vecteur Vecteur::subvec(int i, int j)
{
  Vecteur subv(j-i+1);
  for(int l = 0; l < j-i+1; l++)
  {
    subv[l] = (*this)[i + l];
  }
  return subv;
}

Vecteur::Vecteur(const Vecteur &v)
{
  dim = v.dim;
  tab = new float[dim];
  copy(v.tab,v.tab + dim, tab);
}

Vecteur & Vecteur::operator=(const Vecteur &v)
{
  if(this != &v)
  {
    delete[] tab;
    dim = v.dim;
    tab = new float[dim];
    copy(v.tab,v.tab + dim, tab);
  }
  return *this;
}

Vecteur::~Vecteur()
{
  delete[] tab;
}

void Vecteur::submod(int i, int j, Vecteur x)
{
  // takes two indices i and j, a vector x of length j-i+1 and
  //assings x to the subvector i:j of our original vector
  for(int l = 0; l < j-i+1; l++)
  {
    tab[l + i] = x[l];
  }
}

float Vecteur::sum()
{
  //computes the sum of the vector's components
  float s = 0;
  for(int i = 0; i < dim; i++)
  {
    s += tab[i];
  }
  return s;
}

Vecteur ones(int d)
{//returns a d-dimensional vector of ones
  if(d<0) throw runtime_error("la dimension doit etre un entier positif ");
  Vecteur v(d);
  for(int i = 0; i < d; i++)
  {
    v[i] = 1;
  }
  return v;
}

float max(Vecteur v)
{ // retourne le max des elts de v
  float m = v[0];
  int n = v.len();
  for(int i = 0; i < n; i++)
  {
    if(v[i] > m) {m = v[i];}
  }
  return m;
}

int argmax(Vecteur v)
{
  //renvoi le plus petit argmax d v
  float m = v[0];
  int index = 0;
  for(int i = 0; i < v.len(); i++)
  {
    int n = v.len();
    if(v[i] > m) {index = i;}
  }
  return index;
}
