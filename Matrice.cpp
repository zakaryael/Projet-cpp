#include "Matrice.h"
#include <iostream>
#include <cmath>
using namespace std;

Matrice::Matrice(int n, int m)
{
  dims[0] = n;
  dims[1] = m;
  mat = new Vecteur[m];
  for(int j = 0; j < m; j++)
  {
    mat[j] = Vecteur(n);
  }
}

Matrice::Matrice(int m, Vecteur* v)
{
  //check all vectors are the same dimension? and both tables have the same number of vectors
  int n = (*v).dim;
  mat = new Vecteur[m];
  dims[0] = n; dims[1] = m;
  copy(v,v + m, mat);
}

Matrice::Matrice(Vecteur v)
{
  int d = v.dim;
  dims[0] = d; dims[1] = d;
  mat = new Vecteur[d];
  for(int j = 0; j < d; j++)
  {
    mat[j] = Vecteur(d);
    mat[j][j] = v[j];
  }
}

Matrice::Matrice()
{
  dims[0] = 0;
  dims[1] = 1;
  mat =  nullptr;
}

Matrice::~Matrice()
{
  delete[] mat;
}

Matrice::Matrice(const Matrice &M)
{
  dims[0] = M.dims[0];
  dims[1] = M.dims[1];
  mat = new Vecteur[dims[1]];
  copy(M.mat, M.mat + dims[1], mat);
}

void Matrice::affiche()
{
  cout<<"______________ \n";
  cout << "dim: "<< dims[0]<< "x"<< dims[1] << "." << endl;
  cout << "matrice: \n";
  cout << "[";
  for(int i = 0; i < dims[0]; i++)
  {
    cout << "[";
    for(int j = 0; j < dims[1]; j++)
    {
      cout << mat[j][i]; if(j < dims[1]-1) cout << ",";
    }
    cout << "]";
    if(i < dims[0]-1) cout <<endl;
  }
  cout << "]" << endl;
  cout<<"______________ \n";

}

Matrice & Matrice::operator=(const Matrice &M)
{
  if(this != &M)
  {
    delete[] mat;
    dims[0] = M.dims[0]; dims[1] = M.dims[1];
    mat = new Vecteur[dims[1]];
    copy(M.mat, M.mat + dims[1], mat);
  }
  return *this;
}

Vecteur& Matrice::operator[](int j)
{
  //check for j in 0..m
  return mat[j];
}
Matrice Matrice::operator+(Matrice M)
{
  //check for dimensions compatibility...later
  Matrice sum(dims[0], dims[1]);
  for(int j = 0; j < dims[1]; j++)
  {
    sum[j] = mat[j] + M.mat[j];
  }
  return sum;
}

Matrice Matrice::operator-(Matrice M)
{
  Matrice diff(dims[0], dims[1]);
  for(int j = 0; j < dims[1]; j++)
  {
    diff[j] = mat[j] - M.mat[j];
  }
  return diff;
}

Matrice operator*(float alpha, Matrice M)
{
  Matrice prod(M.dims[0], M.dims[1]);
  for(int j = 0; j < M.dims[1]; j++)
  {
    prod[j] = alpha * M.mat[j];
  }
  return prod;
}

Vecteur Matrice::mvprod(Vecteur v)
{
  //check dimensions as always
  Vecteur prod(dims[0]);
  for(int i = 0; i < dims[0]; i++)
  {
    for(int j = 0; j < dims[1]; j++)
    {
      prod[i] += mat[j][i] * v[j];
    }

  }
  return prod;
}

Matrice Matrice::transpose()
{
  Matrice t(dims[1], dims[0]);
  for(int j = 0; j < dims[1]; j++)
  {
    for(int i = 0; i < dims[0]; i++)
    {
      t[i][j] = mat[j][i];
    }
  }
  return t;
}

Matrice Matrice::submat(int i_l,int j_l, int i_c, int j_c)
{
  Matrice sub(j_l - i_l + 1, j_c - i_c + 1);
  for(int j = 0; j < j_c - i_c + 1; j++)
  {
    sub[j] = ((*this)[i_c + j]).subvec(i_l, j_l);
  }
  return sub;
}

Matrice Matrice::operator*(Matrice M)
{
  int n = dims[0]; int m = M.dims[1]; // nombre de lignes/colonnes de la matrice produit
  Matrice prod(n, m);
  for(int j = 0; j < m; j ++)
  {
    for(int i = 0; i < n; i++)
    {
      for(int h = 0; h < dims[1]; h++)
      {
        prod[j][i] += mat[h][i] * M[j][h];
      }
    }
  }
  return prod;
}

float Matrice::norm()
{
  float norm_squared = 0;
  for(int j = 0; j < dims[1]; j++)
  {
    norm_squared += mat[j].dot(mat[j]);
  }
  return sqrt(norm_squared);
}

Matrice outer(Vecteur u, Vecteur v)
{
  //creation des deux matrices colonnes correspondantes:
  Vecteur Ut[] = {u}; Matrice U(1, Ut);
  Vecteur Vt[] = {v}; Matrice V(1, Vt);
  //retourner leur prduit externe:
  return U * V.transpose();
}

Matrice identity(int d)
{
  //Matrice M = Matrice(d,d);
  return Matrice(ones(d));
}


void Matrice::submod(int i_l, int j_l, int i_c, int j_c, Matrice M)
{
  for(int j = 0; j < j_c - i_c + 1; j++)
  {
    ((*this).mat[i_l+j]).submod(i_l, j_l, M[j]);
  }
}

bool Matrice::isdiag()
{
  for(int j = 0; j < dims[1]; j++)
  {
    for(int i = 0; i < dims[0]; i++)
    {
      if((j != i) && (mat[j][i] != 0)) //(abs(mat[j][i]) > pow(10, -6)))
      {
        return false;
      }
    }
  }
  return true;
}

bool Matrice::isnull()
{
  for(int j = 0; j < dims[1]; j++)
  {
    for(int i = 0; i < dims[0]; i++)
    {
      if(mat[j][i] != 0) //(abs(mat[j][i]) > pow(10, -6)))
      {
        return false;
      }
    }
  }
  return true;
}
