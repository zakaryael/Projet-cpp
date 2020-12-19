#include "svd.h"
#include <cmath>
#include <iostream>
using namespace std;

void givens(float x, float z, float& c, float& s)
{
  if(z == 0)
  {
    c = 1;
    s = 0;
  }
  else
  {
    if(abs(z) > abs(x))
    {
      float tau = - x / z;
      s = 1 / sqrt(1 + pow(tau, 2));
      c = s * tau;
    }
    else
    {
      float tau = -z / x;
      c = 1 / sqrt(1 + pow(tau, 2));
      s = c * tau;
    }
  }
}

Vecteur householder(Vecteur x, float& beta)
{
  int n = x.len();
  Vecteur y = x.subvec(1, n-1);
  float sigma = pow(y.norm(), 2);
  Vecteur v(n);
  v[0] = 1;
  v.submod(1, n-1, y);
  if(sigma == 0 & x[0] >= 0)
  {
    beta = 0;
  }
  else
  {
    float mu = sqrt(pow(x[0], 2) + sigma);
    if(x[0] <= 0)
    {
      v[0] = x[0] - mu;
    }
    else
    {
      v[0] = -sigma / (x[0] + mu);
    }
    beta = 2 * pow(v[0], 2) / (sigma + pow(v[0], 2));
    v = (1 / v[0]) * v;
  }
  return v;
}

int signe(float d)
{
  if(d >= 0) return 1;
  else return -1;
}

Matrice reductridiag(Matrice &D)
{
  int l = D.shape()[0];
  float d = (D[l-2][l-2] - D[l-1][l-1]) / 2;
  float mu = D[l-1][l-1] - pow(D[l-2][l-1], 2) / (d + signe(d) * sqrt(pow(d, 2) + pow(D[l-2][l-1], 2)));
  float x = D[0][0] - mu;
  float z = D[0][1];
  Matrice Z (identity(l));
  for(int k = 1; k <= l-1; k++)
  {
    float c,s;
    givens(x, z, c, s);
    //application de la rotation de givens a droite
    for(int j = 1; j <= l; j++)
    {
      float tau1 = D[k-1][j-1];
      float tau2 = D[k][j-1];

      D[k-1][j-1] = c * tau1 - s * tau2;
      D[k][j-1] = s * tau1 + c * tau2;
      tau1 = Z[k-1][j-1];
      tau2 = Z[k][j-1];
      Z[k-1][j-1] = c * tau1 - s * tau2;
      Z[k][j-1] = s * tau1 + c * tau2;
    }
    //applicaiton de la rotation de givens a gauche
    for(int j = 1; j <= l; j++)
    {
      float tau1 = D[j-1][k-1];
      float tau2 = D[j-1][k];
      D[j-1][k-1] = c * tau1 - s * tau2;
      D[j-1][k] = s * tau1 + c * tau2;
    }
    if(k < l-1)
    {
      x = D[k-1][k];
      z = D[k-1][k+1];
    }
  }
  return Z;
}
void qrsym(Matrice &A, Matrice &Q)
{

  int n = A.shape()[0];
  Q = identity(n);

  //Tridiagonalisation de la matrice: A = QT
  for(int k = 1; k <= n-2; k++)
  {
    float beta;
    Vecteur x ( A.submat(k, n-1, k-1, k-1)[0] );
    Vecteur v ( householder(x, beta) );
    Vecteur p ( beta * A.submat(k, n-1, k, n-1).mvprod(v) );
    Vecteur w ( p - (beta / 2) * (p.dot(v)) * v );
    A[k-1][k] = x.norm();
    A[k][k-1] = A[k-1][k];
    Matrice M ( A.submat(k, n-1, k, n-1) - outer(v, w) - outer(w, v) );
    A.submod(k, n-1, k, n-1, M);
    Matrice Qsub ( Q.submat(k, n-1, k, n-1) );
    M = (identity(n-k) - beta * outer(v, v)) * Qsub;
    Q.submod(k, n-1, k, n-1, M);
  }
  Matrice T(n, n);
  for(int j = n; j >= 1; j--)
  {
    T[j-1][j-1] = A[j-1][j-1];
    if(j > 1)
    {
      T[j-1][j-2] = A[j-2][j-1];
      T[j-2][j-1] = T[j-1][j-2];
    }
  }
  //diagoalisation de T et mise a jour de Q:
  int counter1 = 0;
while(!T.isdiag() && counter1 < 10) //c'est plus efficace d'implementer une autre fct de test de diagonalite pour les matrice tridiag
{
  for(int i = 1; i < n; i++)
  {
    if(abs(T[i][i-1]) + abs(T[i-1][i]) <= (abs(T[i-1][i-1]) + abs(T[i][i])) * pow(10, -9))
    {
      T[i][i-1] = 0;
      T[i-1][i] = 0;
    }
  }
  int p = 0;
  for(; T.submat(0, p, 0, p).isdiag() && T.submat(p+1, n-1, p+1, n-1).isnull()  ;p++);
//les test sont a raffiner car on sait deja que T est tridiagonale, pour gagner en efficacité.
  int q = 0;
  for(; T.submat(n-1-q, n-1, n-1-q, n-1).isdiag()  && T.submat(0, n-q-2, 0, n-q-2).isnull() && q < 10 ;q++);
  if(p + q < n)
  {
    Matrice T2 (T.submat(p, n-q-1, p, n-q-1));
    Matrice Z ( reductridiag(T2));
    Matrice T_hat (T);
    T_hat.submod(p, n-q-1, p, n-q-1, T2);
    Matrice R (identity(n));
    R.submod(p, n-q-1, p, n-q-1, Z);

    Q = Q * R;
    T = 0.5 * (T_hat + T_hat.transpose());
  }
}
}

Matrice qrpivot(Matrice &A, Matrice &Q)
{
  int n = A.shape()[1]; int m = A.shape()[0];
  Matrice Pi(identity(m));
  Vecteur c(n);
  for(int j = 0; j < n; j++)
  {
    c[j] = A[j].dot(A[j]);
  }
  int r = 0;
  float tau = max(c);
  for(int r = 0; r < n && tau > 0; r++)
  {
    int k = r + argmax(c.subvec(r, n-1));
    if(k != r)
    {
      Vecteur temp1(A[r]);
      A[r] = A[k]; A[k] = temp1;
      float temp2 = c[r];
      c[r] = c[k]; c[k] = temp2;
      Vecteur temp3(Pi[r]);
      Pi[r] = Pi[k]; Pi[k] = temp3;
    }
    float beta = 2;
    Vecteur x(A.submat(r, m-1, r, r)[0]);
    Vecteur v (householder(x, beta));
    Matrice A_sub(A.submat(r, m-1, r, m-1));
    Matrice M(A_sub - beta * outer(v, v) * A_sub);
    A.submod(r, m-1, r, m-1, M);

    (A[r]).submod(r, m-1, v.subvec(1, m-r-1));//pb ici dans l'algo
    for(int i = r; i < n; i++)
    {
     c[i] = c[i] - pow(A[r][i], 2);
    }

    if(r < n-1) tau = max(c.subvec(r+1, n-1));//pb ici dans l'olgo
    else tau = 0;
  }
  //Calcul de Q
  Q = identity(m);
  Vecteur v(m);
  for(int j = m; j > 0; j--)
  {
    v[j-1] = 1;
    v.submod(j, m-1, A[j-1].subvec(j, m-1));
    float beta = 2 / v.dot(v);
    Vecteur x(v.subvec(j-1, m-1));
    float b = 2 / x.dot(x);
    Matrice Q_sub(Q.submat(j-1, m-1, j-1, m-1));
    Matrice M(Q_sub - b * outer(x, x) * Q_sub);
  }
  return Pi;
}

void svd(Matrice &A, Matrice *M)
{

}
