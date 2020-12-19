#include "Vecteur.h"

#include <iostream>
#include <cmath>
using namespace std;

int main(int argc, char const *argv[]) {
  /* code */
  //1 creation et affichage de u et v:
  float u_t[3] = {1, 1, 1};
  Vecteur u(3, u_t);
  cout << "Affichage du vecteur u: \n";
  u.affiche();

  float v_t[4] = {3, 4, 0, 0};
  Vecteur v(4, v_t);
  cout << "\nAffichage du vecteur v: \n";
  v.affiche();

  //2.copie du vecteur u dans t:
  Vecteur t(u);
  cout << "\n \nAffichage du vecteur t: \n";
  t.affiche();

  //3.modification de u et affichage de u et t:
  u[2] = 0;
  cout << "\n \nAffichage du vecteur u modifié: \n";
  u.affiche();
  cout << "\nAffichage du vecteur t après modification de u: \n \n";
  t.affiche();

  //4.Calcul du produit scalair v'v:
  float produit = v.dot(v);
  cout << "le produit vTv: " << produit << endl;
  float norm_v = v.norm();
  //5.calcul et affichage de v / norm(v):
  float alpha = 1 / norm_v;
  Vecteur unit = v * alpha;
  cout << "\n \n Affichage de v / norm(v): \n";
  unit.affiche();
  assert(unit.norm() == 1); // s'assurer que u est un vecteur unitair

  //6.clacul de w et affichage de v et w:

  Vecteur w = v.subvec(1,3);
  cout << "\n \n Affichage de v: \n";
  v.affiche();
  cout << "\n \n Affichage de w: \n";
  w.affiche();

  //7. calcul et affichage de u+w et u-w

  Vecteur uplusw = u + w;
  Vecteur umoinsw = u - w;
  cout << "\n \n Affichage de u+w: \n";
  uplusw.affiche();
  cout << "\n \n Affichage de u-w: \n";
  umoinsw.affiche();
  return 0;
}
