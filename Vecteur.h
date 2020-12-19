#ifndef VECTEUR_H
#define VECTEUR_H

class Vecteur
{
public:
  Vecteur(int d, float* t);
  Vecteur(int d);
  Vecteur();
  void affiche();
  Vecteur(const Vecteur &); // constructeur de copie
  float& operator[](int i);
  Vecteur operator+(Vecteur  v);
  Vecteur operator-(Vecteur  v);
  Vecteur operator*(float alpha); // add a friend func or replace with
  friend Vecteur operator*(float alpha, Vecteur v);
  float dot(Vecteur  v);
  float norm();
  Vecteur subvec(int i, int j);
  Vecteur & operator=(const Vecteur &);
  ~Vecteur();
  float sum();
  int len(){return dim;} //fonction qui renvoie la dimension du vecteur
  void submod(int i, int j, Vecteur x);


private:
  //members:
  int dim;
  float* tab;
  friend class Matrice;
  friend class Tenseur;
};
Vecteur ones(int d);
float max(Vecteur );
int argmax(Vecteur );
#endif
