#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <iostream>
#include <cmath>

#include "Dense"
#include "read_data.h"
#include "maillage.h"
#include "plic.h"

enum State_Case{
  BORD_HAUT,BORD_BAS, AIR, BORD_DROIT, BORD_GAUCHE, INTERFACE, SOLIDE
};

enum Direction{
  LEFT, BOTTOM, RIGHT, UP
};

class diffusion{

private:
  read_data& _data;
  Cartesien2D& _maillage;
  plic *_plic;
  double dx, dz;


  Eigen::MatrixXd _concentration;
  Eigen::VectorXd _vitesse;


public:

  diffusion(read_data& data,Cartesien2D& maillage);
  ~diffusion(){};
  void resolution();
  void vitesse();
  double fluxGauche(int i, int j);
  double fluxBas(int i, int j);
  double fluxDroite(int i, int j);
  double fluxHaut(int i, int j);
  double longueurArete(int k, int l, enum Direction direction);
  enum State_Case watch(int i, int j);
  void update(plic *plic){_plic = plic;}
  Eigen::MatrixXd GetConcentration() {return _concentration;}
  Eigen::MatrixXd GetVitesse() {return _vitesse;}
};


#endif
