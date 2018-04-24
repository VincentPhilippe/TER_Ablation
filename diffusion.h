#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "objet.h"
#include <iostream>
#include <cmath>

#include "Dense"
#include "read_data.h"
#include "maillage.h"
#include "plic.h"
#include "recul.h"

enum State_Cell{
  BORD_HAUT, BORD_BAS, AIR, BORD_DROIT, BORD_GAUCHE, INTERFACE, SOLIDE
};

enum State_Interf{
  A, S, AS
};

enum Direction{
  LEFT, DOWN, RIGHT, UP
};

class diffusion{

private:
  read_data& _data;
  Cartesien2D& _maillage;
  plic *_plic;
  double dx, dz;


  Eigen::MatrixXd _concentration;
  Eigen::VectorXd _vitesse;
  Eigen::VectorXd _damkohler;


public:

  diffusion(read_data& data,Cartesien2D& maillage);
  ~diffusion(){};

  void resolution();

  double fluxGauche(int i, int j);
  double fluxBas(int i, int j);
  double fluxDroite(int i, int j);
  double fluxHaut(int i, int j);
  double fluxInterf(int i, int j);
  double aireInterf(int i, int j);
  double longueurArete(int k, int l, enum Direction direction);
  enum State_Cell watchCell(int i, int j);
  enum State_Interf watchInterf(int i, int j, enum Direction direction);
  void update(plic *plic){ _plic = plic; }

  Eigen::MatrixXd GetConcentration() { return _concentration; }
  Eigen::MatrixXd GetVitesse() { return _vitesse; }
};


#endif
