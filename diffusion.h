#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <iostream>

#include "Dense"
#include "maillage.h"
#include "ablation.h"

class diffusion{

private:
  read_data& _data;
  Maillage& _maillage;
  plic& _plic;


  Eigen::MatrixXd _concentration;
  Eigen::VectorXd _vitesse;


public:

  Diffusion(read_data& data, Maillage2DCarre& maillage);
  void Resolution();
  void Vitesse();
  ~diffusion(){};
  bool Watch(int i, int j);
  void update(plic& plic){_plic = plic};
  Eigen::MatrixXd GetConcentration() {return _concentration};
  Eigen::MatrixXd GetVitesse() {return _vitesse};
};


#endif
