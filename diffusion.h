#include <iostream>

#include "Dense"
#include "maillage.h"
#include "ablation.h"

class diffusion{

private:
  Maillage2DCarre& _maillage;

  Eigen::MatrixXd _concentration;
  Eigen::VectorXd _vitesse;


public:

  diffusion(Maillage2DCarre& maillage);
  void resolution();
  void vitesse();
  ~diffusion(){};

  Eigen::MatrixXd GetConcentration() {return _concentration};
  Eigen::MatrixXd GetVitesse() {return _vitesse};
};
