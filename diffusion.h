#include <iostream>

#include "dense"
#include "maillage.h"
#include "read_data.h"
#include "rebuild_surface.h"

class diffusion{

private:
  Maillage2DCarre& _maillage;
  ReadData& _data;
  Eigen::MatrixXd _ninterf;
  Eigen::MatrixXd _interface;


  Eigen::MatrixXd _concentration;
  Eigen::VectorXd _vitesse;


public:

  diffusion(Maillage2DCarre& maillage, ReadData& data); // Construteur à l'état initial
  diffusion(Maillage2DCarre& maillage, ReadData& data, Plic& _surface); // Constructeur à l'étape 2
  void resolution();
  void vitesse();

  Eigen::MatrixXd GetConcentration() {return _concentration};
  Eigen::MatrixXd GetVitesse() {return _vitesse};
};
