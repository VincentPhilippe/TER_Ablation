#include <iostream>

#include "Dense"
#include "maillage.h"
#include "read_data.h"
#include "rebuild_surface.h"

class diffusion{

private:
  Maillage2DCarre& _maillage;
  ReadData& _data;
  Plic _surface;

  Eigen::MatrixXd _concentration;
  Eigen::VectorXd _vitesse;


public:

  diffusion(Maillage2DCarre& maillage, ReadData& data, Plic& _surface);
  void resolution();
  void vitesse();
  ~diffusion(){};

  Eigen::MatrixXd GetConcentration() {return _concentration};
  Eigen::MatrixXd GetVitesse() {return _vitesse};
};
