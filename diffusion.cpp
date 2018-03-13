#include <iostream>

#include "Dense"
#include "maillage.h"
#include "read_data.h"
#include "rebuild_surface.h"

using namespace std;
using namespace Eigen;

diffusion::diffusion(Maillage2DCarre& maillage, ReadData& data)
{
  _maillage = maillage;
  _data = data;
  _surface = surface;

  _concentration = data.GetC0();
  _vitesse = MatrixXd::Zero(maillage.GetNx());
}

void diffusion::resolution() //Résolution de dC/dt = d2C/dx2
{
  double dt = 0.4*dx, erreur = 10;
  int n;
  MatrixXd C1;
  C1 = MatrixXd::Zero(maillage.GetNz(), maillage.GetNx());

  while(erreur > 1e-5, n < 100000)
  {
    
  }
}
