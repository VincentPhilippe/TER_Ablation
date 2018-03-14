#include <iostream>

#include "Dense"
#include "maillage.h"
#include "ablation.h"

using namespace std;
using namespace Eigen;

diffusion::diffusion(Maillage2DCarre& maillage)
{
  _maillage = maillage;
  _data = data;
  _surface = surface;

  _concentration = data.GetC0();
  _vitesse = MatrixXd::Zero(maillage.GetNx());
}

void diffusion::resolution() //Résolution de dC/dt = d2C/dx2
{
  double dt = 0.4*_maillage.Getdx(), erreur = 10;
  int n;
  MatrixXd C1;
  C1 = MatrixXd::Zero(_maillage.GetNz(), _maillage.GetNx());

  while(erreur > 1e-5, n < 100000)
  {

  }
}
