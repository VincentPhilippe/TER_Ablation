#include <iostream>

#include "Dense"
#include "maillage.h"


using namespace std;
using namespace Eigen;

diffusion::Diffusion(read_data& data, Maillage& maillage)
{
  _maillage = maillage;
  _data = data;

  _concentration = data.GetC0();
  _vitesse = MatrixXd::Zero(maillage.GetNx());
}

void diffusion::Resolution() //Résolution de dC/dt = d2C/dx2
{
  double dt = 0.4*_maillage.Getdx(), erreur = 10;
  int n=0;
  MatrixXd C1;
  C1 = MatrixXd::Zero(_maillage.GetNz(), _maillage.GetNx());

}

void diffusion::Vitesse() // Permet de donner la vitesse normale en chaque point de la surface
{

}

string diffusion::Watch(int i, int j)
{
  if(i == 0)
  {
    return('')
  }
}
