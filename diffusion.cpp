#include "diffusion.h"


using namespace std;
using namespace Eigen;

diffusion::Diffusion(read_data& data, Maillage& maillage)
{
  _maillage = &maillage;
  _data = &data;

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

enum State_Case diffusion::Watch(int i, int j) // Regarde l'état de la case (i,j)
{
  if(i == 0){
    return(BORD_HAUT);
  }

  if(i == _maillage.GetNx()+1){
    return(BORD_BAS);
  }

  if(j == 0){
    return(BORD_GAUCHE);
  }

  if (j == _maillage.GetNz()+1){
    return(BORD_DROIT);
  }

  if(_plic->Get_ninterf[i,j] != 0){
    return(INTERFACE);
  }
  if()

  return(AIR);

}
