#include "diffusion.h"


using namespace std;
using namespace Eigen;

diffusion::diffusion(read_data& data, Maillage& maillage)
 : _data(data), _maillage(maillage)
{
  dx = _data.Get_dx();
  dz = _data.Get_dz();

  _concentration = _data.Get_C0();
  _vitesse = MatrixXd::Zero(maillage.GetNx());
}

void diffusion::resolution() //Résolution de dC/dt = d2C/dx2
{
  double dt = 0.4*(dx+dz), erreur = 10, flux;
  int n=0, i, j;
  MatrixXd C1;
  C1 = MatrixXd::Zero(_maillage.GetNz(), _maillage.GetNx());

  while(erreur>10e-9)
  {
    for(int i = 1; i < _maillage.GetNx(); i++){
      j = 1;
      flux = 0;
      while(plic->Get_ninterf()[i,j] == 0){
        flux += fluxGauche(i,j);
        flux += fluxBas(i,j);
        flux += fluxDroite(i,j);
        flux += fluxHaut(i,j);

        C1[i,j] = _concentration[i,j] + (dt/dx*dz)*flux;

        erreur += abs(C1[i,j] - _concentration[i,j]);

        j++;

      }
    }
  }

    _concentration = C1;

}

double diffusion::fluxGauche(int i, int j)
{
  double flux;
  switch(watch(i-1,j))
  {
    case BORD_GAUCHE:
      flux = -(_concentration[Nx,j]-_concentration[i,j])/dx;
      flux *= longueurArete(i,j,LEFT);
      break;

    default :
      flux = -(_concentration[i-1,j]-_concentration[i,j])/dx;
      flux *= longueurArete(i,j,LEFT);
      break;
  }
}

double diffusion::fluxBas(int i, int j)
{
  double flux;
  switch(watch(i,j+1))
  {
    case BORD_BAS:
      flux = 0;
      break;

    default :
      flux = -(_concentration[i,j+1]-_concentration[i,j])/dz;
      flux *= longueurArete(i,j,BOTTOM);
      break;
  }
}

double diffusion::fluxDroite(int i, int j)
{
  double flux;
  switch(watch(i+1,j))
  {
    case BORD_DROIT:
      flux = -(_concentration[0,j]-_concentration[i,j])/dx;
      flux *= longueurArete(i,j,RIGHT);
      break;

    default :
      flux = -(_concentration[i+1,j]-_concentration[i,j])/dx;
      flux *= longueurArete(i,j,RIGHT);
      break;
  }
}

double diffusion::fluxHaut(int i, int j)
{
  double flux;
  switch(watch(i,j-1))
  {
    case BORD_HAUT:
      flux = -(1-_concentration[i,j])/dz;
      flux *= longueurArete(i,j,UP);
      break;

    default :
      flux = -(_concentration[i,j-1]-_concentration[i,j])/dz;
      flux *= longueurArete(i,j,UP);
      break;
  }
}

double diffusion::longueurArete(int k, int l, enum Direction direction)
{



  }
}

void diffusion::vitesse() // Permet de donner la vitesse normale en chaque point de la surface
{

}

enum State_Case diffusion::watch(int i, int j) // Regarde l'état de la case (i,j)
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

  if(_plic->Get_ninterf()[i,j] > 0){
    return(INTERFACE);
  }
  if(_plic->Get_ninterf()[i,j] = -1){
    return(SOLIDE);
  }

  return(AIR);

}
