#include "diffusion.h"


using namespace std;
using namespace Eigen;

 : _data(data), _maillage(maillage)
{
  dx = _data.Get_dx();
  dz = _data.Get_dz();

  _concentration = _data.Get_C0();
}

void diffusion::resolution() //Résolution de dC/dt = d2C/dx2
{
  double dt = 0.4*(dx+dz), erreur = 10, flux;
  MatrixXd C1;
  C1 = MatrixXd::Zero(_maillage.GetNz(), _maillage.GetNx());

  while(erreur>10e-9)
  {
    for(int i = 1; i < _maillage.GetNx(); i++){
      j = 1;
      flux = 0;
        flux += fluxGauche(i,j);
        flux += fluxBas(i,j);
        flux += fluxDroite(i,j);
        flux += fluxHaut(i,j);



        j++;

      }
    }
  }

    _concentration = C1;

}

double diffusion::fluxGauche(int i, int j)
{
  double flux;
  switch(watchCell(i-1,j))
  {
    case BORD_GAUCHE:
      flux *= longueurArete(i,j,LEFT);
      break;

    default :
      flux *= longueurArete(i,j,LEFT);
      break;
  }
}

double diffusion::fluxBas(int i, int j)
{
  double flux;
  switch(watchCell(i,j+1))
  {
    case BORD_BAS:
      flux = 0;
      break;

    default :
<<<<<<< HEAD
      flux = -(_concentration[i,j+1]-_concentration[i,j])/dz;
      flux *= longueurArete(i,j,DOWN);
=======
      flux = -(_concentration(i,j+1)-_concentration(i,j))/dz;
      flux *= longueurArete(i,j,BOTTOM);
>>>>>>> 4c3dae7aa2af92daa39691ab5488aa149f643f0a
      break;
  }
}

double diffusion::fluxDroite(int i, int j)
{
  double flux;
  switch(watchCell(i+1,j))
  {
    case BORD_DROIT:
      flux *= longueurArete(i,j,RIGHT);
      break;

    default :
      flux *= longueurArete(i,j,RIGHT);
      break;
  }
}

double diffusion::fluxHaut(int i, int j)
{
  double flux;
  switch(watchCell(i,j-1))
  {
    case BORD_HAUT:
      flux *= longueurArete(i,j,UP);
      break;

    default :
      flux *= longueurArete(i,j,UP);
      break;
  }
}

double diffusion::longueurArete(int k, int l, enum Direction direction)
{
<<<<<<< HEAD
  
=======



>>>>>>> 4c3dae7aa2af92daa39691ab5488aa149f643f0a
}

void diffusion::vitesse() // Permet de donner la vitesse normale en chaque point de la surface
{

}

enum State_Cell diffusion::watchCell(int i, int j) // Regarde l'état de la case (i,j)
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

    return(INTERFACE);
  }
    return(SOLIDE);
  }

  return(AIR);

}
