#include "diffusion.h"


using namespace std;
using namespace Eigen;

diffusion::diffusion(read_data& data, Cartesien2D& maillage)
 : _data(data), _maillage(maillage)
{
  dx = _data.Get_dx();
  dz = _data.Get_dz();

  _concentration = _data.Get_C0();
  _vitesse = VectorXd::Zero(_maillage.GetNx());
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
      while(_plic->Get_interface()(i,j) == 0){
        flux += fluxGauche(i,j);
        flux += fluxBas(i,j);
        flux += fluxDroite(i,j);
        flux += fluxHaut(i,j);

        C1(i,j) = _concentration(i,j) + (dt/dx*dz)*flux;
        erreur += abs(C1(i,j) - _concentration(i,j));
        j++;
      }
// Condition limite interface




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
      flux = -(_concentration(_maillage.GetNx(),j)-_concentration(i,j))/dx;
      flux *= longueurArete(i,j,LEFT);
      break;

    default :
      flux = -(_concentration(i-1,j)-_concentration(i,j))/dx;
      flux *= longueurArete(i,j,LEFT);
      break;
  }
  return(flux);
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
=======
      flux = -(_concentration(i,j+1)-_concentration(i,j))/dz;
      flux *= longueurArete(i,j,DOWN);
>>>>>>> 2425892898eab55c23d397bc2f06acb7bf772fbf
      flux = -(_concentration(i,j+1)-_concentration(i,j))/dz;
      flux *= longueurArete(i,j,DOWN);
      break;
  }
  return(flux);
}

double diffusion::fluxDroite(int i, int j)
{
  double flux;
  switch(watchCell(i+1,j))
  {
    case BORD_DROIT:
      flux = -(_concentration(0,j)-_concentration(i,j))/dx;
      flux *= longueurArete(i,j,RIGHT);
      break;

    default :
      flux = -(_concentration(i+1,j)-_concentration(i,j))/dx;
      flux *= longueurArete(i,j,RIGHT);
      break;
  }
  return(flux);
}

double diffusion::fluxHaut(int i, int j)
{
  double flux;
  switch(watchCell(i,j-1))
  {
    case BORD_HAUT:
      flux = -(1-_concentration(i,j))/dz;
      flux *= longueurArete(i,j,UP);
      break;

    default :
      flux = -(_concentration(i,j-1)-_concentration(i,j))/dz;
      flux *= longueurArete(i,j,UP);
      break;
  }
  return(flux);
}

double diffusion::longueurArete(int i, int j, enum Direction direction)
{
  double x1, z1, x2, z2;
  int num_cell = int()_plic->Get_ninterface()(i,j);

  switch(watchInterf(i, j, direction))
  {
<<<<<<< HEAD
    case A:
      if(direction == UP || direction == DOWN)
        return(dx);
      return(dz);
      break;
=======
      num_cell = _plic->Get_interface()(i,j);
      liste(0,k) = _plic->Get_interface()(0,num_cell-1);
      liste(1,k) = _plic->Get_interface()(1,num_cell-1);
      k++;
      liste(0,k) = _plic->Get_interface()(2,num_cell-1);
      liste(1,k) = _plic->Get_interface()(3,num_cell-1);
      k++;
>>>>>>> 2425892898eab55c23d397bc2f06acb7bf772fbf

    case S:
      return(0);
      break;

    default:
      x1 = _plic->Get_interface()(0,numcell-1);
      z1 = _plic->Get_interface()(1,numcell-1);
      x2 = _plic->Get_interface()(2,numcell-1);
      z2 = _plic->Get_interface()(3,numcell-1);
      if(direction == LEFT)
      {
<<<<<<< HEAD
        if(x1 == 0)
        {
          return(dz - z1);
        }
        return(dz - z2);
      }

      if(direction == RIGHT)
      {
        if(x1 == dx)
        {
          return(dz - z1);
        }
        return(dz - z2);
      }

      if(direction == UP)
      {
        if(z1 == dz)
        {
          return(dx - x1);
        }
        return(dx - x2);
      }

      if(direction ==DOWN)
      {
        if(z1 == 0)
        {
          return(dx - x1);
        }
        return(dx - x2);
      }
      break;
=======
        num_cell = _plic->Get_interface()(i,j);
        liste(0,k) = _plic->Get_interface()(0,num_cell-1);
        liste(1,k) = _plic->Get_interface()(1,num_cell-1);
        k++;
        liste(0,k) = _plic->Get_interface()(2,num_cell-1);
        liste(1,k) = _plic->Get_interface()(3,num_cell-1);
        k++;
      }

>>>>>>> 2425892898eab55c23d397bc2f06acb7bf772fbf
  }

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

  if(_plic->Get_interface()(i,j) > 0){
    return(INTERFACE);
  }
  if(_plic->Get_interface()(i,j) == -1){
    return(SOLIDE);
  }

  return(AIR);

}

enum State_Interf diffusion::watchInterf(int i, int j, enum Direction direction)
{
  Vector< double, 2 > point1, point2;
  int num_cell = int()_plic->Get_ninterface()(i,j);

  if(num_cell == 0)
  {
    return(A);
  }
  if(num_cell < 0)
  {
    return(S);
  }
  else
  {
    point1(0) = _plic->Get_interface()(0,numcell-1);
    point1(1) = _plic->Get_interface()(1,numcell-1);

    point2(0) = _plic->Get_interface()(2,numcell-1);
    point2(1) = _plic->Get_interface()(3,numcell-1);

    switch(direction)
    {
      case LEFT:
        if(point1(0) == 0 || point2(0) == 0)
        {
          return(AS);
        }
        if(point1(1) == dz || point2(1) == dz)
        {
          return(S);
        }
        return(A);
        break;

        case RIGHT:
          if(point1(0) == dx || point2(0) == dx)
          {
            return(AS);
          }
          if(point1(1) == dz || point2(1) == dz)
          {
            return(S);
          }
          return(A);
          break;

        case UP:
          if(point1(1) == dz || point2(1) == dz)
          {
            return(AS);
          }
          return(A);
          break;

        case UP:
          if(point1(1) == dz || point2(1) == dz)
          {
            return(AS);
          }
          return(A);
          break;

        case DOWN:
          if(point1(1) == 0 || point2(1) == 0)
          {
            return(AS);
          }
          return(S);
    }
  }
}
