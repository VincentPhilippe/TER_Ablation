#include "diffusion.h"


using namespace std;
using namespace Eigen;

diffusion::diffusion(read_data& data, Cartesien2D& maillage)
 : _data(data), _maillage(maillage)
{
  dx = _data.Get_dx();
  dz = _data.Get_dz();

  _concentration = _data.Get_C0();
  _damkohler = _data.Get_Da();
}

void diffusion::resolution() //Résolution de dC/dt = d2C/dx2
{
  MatrixXd interf = _plic->Get_ninterface();
  double dt = 0.4*(dx+dz), erreur = 10, flux, a;
  int n=0, j, num_cell;
  MatrixXd C1;
  C1 = MatrixXd::Zero(_maillage.GetNz(), _maillage.GetNx());
  _vitesse = VectorXd::Zero(interf.maxCoeff());

  cout << "ok" << endl;

  while(erreur>10e-9 && n<10000)
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
      // Condition limite interface : calcul des 4 flux + flux interface~ -Da * C
      while(_plic->Get_interface()(i,j+1) != -1 && j <= _maillage.GetNz())
      {
        num_cell = (int)(_plic->Get_ninterface()(i,j));
        flux = 0;
        flux += fluxGauche(i,j);
        flux += fluxBas(i,j);
        flux += fluxDroite(i,j);
        flux += fluxHaut(i,j);
        flux += fluxInterf(i,j);
        a = aireInterf(i,j);
        C1(i,j) = _concentration(i,j) + (dt/a)*flux;

        _vitesse(num_cell-1) = C1(i,j);

        j++;
      }
    }
    n++;
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

double diffusion::fluxInterf(int i, int j)
{
  double x1, z1, x2, z2;
  double Da, l;
  int num_cell = (int)(_plic->Get_ninterface()(i,j));

  x1 = (_plic->Get_interface())(num_cell-1, 0);
  z1 = (_plic->Get_interface())(num_cell-1, 1);
  x2 = (_plic->Get_interface())(num_cell-1, 2);
  z2 = (_plic->Get_interface())(num_cell-1, 3);

  l = sqrt( (x1 - x2)*(x1 - x2) + (z1 - z2)*(z1 - z2) );

  Da = _damkohler(i);

  return(-Da*_concentration(i,j)*l);
}

double diffusion::aireInterf(int i, int j)
{
  double c, d, e, f, aire;
  int num_cell = (int)(_plic->Get_ninterface()(i,j));

  if((_plic->Get_normal())(num_cell-1, 0) >= 0)
  {
    c = (_plic->Get_interface())(num_cell-1, 0);
    d = (_plic->Get_interface())(num_cell-1, 1);
    e = (_plic->Get_interface())(num_cell-1, 2);
    f = (_plic->Get_interface())(num_cell-1, 3);

    aire = abs(d*e) + abs((c-e)*d) + abs((f-d)*e) + 0.5*abs((c-e)*(f-d));
    return(dx*dz - aire);
  }
  else
  {
    c = (_plic->Get_interface())(num_cell-1, 0);
    d = (_plic->Get_interface())(num_cell-1, 1);
    e = (_plic->Get_interface())(num_cell-1, 2);
    f = (_plic->Get_interface())(num_cell-1, 3);

    aire = 0.5*abs((f-d)*(e-c)) + abs(e-c)*d + abs(dx-e)*d + abs((f-d)*(dx-e));
    return(dx*dz - aire);
  }
}

double diffusion::longueurArete(int i, int j, enum Direction direction)
{
  double x1, z1, x2, z2;
  int num_cell = (int)(_plic->Get_ninterface()(i,j));

  switch(watchInterf(i, j, direction))
  {

    case A:
      if(direction == UP || direction == DOWN)
        return(dx);
      return(dz);
      break;

    case S:
      return(0);
      break;

    default:

      x1 = (_plic->Get_interface())(num_cell-1, 0);
      z1 = (_plic->Get_interface())(num_cell-1, 1);
      x2 = (_plic->Get_interface())(num_cell-1, 2);
      z2 = (_plic->Get_interface())(num_cell-1, 3);


      if(direction == LEFT)
      {
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

      }

      return(0);
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
  VectorXd point1, point2;
  int num_cell = (int)(_plic->Get_ninterface()(i,j));

  point1.resize(2);
  point2.resize(2);
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
    point1(0) = (_plic->Get_interface())(num_cell-1, 0);
    point1(1) = (_plic->Get_interface())(num_cell-1, 1);

    point2(0) = (_plic->Get_interface())(num_cell-1, 2);
    point2(1) = (_plic->Get_interface())(num_cell-1, 3);

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

        case DOWN:
          if(point1(1) == 0 || point2(1) == 0)
          {
            return(AS);
          }
          return(S);
    }
  }

  return(S);

}
