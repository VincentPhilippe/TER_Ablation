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
  double cfl, dt , e = 10, flux, a;
  int n=0, i, num_cell;
  MatrixXd C1, erreur;
  VectorXd aire = VectorXd::Zero(interf.maxCoeff()+1);;

  C1 = _concentration;

  for(int i=0; i< _maillage.GetNz(); i++){
    for(int j=0; j< _maillage.GetNx(); j++){
      if(_plic->Get_ninterface()(i,j) >= 0){
        aire(_plic->Get_ninterface()(i,j)) = aireInterf(i,j);
      }
    }
  }
  cfl = 0.4;

  dt = cfl*aire.minCoeff();
  _vitesse = VectorXd::Zero(interf.maxCoeff()+1);

  cout<<"CONCENTRATION AVANT"<<endl<<_concentration<<endl;


  while(e>10e-5 && n<10000)
  {
    for(int j = 0; j < _maillage.GetNx(); j++){

      i = 0;
      flux = 0;

      while((_plic->Get_ninterface())(i,j) == -2){
        flux += fluxGauche(i,j);
        flux += fluxBas(i,j);
        flux += fluxDroite(i,j);
        flux += fluxHaut(i,j);
        C1(i,j) = _concentration(i,j) + (dt/(dx*dz))*flux;

        i++;
      }

      // Condition limite interface : calcul des 4 flux + flux interface~ -Da * C
      while(_plic->Get_ninterface()(i,j) >= 0 && i <= _maillage.GetNz() )
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


        _vitesse(num_cell) = C1(i,j);

        i++;
      }
    }
    erreur = (C1-_concentration).cwiseAbs();
    e = erreur.maxCoeff();

    cout<<"~~~~~~ERREUR= "<<e<<"~~~~~~~~~~"<<endl;

    _concentration = C1;
    if(e>100){
      cout<<"ERREUR TROP GRANDE"<<endl;
      exit(0);
    }

  }
  if(n >= 10000)
    cout<<"SORTIE CAR TROP D'ITERATION"<<endl;



  saveCFluid();
}

double diffusion::fluxGauche(int i, int j)
{
  double flux = 0;
  if(j == 0)
  {
      flux = (_concentration(i,j)-_concentration(i,_maillage.GetNx()-1))/dx;
      flux *= longueurArete(i,j,LEFT);
      return(-flux);
  }
  /*if(_plic->Get_ninterface()(i,j-1) > 0 && _plic->Get_ninterface()(i,j) == 0)
  {
    flux = -_damkohler(j-1)*_concentration(i,j);
    flux *= dz - longueurArete(i,j,LEFT);
  }*/

  flux += (_concentration(i,j)-_concentration(i,j-1))/dx;
  flux *= longueurArete(i,j,LEFT);

  return(-flux);
}

double diffusion::fluxBas(int i, int j)
{
  double flux;
  if(i == _maillage.GetNz()-1)
  {
      flux = 0;
      return(flux);
  }
  flux = (_concentration(i+1,j)-_concentration(i,j))/dz;
  flux *= longueurArete(i,j,DOWN);

  return(flux);
}

double diffusion::fluxDroite(int i, int j)
{
  double flux;
  if(j == _maillage.GetNx()-1)
  {
      flux = (_concentration(i,0)-_concentration(i,j))/dx;
      flux *= longueurArete(i,j,RIGHT);
      return(flux);
  }
  /*if(_plic->Get_ninterface()(i,j+1) > 0 && _plic->Get_ninterface()(i,j) == 0)
  {
    flux = -_damkohler(j+1)*_concentration(i,j);
    flux *= dz - longueurArete(i,j,RIGHT);
  }*/

  flux += (_concentration(i,j+1)-_concentration(i,j))/dx;
  flux *= longueurArete(i,j,RIGHT);

  return(flux);
}

double diffusion::fluxHaut(int i, int j)
{
  double flux;
  if(i == 0)
  {
      flux = (_concentration(i,j) - 1)/dz;
      flux *= longueurArete(i,j,UP);
      return(-flux);
  }

  flux = (_concentration(i,j)-_concentration(i-1,j))/dz;
  flux *= longueurArete(i,j,UP);

  return(-flux);
}

double diffusion::fluxInterf(int i, int j)
{
  double x1, z1, x2, z2;
  double Da, l;
  int num_cell = (int)(_plic->Get_ninterface()(i,j));



  x1 = (_plic->Get_interface())(num_cell, 0);
  z1 = (_plic->Get_interface())(num_cell, 1);
  x2 = (_plic->Get_interface())(num_cell, 2);
  z2 = (_plic->Get_interface())(num_cell, 3);
  l = sqrt( (x1 - x2)*(x1 - x2) + (z1 - z2)*(z1 - z2) );

  Da = _damkohler(j);

  return(-Da*_concentration(i,j)*l);
}

double diffusion::aireInterf(int i, int j)
{
  double Mx, Mz, Nx, Nz, aire;
  int num_cell = (int)(_plic->Get_ninterface()(i,j));

  if((_plic->Get_normal())(num_cell, 0) >= 0)
  {
    Mx = (_plic->Get_interface())(num_cell, 0);
    Mz = (_plic->Get_interface())(num_cell, 1);
    Nx = (_plic->Get_interface())(num_cell, 2);
    Nz = (_plic->Get_interface())(num_cell, 3);

    aire = abs(Mz*Nx) + abs((Mx-Nx)*Mz) + abs((Nz-Mz)*Nx) + 0.5*abs((Mx-Nx)*(Nz-Mz));
    return(dx*dz - aire);
  }
  else
  {

    Mx = (_plic->Get_interface())(num_cell, 0);
    Mz = (_plic->Get_interface())(num_cell, 1);
    Nx = (_plic->Get_interface())(num_cell, 2);
    Nz = (_plic->Get_interface())(num_cell, 3);

    aire = 0.5*abs((Mz-Nz)*(Mx-Nx)) + abs(Nz*Mx) + abs(Nz*(dx-Mx)) + abs((Mz-Nz)*(dx-Mx));
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

      x1 = (_plic->Get_interface())(num_cell, 0);
      z1 = (_plic->Get_interface())(num_cell, 1);
      x2 = (_plic->Get_interface())(num_cell, 2);
      z2 = (_plic->Get_interface())(num_cell, 3);


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

enum State_Interf diffusion::watchInterf(int i, int j, enum Direction direction)
{
  VectorXd point1, point2;
  int num_cell = (int)(_plic->Get_ninterface()(i,j));

  point1.resize(2);
  point2.resize(2);
  if(num_cell == -2)
  {
    return(A);
  }
  if(num_cell == -1)
  {
    return(S);
  }
  else
  {
    point1(0) = (_plic->Get_interface())(num_cell, 0);
    point1(1) = (_plic->Get_interface())(num_cell, 1);

    point2(0) = (_plic->Get_interface())(num_cell, 2);
    point2(1) = (_plic->Get_interface())(num_cell, 3);

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

void diffusion::saveCFluid()
{

  string name_file = "CFluid.dat";
  ofstream concentration;
	concentration.open(name_file, ios::out);
	concentration.precision(7);
  concentration << _concentration<<'\n';
  concentration.close();

}
