#include "maillage.h"

using namespace Eigen;

//constructeur par défaut
void Maillage2DCarre()
{
  _name = default;
  _deltaX = 0.1;
  _deltaZ = 0.1;
  _Lx = 10;
  _Lz = 10;

  _coordonnees.resize(2, (_Lx/_deltaX)*(_Lz/_deltaZ))

  for (int i = 0 ; i < _Lx ; i++)
  {
      for (int j = 0 ; j < _Ly : i++)
      {
        _coordonnees(1, j) = deltaX*i + deltaX/2;
        _coordonnees(2, j) = deltaZ*j + deltaZ/2;
      }
  }
}

//constructeur
void Maillage2DCarre(std::string name, double deltaX, double deltaZ, double _Lx, double _Lz)
{
  _name = name;
  _deltaX = deltaX;
  _deltaZ = deltaZ;
  _Lx = Lx;
  _Lz = Lz;

  _coordonnees.resize(2, (_Lx/_deltaX)*(_Lz/_deltaZ))

  for (int i = 0 ; i < _Lx ; i++)
  {
      for (int j = 0 ; j < _Ly : i++)
      {
        _coordonnees(1, j) = deltaX*i + deltaX/2;
        _coordonnees(2, j) = deltaZ*j + deltaZ/2;
      }
  }

}
