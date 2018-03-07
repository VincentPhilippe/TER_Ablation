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

  int _nbX = _Lx/_deltaX;
  int _nbZ = _Lz/_deltaZ

  _coordX.resize(_nbX);
  _coordZ.resize(_nbZ);

  for (int i = 0 ; i < _nbX ; i++)
  {
    _coordX(i) = deltaX*i + deltaX/2;
  }

  for (int i = 0 ; i < _nbZ ; i++)
  {
    _coordZ(i) = deltaZ*i + deltaZ/2;
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

  int _nbX = _Lx/_deltaX;
  int _nbZ = _Lz/_deltaZ

  _coordX.resize(_nbX);
  _coordZ.resize(_nbZ);

  for (int i = 0 ; i < _nbX ; i++)
  {
    _coordX(i) = deltaX*i + deltaX/2;
  }

  for (int i = 0 ; i < _nbZ ; i++)
  {
    _coordZ(i) = deltaZ*i + deltaZ/2;
  }

}
