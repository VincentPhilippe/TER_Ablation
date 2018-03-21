#include "maillage.h"

using namespace Eigen;

// constructeur maillage par défaut
Maillage::Maillage(){}

//destructeur maillage
Maillage::~Maillage(){}

//constructeur maillage cartésien 2D par défaut
Cartesien2D::Cartesien2D() : Maillage()
{
  _name = "default";
  _deltaX = 0.1;
  _deltaZ = 0.1;
  _Lx = 10;
  _Lz = 10;

  _nbX = _Lx/_deltaX;
  _nbZ = _Lz/_deltaZ;

  _coordX.resize(_nbX);
  _coordZ.resize(_nbZ);
  _indices.resize(_nbX*_nbZ);

  for (int i = 0 ; i < _nbX*_nbZ ; i++)
  {
    _indices(i) = i;
  }

  for (int i = 0 ; i < _nbZ ; i++)
  {
    for (int j = 0 ; j < _nbX ; j++)
    {
      _coordX(i*_nbX+j) = _deltaX*j + _deltaX/2;
      _coordZ(i*_nbX+j) = _deltaZ*i + _deltaZ/2;
    }
  }
};

//constructeur maillage 2D cartésien
Cartesien2D::Cartesien2D(std::string name, double deltaX, double deltaZ, double Lx, double Lz) : Maillage()
{
  _name = name;
  _deltaX = deltaX;
  _deltaZ = deltaZ;
  _Lx = Lx;
  _Lz = Lz;

  _nbX = _Lx/_deltaX;
  _nbZ = _Lz/_deltaZ;

  _coordX.resize(_nbX*_nbZ);
  _coordZ.resize(_nbX*_nbZ);
  _indices.resize(_nbX*_nbZ);

  for (int i = 0 ; i < _nbX*_nbZ ; i++)
  {
    _indices(i) = i;
  }

  for (int i = 0 ; i < _nbZ ; i++)
  {
    for (int j = 0 ; j < _nbX ; j++)
    {
      _coordX(i*_nbX+j) = _deltaX*j + _deltaX/2;
      _coordZ(i*_nbX+j) = _deltaZ*i + _deltaZ/2;
    }
  }
};

//constructeur maillage cartésien 3D par défaut
Cartesien3D::Cartesien3D() : Maillage()
{
  _name = "default";
  _deltaX = 0.1;
  _deltaY = 0.1;
  _deltaZ = 0.1;
  _Lx = 10;
  _Ly = 10;
  _Lz = 10;

  _nbX = _Lx/_deltaX;
  _nbY = _Ly/_deltaY;
  _nbZ = _Lz/_deltaZ;

  _coordX.resize(_nbX);
  _coordY.resize(_nbY);
  _coordZ.resize(_nbZ);
  _indices.resize(_nbX*_nbZ*_nbY);

  for (int i = 0 ; i < _nbX*_nbZ*_nbY ; i++)
  {
    _indices(i) = i;
  }

  for (int i = 0 ; i < _nbZ ; i++)
  {
    for (int j = 0 ; j < _nbY ; j++)
    {
      for (int k = 0 ; k < _nbX ; k++)
      {
      _coordX(i*_nbX*_nbY+k) = _deltaX*k + _deltaX/2;
      _coordY(i*_nbX*_nbY+k) = _deltaY*j + _deltaY/2;
      _coordZ(i*_nbX*_nbY+k) = _deltaZ*i + _deltaZ/2;
      }
    }
  }
};

//constructeur maillage 3D cartésien
Cartesien3D::Cartesien3D(std::string name, double deltaX, double deltaY, double deltaZ, double Lx, double Ly, double Lz) : Maillage()
{
  _name = name;
  _deltaX = deltaX;
  _deltaY = deltaY;
  _deltaZ = deltaZ;
  _Lx = Lx;
  _Ly = Ly;
  _Lz = Lz;

  _nbX = _Lx/_deltaX;
  _nbY = _Ly/_deltaY;
  _nbZ = _Lz/_deltaZ;

  _coordX.resize(_nbX*_nbZ*_nbY);
  _coordZ.resize(_nbX*_nbZ*_nbY);
  _coordZ.resize(_nbX*_nbZ*_nbY);
  _indices.resize(_nbX*_nbZ*_nbY);

  for (int i = 0 ; i < _nbZ ; i++)
  {
    for (int j = 0 ; j < _nbY ; j++)
    {
      for (int k = 0 ; k < _nbX ; k++)
      {
      _coordX(i*_nbX*_nbY+k*_nbY) = deltaX*k + deltaX/2;
      _coordY(i*_nbX*_nbY+k*_nbY) = deltaY*j + deltaY/2;
      _coordZ(i*_nbX*_nbY+k*_nbY) = deltaZ*i + deltaZ/2;
      }
    }
  }
};

//destructeurs
Cartesien2D::~Cartesien2D(){};
Cartesien3D::~Cartesien3D(){};
