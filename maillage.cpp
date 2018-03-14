#include "maillage.h"

using namespace Eigen;

// constructeur maillage par défaut
Maillage::Maillage(){}

//destructeur maillage
Maillage::~Maillage(){}

//constructeur maillage cartésien par défaut
Cartesien::Cartesien() : Maillage()
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

  for (int i = 0 ; i < _nbX ; i++)
  {
    _coordX(i) = _deltaX*i + _deltaX/2;
  }

  for (int i = 0 ; i < _nbZ ; i++)
  {
    _coordZ(i) = _deltaZ*i + _deltaZ/2;
  }

};

//constructeur maillage cartésien
Cartesien::Cartesien(std::string name, double deltaX, double deltaZ, double Lx, double Lz) : Maillage()
{
  _name = name;
  _deltaX = deltaX;
  _deltaZ = deltaZ;
  _Lx = Lx;
  _Lz = Lz;

  _nbX = _Lx/_deltaX;
  _nbZ = _Lz/_deltaZ;

  _coordX.resize(_nbX);
  _coordZ.resize(_nbZ);
  _indices.resize(_nbX*_nbZ);

  for (int i = 0 ; i < _nbX*_nbZ ; i++)
  {
    _indices(i) = i;
  }

  for (int i = 0 ; i < _nbX ; i++)
  {
    _coordX(i) = deltaX*i + deltaX/2;
  }

  for (int i = 0 ; i < _nbZ ; i++)
  {
    _coordZ(i) = deltaZ*i + deltaZ/2;
  }

};

//destructeur
Cartesien::~Cartesien(){};
