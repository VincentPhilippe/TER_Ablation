#ifndef _READ_DATA_H

#include <string>
#include <vector>
#include <iostream>
// Définition de la classe

class read_data {
private:
  std::string _file_name;
  double _dx, _dz, _Lx, _Lz, _dt, _tfinal, _D;

  std::string _flux;
  std::string _dim;
  std::string _Da;
  std::string _C0;
  std::string _Surface;

  bool _if_dx;
  bool _if_dz;
  bool _if_Lx;
  bool _if_Lz;
  bool _if_dt;
  bool _if_tfinal;
  bool _if_D;
  bool _if_flux;
  bool _if_dim;
  bool _if_Da;
  bool _if_C0;
  bool _if_Surface;

public: // Méthodes et opérateurs de la classe
  read_data(std::string file_name);
  void read_datafile();
  std::string Get_file_name() const {return _file_name;}
  double Get_dx() const {return _dx;};
  double Get_dz() const {return _dz;};
  double Get_dt() const {return _dt;};
  double Get_tfinal() const {return _tfinal;};
  double Get_Lx() const { return _Lx;}
  double Get_Lz() const { return _Lz;}
  double Get_D() const { return _D;}
  std::string Get_flux() const {return _flux;};
  std::string Get_dim() const {return _dim;};
  std::string Get_Da() const {return _Da;};
  std::string Get_C0() const {return _C0;};
  std::string Get_Surface() const {return _Surface;};
};

#define _READ_DATA_H
#endif
