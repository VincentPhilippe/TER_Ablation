#ifndef _READ_DATA_CPP

#include "read_data.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

read_data::read_data(std::string file_name)
: _file_name(file_name),  _if_dx(false), _if_dz(false), _if_Lx(false),
_if_Lz(false), _if_dt(false), _if_tfinal(false), _if_D(false), _if_flux(false),
_if_dim(false), _if_Da(false), _if_C0(false), _if_Surface(false)
{}

  void read_data::read_datafile()
  {
    ifstream data_file(_file_name.data());
    if (!data_file.is_open())
    {
      cout << "Unable to open file " << _file_name << endl;
      abort();
    }
    else
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Reading data file " << _file_name << endl;
    }

    string file_line;

    while (!data_file.eof())
    {
      getline(data_file, file_line);
      if (file_line.find("dx") != std::string::npos)
      {
        data_file >> _dx; _if_dx = true;
      }

      if (file_line.find("dz") != std::string::npos)
      {
        data_file >> _dz; _if_dz = true;
      }

      if (file_line.find("Lx") != std::string::npos)
      {
        data_file >> _Lx; _if_Lx = true;
      }

      if (file_line.find("Lz") != std::string::npos)
      {
        data_file >> _Lz; _if_Lz = true;
      }

      if (file_line.find("dt") != std::string::npos)
      {
        data_file >> _dt; _if_dt = true;
      }

      if (file_line.find("tfinal") != std::string::npos)
      {
        data_file >> _tfinal; _if_tfinal = true;
      }

      if (file_line.find("D") != std::string::npos)
      {
        data_file >> _D; _if_D = true;
      }

      if (file_line.find("Flux") != std::string::npos)
      {
        _if_flux = true;
      }

      if (file_line.find("Dim") != std::string::npos)
      {
        _if_dim = true;
      }

      if (file_line.find("Da") != std::string::npos)
      {
        _if_Da = true;
      }

      if (file_line.find("C0") != std::string::npos)
      {
        _if_C0 = true;
      }

      if (file_line.find("Surface") != std::string::npos)
      {
        _if_Surface = true;
      }

  }

}

#define _READ_DATA_CPP
#endif
