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
        _Nx = floor(_Lx/_dx);
        _dx = _Lx/_Nx;
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
        data_file >> _flux; _if_flux = true;
      }

      if (file_line.find("Dim") != std::string::npos)
      {
        data_file >> _dim; _if_dim = true;
      }

      if (file_line.find("Da") != std::string::npos)
      {
        std::string Da_string;
        data_file >> Da_string; _if_Da = true;
        if (Da_string == "uniforme")
        {
          double Da_unif, Nx_temp;
          data_file >> Da_unif;
          _Da.resize(_Nx+1,2);
          _Da(0,0)=_Nx;
          for (int i=1; i<_Nx+1; i++)
          {
            _Da(i,0) = i*_dx;
            _Da(i,1) = Da_unif;
          }
        }
        else if (Da_string == "step")
        {
          double Da_step1, Da_step2, Nx_temp;
          data_file >> Da_step1 >> Da_step2;
          _Da.resize(_Nx+1,2);
          _Da(0,0)=_Nx;
          for (int i=1; i<floor(_Nx/3); i++)
          {
            _Da(i,0) = i*_dx;
            _Da(i,1) = Da_step1;
          }
          for (int i=floor(_Nx/3); i<2*floor(_Nx/3); i++)
          {
            _Da(i,0) = i*_dx;
            _Da(i,1) = Da_step2;
          }
          for (int i=2*floor(_Nx/3); i<_Nx+1; i++)
          {
            _Da(i,0) = i*_dx;
            _Da(i,1) = Da_step1;
          }
        }
        else if (Da_string == "retrieve")
        {
          std::string _Da_file_name;
          data_file >> _Da_file_name;
          ifstream Da_file(_Da_file_name.data());
          if (!Da_file.is_open())
          {
            cout << "Unable to open Da file " << _Da_file_name << endl;
            abort();
          }
          else
          {
            cout << "-------------------------------------------------" << endl;
            cout << "Reading Da data file " << _Da_file_name << endl;
          }
          int Da_size(0);
          double _x, _Da_x;
          Da_file >> Da_size;
          _Da.resize(Da_size,2);
          _Da(0,0)=Da_size;
          for (int i=1; i<_Nx+1; i++)
          {
            Da_file >> _x >> _Da_x;
            _Da(i,0)=_x;
            _Da(i,1)=_Da_x;
          }
          Da_file.close();
        }
        else
        {
          cout << "Try again." << endl;
          abort();
        }
      }

      if (file_line.find("C0") != std::string::npos)
      {
        std::string C0_string;
        data_file >> C0_string; _if_C0 = true;
        if (C0_string == "uniforme")
        {

        }
        else if (C0_string == "retrieve")
        {

        }
        else
        {
          cout << "Try again." << endl;
          abort();
        }
      }

      if (file_line.find("Surface") != std::string::npos)
      {
        std::string Surface_string;
        data_file >> Surface_string; _if_Surface = true;
        if (Surface_string == "uniforme")
        {

        }
        else if (Surface_string == "retrieve")
        {

        }
        else
        {
          cout << "Try again." << endl;
          abort();
        }
      }

  }

}

#define _READ_DATA_CPP
#endif
