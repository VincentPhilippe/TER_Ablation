#ifndef _READ_DATA_CPP

#include "read_data.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

read_data::read_data(std::string file_name)
: _file_name(file_name),  _if_dx(false), _if_dz(false), _if_Lx(false),
_if_Lz(false), _if_dt(false), _if_tfinal(false), _if_Diff(false), _if_flux(false),
_if_dim(false), _if_Da(false), _if_C0(false), _if_C_Solide(false)
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
      if (file_line.find("dx") != std::string::npos){data_file >> _dx; _if_dx = true;}

      if (file_line.find("dz") != std::string::npos){data_file >> _dz; _if_dz = true;}

      if (file_line.find("Lx") != std::string::npos)
      {
        data_file >> _Lx; _if_Lx = true;
        _Nx = floor(_Lx/_dx)+1;
        _dx = _Lx/(_Nx-1);
      }

      if (file_line.find("Lz") != std::string::npos)
      {
        data_file >> _Lz; _if_Lz = true;
        _Nz = floor(_Lz/_dz)+1;
        _dz = _Lz/(_Nz-1);
      }

      if (file_line.find("dt") != std::string::npos){data_file >> _dt; _if_dt = true;}

      if (file_line.find("tfinal") != std::string::npos){data_file >> _tfinal; _if_tfinal = true;}

      if (file_line.find("Diff") != std::string::npos){data_file >> _Diff; _if_Diff = true;}

      if (file_line.find("Flux") != std::string::npos){data_file >> _flux; _if_flux = true;}

      if (file_line.find("Dim") != std::string::npos){data_file >> _dim; _if_dim = true;}

      if (file_line.find("Da") != std::string::npos)
      {
        std::string Da_string;
        data_file >> Da_string; _if_Da = true;
        if (_dim == "2D")
        {
          if (Da_string == "uniforme")
          {
            double Da_unif;
            data_file >> Da_unif;
            _Da.resize(_Nx+1,2);
            _Da(0,0)=_Nx;
            for (int i=1; i<_Nx+1; i++)
            {
              _Da(i,0) = (i-1)*_dx;
              _Da(i,1) = Da_unif;
            }
          }
          else if (Da_string == "step")
          {
            double Da_step1, Da_step2;
            data_file >> Da_step1 >> Da_step2;
            _Da.resize(_Nx+1,2);
            _Da(0,0)=_Nx;
            for (int i=1; i<=floor(_Nx/3); i++)
            {
              _Da(i,0) = (i-1)*_dx;
              _Da(i,1) = Da_step1;
            }
            for (int i=floor(_Nx/3)+1; i<=2*floor(_Nx/3); i++)
            {
              _Da(i,0) = (i-1)*_dx;
              _Da(i,1) = Da_step2;
            }
            for (int i=2*floor(_Nx/3)+1; i<_Nx+1; i++)
            {
              _Da(i,0) = (i-1)*_dx;
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
            _Da.resize(Da_size+1,2);
            _Da(0,0)=Da_size;
            for (int i=1; i<_Nx+1; i++)
            {
              Da_file >> _x >> _Da_x;
              _Da(i,0)=_x;
              _Da(i,1)=_Da_x;
            }
            cout << "End of reading Da data file" << endl;
            cout << "-------------------------------------------------" << endl;
            Da_file.close();
          }
          else
          {
            cout << "Try again." << endl;
            abort();
          }
        }
        else if (_dim == "3D")
        {
          if (Da_string == "uniforme")
          {
            double Da_unif;
            data_file >> Da_unif;
            _Da.resize(_Nx+1,_Nz+1);
            _Da(0,1)=_Nx;
            _Da(1,0)=_Nz;
            for (int i=1; i<_Nx+1; i++)
            {
              for (int j=1; j<_Nz+1; j++)
              {
              _Da(i,0) = (i-1)*_dx;
              _Da(0,j) = (j-1)*_dz;
              _Da(i,j) = Da_unif;
              }
            }
          }
          else if (Da_string == "step")
          {
            double Da_step1, Da_step2;
            data_file >> Da_step1 >> Da_step2;
            _Da.resize(_Nx+1,_Nz+1);
            _Da(0,1)=_Nx;
            _Da(1,0)=_Nz;
            for (int i=1; i<=floor(_Nx/3); i++)
            {
              for (int j=1; j<=floor(_Nz/3); j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step1;
              }
              for (int j=floor(_Nz/3)+1; j<=2*floor(_Nz/3); j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step2;
              }
              for (int j=2*floor(_Nz/3)+1; j<_Nz+1; j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step1;
              }
            }
            for (int i=floor(_Nx/3)+1; i<=2*floor(_Nx/3); i++)
            {
              for (int j=1; j<=floor(_Nz/3); j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step2;
              }
              for (int j=floor(_Nz/3)+1; j<=2*floor(_Nz/3); j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step1;
              }
              for (int j=2*floor(_Nz/3)+1; j<_Nz+1; j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step2;
              }
            }
            for (int i=2*floor(_Nx/3)+1; i<_Nx+1; i++)
            {
              for (int j=1; j<=floor(_Nz/3); j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step1;
              }
              for (int j=floor(_Nz/3)+1; j<=2*floor(_Nz/3); j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step2;
              }
              for (int j=2*floor(_Nz/3)+1; j<_Nz+1; j++)
              {
                _Da(i,0) = (i-1)*_dx;
                _Da(0,j) = (j-1)*_dz;
                _Da(i,j) = Da_step1;
              }
            }
          }
          else if (Da_string == "retrieve")
          {
            ////////////////////////////////////////////////////////////////////////////////////
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
            _Da.resize(Da_size+1,2);
            _Da(0,0)=Da_size;
            for (int i=1; i<_Nx+1; i++)
            {
              Da_file >> _x >> _Da_x;
              _Da(i,0)=_x;
              _Da(i,1)=_Da_x;
            }
            cout << "End of reading Da data file" << endl;
            cout << "-------------------------------------------------" << endl;
            Da_file.close();
          }
          else
          {
            cout << "Try again." << endl;
            abort();
          }
        }
      }

      if (file_line.find("C0") != std::string::npos)
      {
        std::string C0_string;
        data_file >> C0_string; _if_C0 = true;
        if (_dim == "2D")
        {
            if (C0_string == "uniforme")
            {
              double C0_unif;
              data_file >> C0_unif;
              _C0.resize(_Nx+1,2);
              _C0(0,0)=_Nx;
              for (int i=1; i<_Nx+1; i++)
              {
                _C0(i,0) = (i-1)*_dx;
                _C0(i,1) = C0_unif;
              }
            }
            else if (C0_string == "step")
            {
              double C0_step1, C0_step2;
              data_file >> C0_step1 >> C0_step2;
              _C0.resize(_Nx+1,2);
              _C0(0,0)=_Nx;
              for (int i=1; i<=floor(_Nx/3); i++)
              {
                _C0(i,0) = (i-1)*_dx;
                _C0(i,1) = C0_step1;
              }
              for (int i=floor(_Nx/3)+1; i<=2*floor(_Nx/3); i++)
              {
                _C0(i,0) = (i-1)*_dx;
                _C0(i,1) = C0_step2;
              }
              for (int i=2*floor(_Nx/3)+1; i<_Nx+1; i++)
              {
                _C0(i,0) = (i-1)*_dx;
                _C0(i,1) = C0_step1;
              }
            }
            else if (C0_string == "retrieve")
            {
              std::string _C0_file_name;
              data_file >> _C0_file_name;
              ifstream C0_file(_C0_file_name.data());
              if (!C0_file.is_open())
              {
                cout << "Unable to open C0 file " << _C0_file_name << endl;
                abort();
              }
              else
              {
                cout << "-------------------------------------------------" << endl;
                cout << "Reading C0 data file " << _C0_file_name << endl;
              }
              int C0_size(0);
              double _x, _C0_x;
              C0_file >> C0_size;
              _C0.resize(C0_size+1,2);
              _C0(0,0)=C0_size;
              for (int i=1; i<_Nx+1; i++)
              {
                C0_file >> _x >> _C0_x;
                _C0(i,0)=_x;
                _C0(i,1)=_C0_x;
              }
              cout << "End of reading C0 data file" << endl;
              cout << "-------------------------------------------------" << endl;
              C0_file.close();
            }
            else
            {
              cout << "Try again." << endl;
              abort();
            }
          }
        else if (_dim == "3D")
        {}
      }

      if (file_line.find("C_Solide") != std::string::npos)
      {
        std::string C_Solide_string;
        data_file >> C_Solide_string; _if_C_Solide = true;
        if (_dim == "2D")
        {
              if (C_Solide_string == "uniforme")
              {
                double C_Solide_unif;
                data_file >> C_Solide_unif;
                _C_Solide.resize(_Nz,_Nx);
                for (int j=0; j<_Nx; j++)
                {
                  int i=0.;
                  while (i*_dz<C_Solide_unif)
                  {
                    _C_Solide(i,j)=1.;
                    i++;
                  }
                  _C_Solide(i,j)=0.5;
                  i++;
                  while (i*_dz<_Lz)
                  {
                    _C_Solide(i,j)=0.;
                    i++;
                  }
                }
              }
              else if (C_Solide_string == "step")
              {
                double C_Solide_step1, C_Solide_step2;
                data_file >> C_Solide_step1 >> C_Solide_step2;
                _C_Solide.resize(_Nz,_Nx);
                for (int j=1; j<=floor(_Lx/3); j++)
                {
                  int i=0.;
                  while (i*_dz<C_Solide_step1)
                  {
                    _C_Solide(i,j)=1.;
                  }
                  _C_Solide(i,j)=0.5;
                  while (i*_dz<_Nz)
                  {
                    _C_Solide(i,j)=0.;
                  }
                  i++;
                }
                for (int j=floor(_Lx/3)+1; j<=2*floor(_Lx/3); j++)
                {
                  int i=0.;
                  while (i*_dz<C_Solide_step2)
                  {
                    _C_Solide(i,j)=1.;
                  }
                  _C_Solide(i,j)=0.5;
                  while (i*_dz<_Nz)
                  {
                    _C_Solide(i,j)=0.;
                  }
                  i++;
                }
                for (int j=2*floor(_Lx/3)+1; j<_Lx+1; j++)
                {
                  int i=0.;
                  while (i*_dz<C_Solide_step1)
                  {
                    _C_Solide(i,j)=1.;
                  }
                  _C_Solide(i,j)=0.5;
                  while (i*_dz<_Lz)
                  {
                    _C_Solide(i,j)=0.;
                  }
                  i++;
                }
              }
              else if (C_Solide_string == "retrieve")
              {
                std::string _C_Solide_file_name;
                data_file >> _C_Solide_file_name;
                ifstream C_Solide_file(_C_Solide_file_name.data());
                if (!C_Solide_file.is_open())
                {
                  cout << "Unable to open C_Solide file " << _C_Solide_file_name << endl;
                  abort();
                }
                else
                {
                  cout << "-------------------------------------------------" << endl;
                  cout << "Reading C_Solide data file " << _C_Solide_file_name << endl;
                }
                _C_Solide.resize(_Nz,_Nx);
                for (int i=0; i<_Nz; i++)
                {
                  for (int j=0; j<_Nx; j++)
                  {
                    C_Solide_file >> _C_Solide(i,j);
                  }
                }
                cout << "End of reading C_Solide data file" << endl;
                cout << "-------------------------------------------------" << endl;
                C_Solide_file.close();
              }
              else
              {
                cout << "Try again." << endl;
                abort();
              }
            }
        else if (_dim == "3D")
        {}
      }

    }

  }

#define _READ_DATA_CPP
#endif
