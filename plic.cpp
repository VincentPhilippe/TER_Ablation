#include "plic.h"
#include <iostream>

using namespace std;

//constructeur
plic::plic()
  :
  {
    _interface.resize(1,4);
  }

  plic::~plic()
  {}

double plic::grad_x(const int i,const int j)
{
  if (i==0)
  {
    return _phi(i+1,j)/2;
  }
  if (i==lon)
  {
    return _phi(i-1,j)/2;
  }
  else
  {
    return (_phi(i+1,j)+_phi(i-1,j))/2;
  }
}

double plic::grad_y(const int i,const int j)
{
  {
    if (j==0)
    {
      return _phi(i,j+1)/2;
    }
    if (j==lar)
    {
      return _phi(i,j-1)/2;
    }
    else
    {
      return (_phi(i,j+1)+_phi(i,j-1))/2;
    }
  }
}

void plic::interf()
{
    _phi=_recul->Get_C_solide();
    int lon=_phi.cols();
    int lar=_phi.rows();
    int k=0;
    _ninterf.resize(lon,lar);
    for (int i=1;i<lon;i++)
     {
        for (int j=1;j<lar;j++)
        {
          p=_phi(i,j);
          if ((p>0.) && (p<1.))   //si on est sur l'interface
          {
              k++;
              _ninterf(i,j)=k;
          }
          else
          {
            _ninterf(i,j)=0;
          }
        }
      }
    //ifinterf.resize(_squares.size());  //ifinterf(i)=1 si _squares(i) interface, 0 sinon
    _interface.resize(k,4) ;
    //k=0;
    _pointsupl=0;
    for (int i=1;i<lon;i++)
     {
        for (int j=1;j<lar;j++)
        {
            p=_phi(i,j);
            if ((p>0.) && (p<1.))   //si on est sur l'interface
            {
                //k++;
                //_ninterf(i)(j)=k;

                //Calcul du gradient
                nx=grad_x(i,j)/sqrt(grad_x(i,j)^2+grad_y(i,j)^2);
                ny=grad_y(i,j)/sqrt(grad_x(i,j)^2+grad_y(i,j)^2);
                nxx=abs(nx);


                //interface
                //_pointsupl+=2;
                if (p<=ny/(2*nxx))
                {
                    typinterf(i,j)=3;  //triangle vers la droite

                    _interface(k,1)=sqrt(2*p*ny/nxx);
                    _interface(k,2)=0;
                    _interface(k,3)=0;
                    _interface(k,4)=2*p/a(1);

                    if (_interface(k,4)==1)then   // si l'un des nouveaux points tombe sur l'angle du carré, on ne le compte pas comme point supplémentaire
                    {
                        //_pointsupl-=1
                    }
                }
                else if (p>=1-ny/(2*nxx))
                {
                    typinterf(i,j)=5;  //pentagone vers la droite

                    _interface(k,1)=1;
                    _interface(k,2)=1-sqrt(2*(1-p)*nxx/ny);
                    _interface(k,3)=1-2*(1-p)/(1-a(1));
                    _interface(k,4)=1;

                    if (_interface(k,1)==1)then
                    {
                        //_pointsupl-=1   // si l'un des nouveaux points tombe sur l'angle du carré, on ne le compte pas comme point supplémentaire
                    }
                }
                else
                {
                    if (nxx<ny)
                    {
                        typinterf(i,j)=4; //quadrillatère vers la droite

                        _interface(k,1)=1;
                        _interface(k,2)=p-nxx/(2*ny);
                        _interface(k,3)=0;
                        _interface(k,4)=p+nxx/(2*ny);
                    }
                    else
                    {
                        typinterf(i,j)=-4; //quadrillatère vers le haut


                        _interface(k,1)=p+ny/(2*nxx);
                        _interface(k,2)=0;
                        _interface(k,3)=p-ny/(2*nxx);
                        _interface(k,4)=1;
                    }
                }
                if (nx<0)
                {
                    if (typinterf(i,j)==-4)then
                    {
                        typinterf(i,j)-=10;
                    }
                    else
                    {
                        typinterf(i,j)+=10;  //vers la gauche
                    }

                    _interface(k,2)=1-_interface(k,2);
                    _interface(k,4)=1-_interface(k,4);
                }
            }
          }
          else
          {
            //_ninterf(i,j)=0;
            typinterf(i,j)=0;
          }
        }
        return 0;
}


/*

// Sauvegarde la solution
void plic::SaveSol( int n)
{
	string name_file = _results + "/solution_" + std::to_string(n) + ".vtk";

  int nb_vert = (lon+1)*(lar+1) + pointsupl;  //nombre de points

  //assert((sol.size() == _triangles.size()) && "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  solution << "# vtk DataFile Version 3.0 " << endl;
  solution << "2D Unstructured Grid" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET UNSTRUCTURED_GRID" << endl;

  solution << "POINTS " << nb_vert << " float " << endl;   //ajouter des points sur l'interface en fonction du type d'interface
  for (int i = 0 ; i < lon+1 ; ++i)
  {
      for(int j=0;j<lar+1;j++)
      {
          solution << i*_dx << " " << j*_dx << " 0." << endl;

          if (_ninterf(i,j)!=0)then    // ajout des points ajoutés sur l'interface
          {
              if (_interface(_ninterf(i,j),1)!=1)then
              {
                  solution << i*_dx+_interface(_ninterf(i,j),1) << " " << j*_dx << " 0." << endl;
              }
              if (_interface(_ninterf(i,j),4)!=1)then
              {
                  solution << i*_dx << " " << j*_dx+ _interface(_ninterf(i,j),4)<< " 0." << endl;
              }
          }
      }
  }
  solution << endl;

  solution << "CELLS " << lon*(lar+1) << " " << lon*lar*4 << endl;
  for (int i = 0 ; i < lon ; ++i)
  {
      for (int j=0;j<lar;j++)
      {
        if (typinterf(i,j)==0)then
        {
            solution << 4 << " " << ((_squares(i,j)).GetVertices(),0) << " " << ((_squares(i,j)).GetVertices(),1)<< " " << ((_squares(i,j)).GetVertices(),2) << " " << ((_squares(i,j)).GetVertices(),3) << endl;
        }
        else
        {
            solution << abs(typinterf(i,j))%10 << " "   //AAAAAAAAAAAAAAAAAAAAAAAAAAAAH
            if (typinterf(i)(j)==3)then
            {
                solution << (lar+1)*(i)+j << " " <<
            }
        }
      }
  }
  solution << endl;

  solution << "CELL_TYPES " << lon*(lar+1) << endl;
  for (int i = 0 ; i < lon ; ++i)
  {
      for (int j = 0 ; j < lon ; ++j)
      {
          if (typinterf(i)(j)==0)then
          {
              solution << 8 << endl;
          }
          if (typinterf(i)(j)%10==3)then
          {
              solution << 5 << endl;
          }
          if (typinterf(i)(j)%10==5)then
          {
              solution << 7 << endl;
          }
          if (abs(typinterf(i)(j))%10==4)then
          {
              solution << 9 << endl;
          }
      }
  }
  solution << endl;

  solution << "CELL_DATA " << _squares.size() << endl;
  solution << "SCALARS sol float 1" << endl;
  solution << "LOOKUP_TABLE default" << endl;
	double sum=0;
  for (int i = 0 ; i < lon ; ++i)
  {
      for (int j=0;j<lar;++j)
      {
          if (_phi(i)(j)==0)then
          {
              solution << 0 << endl;
          }
          else if (_phi(i)(j)==1)then
          {
              solution << 1 << endl;
          }
          else   //cas de l'interface
          {

          }
      }
  }
  solution << endl;

	cout<<sqrt(sum)<<endl;
	solution.close();
}

*/
