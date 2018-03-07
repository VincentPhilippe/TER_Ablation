#include "plic.h"
#include <iostream>

using namespace std;

double plic::grad_x(const i,const j)
{
  if (i==0)
  {
    return _phi[i+1]/2;
  }
  if (i==lon)
  {
    return _phi[i-1]/2;
  }
  else
  {
    return (_phi[i+1]+_phi[i-1])/2;
  }
}

double plic::grad_y(const i,const j)
{
  {
    if (j==0)
    {
      return _phi[j+1]/2;
    }
    if (j==lar)
    {
      return _phi[j-1]/2;
    }
    else
    {
      return (_phi[j+1]+_phi[j-1])/2;
    }
  }
}

void plic::interf(const int lon,const int lar)
{
    _ninterf.resize(lon,lar)
    _interface.resize(nb_inter,4)   //nb_inter = max(_ninter), à voir plus tard
    k=0
    for (int i=1;i<lon;i++)
     {
        for (int j=1;j<lar;j++)
        {
            p=_phi(i,j);
            if ((p>0.) && (p<1.))   //si on est sur l'interface
            {
                k++;
                _ninterf[i][j]=k;
                //Calcul du gradient
                nx=grad_x(i,j)/sqrt(grad_x(i,j)^2+grad_y(i,j)^2);
                ny=grad_y(i,j)/sqrt(grad_x(i,j)^2+grad_y(i,j)^2);
                nxx=abs(nx);


                //interface
                if (p<=ny/(2*nxx))
                {
                    _inteface[k][1]=sqrt(2*p*ny/nxx);
                    _inteface[k][2]=0;
                    _inteface[k][3]=0;
                    _inteface[k][4]=2*p/a(1);
                }
                else if (p>=1-ny/(2*nxx))
                {
                    _inteface[k][1]=1;
                    _inteface[k][2]=1-sqrt(2*(1-p)*nxx/ny);
                    _inteface[k][3]=1-2*(1-p)/(1-a(1));
                    _inteface[k][4]=1;
                }
                else
                {
                    if (nxx<ny)
                    {
                        _inteface[k][1]=1;
                        _inteface[k][2]=p-nxx/(2*ny);
                        _inteface[k][3]=0;
                        _inteface[k][4]=p+nxx/(2*ny);
                    }
                    else
                    {
                        _inteface[k][1]=p+ny/(2*nxx);
                        _inteface[k][2]=0;
                        _inteface[k][3]=p-ny/(2*nxx);
                        _inteface[k][4]=1;
                    }
                }
                if (nx<0)
                {
                    _inteface[k][2]=1-_inteface[k][2];
                    _inteface[k][4]=1-_inteface[k][4];
                }
            }
          }
          else
          {
            _ninterf[i][j]=0;
          }
        }
}


// Sauvegarde la solution
void plic::SaveSol( int n)
{
	string name_file = _results + "/solution_" + std::to_string(n) + ".vtk";

  int nb_vert = _vertices.size();

  //assert((sol.size() == _triangles.size()) && "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  solution << "# vtk DataFile Version 3.0 " << endl;
  solution << "2D Unstructured Grid" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET UNSTRUCTURED_GRID" << endl;

  solution << "POINTS " << nb_vert << " float " << endl;
  for (int i = 0 ; i < nb_vert ; ++i)
  {
    solution << ((_vertices[i]).GetCoor())[0] << " " << ((_vertices[i]).GetCoor())[1] << " 0." << endl;
  }
  solution << endl;

  solution << "CELLS " << _triangles.size() << " " << _triangles.size()*4 << endl;
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << 3 << " " << ((_triangles[i]).GetVertices())[0] << " " << ((_triangles[i]).GetVertices())[1]
    << " " << ((_triangles[i]).GetVertices())[2] << endl;
  }
  solution << endl;

  solution << "CELL_TYPES " << _triangles.size() << endl;
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << 5 << endl;
  }
  solution << endl;

  solution << "CELL_DATA " << _triangles.size() << endl;
  solution << "SCALARS sol float 1" << endl;
  solution << "LOOKUP_TABLE default" << endl;
	double sum=0;
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
		sum+=_tri_area(i)*sol[i]*sol[i];
    solution << float(sol[i]) << endl;
  }
  solution << endl;

	cout<<sqrt(sum)<<endl;
	solution.close();
}
