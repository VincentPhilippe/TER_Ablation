#include "plic.h"
#include <fstream>
#include <iostream>

using namespace std;

//constructeur
plic::plic()
  {
    _interface.resize(1,9);
  }

  plic::~plic()
  {}

double plic::grad_x(const int i,const int j, const int k)
{
  if (i==0)
  {
    return _phi(i+1,jpp)/2;
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
    dx=_read_data->Get_dx();
    dz=_read_data->Get_dz();
    _ninterf=_recul->Get_ninterf();
    _phi=_recul->Get_C_solide();
    int lon=_phi.cols();
    int lar=_phi.rows();
    int _kmax=_recul->Get_nbinterface();
    _interface.resize(_kmax,4);
    _normal.resize(_kmax,2);
    tri.resize(3*lon,3); //arbitraire pour le moment, assez grand pour contenir tous les triangles
    quad.resize(3*lon,4);
    penta.resize(3*lon,5);
    pttri.resize(9*lon,3);
    ptquad.resize(12*lon,4);
    ptpenta.resize(15*lon,5);
    trivalcase.resize(3*lon);
    quadvalcase.resize(3*lon);
    pentvalcase.resize(3*lon);
    int num=0;
    int nbtri=0;
    int nbquad=0;
    int nbpenta=0;
    /*
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
    */
    k=0;
    _pointsupl=0;
    for (int i=1;i<lon;i++)
     {
        for (int j=1;j<lar;j++)
        {
            p=_phi(i,j);
            if ((p>0.) && (p<1.))   //si on est sur l'interface
            {
                k++;

                //Calcul du gradient
                nx=grad_x(i,j)/sqrt(grad_x(i,j)*grad_x(i,j)+grad_y(i,j)*grad_y(i,j));
                ny=grad_y(i,j)/sqrt(grad_x(i,j)*grad_x(i,j)+grad_y(i,j)*grad_y(i,j));
                nxx=abs(nx);

                _normal(k,0)=nx;
                _normal(k,1)=ny;


                //interface
                //_pointsupl+=2;
                if (p<=ny/(2*nxx))
                {
                    //typinterf(i,j)=3;  //triangle vers la droite
                    num+=2;
                    _interface(k,1)=dx*sqrt(2*p*ny/nxx);
                    _interface(k,2)=0;
                    _interface(k,3)=0;
                    _interface(k,4)=dz*2*p/_interface(k,1);

                    //on rentre les coordonnées des sommets du triangles


                    if (nx>0)
                    {
                      pttri(nbtri*3,0)=i*dx;
                      pttri(nbtri*3+1,0)=i*dx+_interface(k,1);
                      pttri(nbtri*3+2,0)=i*dx;
                    }
                    else //si orienté vers la gauche
                    {
                      pttri(nbtri*3,0)=(i+1)*dx;
                      pttri(nbtri*3+1,0)=i*dx+dx-_interface(k,1);
                      pttri(nbtri*3+2,0)=(i+1)*dx;
                    }
                    pttri(nbtri*3,1)=(j+1)*dz;
                    pttri(nbtri*3,2)=0.0;
                    pttri(nbtri*3+1,1)=(j+1)*dz;
                    pttri(nbtri*3+1,2)=0.0;
                    pttri(nbtri*3+2,1)=(j+1)*dz-_interface(k,4);
                    pttri(nbtri*3+2,2)=0.0;

                    //case opposée
                    if (nx>0)
                    {
                      ptpenta(nbpenta*5,0)=(i+1)*dx;
                      ptpenta(nbpenta*5+1,0)=i*dx+_interface(k,1);
                      ptpenta(nbpenta*5+2,0)=i*dx;
                      ptpenta(nbpenta*5+3,0)=i*dx;
                      ptpenta(nbpenta*5+4,0)=(i+1)*dx;
                    }
                    else //si orienté vers la gauche
                    {
                      ptpenta(nbpenta*5,0)=(i+1)*dx;
                      ptpenta(nbpenta*5+1,0)=i*dx+dx-_interface(k,1);
                      ptpenta(nbpenta*5+2,0)=(i+1)*dx;
                      ptpenta(nbpenta*5+3,0)=(i+1)*dx;
                      ptpenta(nbpenta*5+4,0)=(i)*dx;
                    }
                    ptpenta(nbpenta*5,1)=(j+1)*dz;
                    ptpenta(nbpenta*5,2)=0.0;
                    ptpenta(nbpenta*5+1,1)=(j+1)*dz;
                    ptpenta(nbpenta*5+1,2)=0.0;
                    ptpenta(nbpenta*5+2,1)=(j+1)*dz-_interface(k,4);
                    ptpenta(nbpenta*5+2,2)=0.0;
                    ptpenta(nbpenta*5+3,1)=(j)*dz;
                    ptpenta(nbpenta*5+3,2)=0.0;
                    ptpenta(nbpenta*5+4,1)=(j)*dz;
                    ptpenta(nbpenta*5+4,2)=0.0;

                    //on assigne les sommets au triangle

                    //1er essai : on définit plusieurs fois les mêmes points
                    //tri(nbtri,0)=;
                    trivalcase(nbtri)=1;
                    pentvalcase(nbpenta)=0;
                    nbtri++;
                    nbpenta++;


                    if (_interface(k,4)==1)   // si l'un des nouveaux points tombe sur l'angle du carré, on ne le compte pas comme point supplémentaire
                    {
                        //_pointsupl-=1
                        cout<<"point sur un angle"<<endl;
                    }
                }
                else if (p>=1-ny/(2*nxx))
                {
                    typinterf(i,j)=5;  //pentagone vers la droite

                    _interface(k,1)=dx*1;
                    _interface(k,2)=dz*(1-sqrt(2*(1-p)*nxx/ny));
                    _interface(k,3)=dx*(1-2*(1-p)/(1-(_interface(k,2))/dz));
                    _interface(k,4)=dz*1;

                    //on rentre les coordonnées des sommets


                    if (nx>0)
                    {
                      ptpenta(nbpenta*5,0)=(i+1)*dx;
                      ptpenta(nbpenta*5+1,0)=i*dx+_interface(k,1);
                      ptpenta(nbpenta*5+2,0)=i*dx+_interface(k,3);
                      ptpenta(nbpenta*5+3,0)=i*dx;
                      ptpenta(nbpenta*5+4,0)=i*dx;
                    }
                    else //si orienté vers la gauche
                    {
                      ptpenta(nbpenta*5,0)=(i)*dx;
                      ptpenta(nbpenta*5+1,0)=i*dx+dx-_interface(k,1);
                      ptpenta(nbpenta*5+2,0)=i*dx+dx-_interface(k,3);
                      ptpenta(nbpenta*5+3,0)=(i+1)*dx;
                      ptpenta(nbpenta*5+4,0)=(i+1)*dx;
                    }
                    ptpenta(nbpenta*5,1)=(j+1)*dz;
                    ptpenta(nbpenta*5,2)=0.0;
                    ptpenta(nbpenta*5+1,1)=(j+1)*dz-_interface(k,2);
                    ptpenta(nbpenta*5+1,2)=0.0;
                    ptpenta(nbpenta*5+2,1)=j*dz;
                    ptpenta(nbpenta*5+2,2)=0.0;
                    ptpenta(nbpenta*5+3,1)=j*dz;
                    ptpenta(nbpenta*5+3,2)=0.0;
                    ptpenta(nbpenta*5+4,1)=(j+1)*dz;
                    ptpenta(nbpenta*5+4,2)=0.0;


                    //case opposée
                    if (nx>0)
                    {
                      pttri(nbtri*3,0)=(i+1)*dx;
                      pttri(nbtri*3+1,0)=(i+1)*dx;
                      pttri(nbtri*3+2,0)=i*dx+_interface(k,3);
                    }
                    else //si orienté vers la gauche
                    {
                      pttri(nbtri*3,0)=(i)*dx;
                      pttri(nbtri*3+1,0)=i*dx;
                      pttri(nbtri*3+2,0)=i*dx+dx-_interface(k,3);
                    }
                    pttri(nbtri*3,1)=(j)*dz;
                    pttri(nbtri*3,2)=0.0;
                    pttri(nbtri*3+1,1)=(j+1)*dz-_interface(k,2);
                    pttri(nbtri*3+1,2)=0.0;
                    pttri(nbtri*3+2,1)=j*dz;
                    pttri(nbtri*3+2,2)=0.0;


                    pentvalcase(nbpenta)=1;
                    trivalcase(nbtri)=0;
                    nbpenta++;
                    nbtri++;

                    if (_interface(k,1)==1)
                    {
                        //_pointsupl-=1   // si l'un des nouveaux points tombe sur l'angle du carré, on ne le compte pas comme point supplémentaire
                    }
                }
                else
                {
                    if (nxx<ny)
                    {
                        typinterf(i,j)=4; //quadrillatère vers le haut

                        _interface(k,1)=dx*1;
                        _interface(k,2)=dz*(p-nxx/(2*ny));
                        _interface(k,3)=0;
                        _interface(k,4)=dz*(p+nxx/(2*ny));
                        //on rentre les coordonnées des sommets
                        if (nx>0)
                        {
                          ptquad(nbquad*4,0)=i*dx;
                          ptquad(nbquad*4+1,0)=(i+1)*dx;
                          ptquad(nbquad*4+2,0)=(i+1)*dx;
                          ptquad(nbquad*4+3,0)=i*dx;
                        }
                        else //si orienté vers la gauche
                        {
                          ptquad(nbquad*4,0)=(i+1)*dx;
                          ptquad(nbquad*4+1,0)=(i)*dx;
                          ptquad(nbquad*4+2,0)=(i)*dx;
                          ptquad(nbquad*4+3,0)=(i+1)*dx;
                        }
                        ptquad(nbquad*4,1)=(j+1)*dz;
                        ptquad(nbquad*4,2)=0.0;
                        ptquad(nbquad*4+1,1)=(j+1)*dz;
                        ptquad(nbquad*4+1,2)=0.0;
                        ptquad(nbquad*4+2,1)=(j+1)*dz-_interface(k,2);
                        ptquad(nbquad*4+2,2)=0.0;
                        ptquad(nbquad*4+3,1)=(j+1)*dz-_interface(k,4);
                        ptquad(nbquad*4+3,2)=0.0;

                        //case opposée
                        if (nx>0)
                        {
                          ptquad(nbquad*4,0)=i*dx;
                          ptquad(nbquad*4+1,0)=(i+1)*dx;
                          ptquad(nbquad*4+2,0)=(i+1)*dx;
                          ptquad(nbquad*4+3,0)=i*dx;
                        }
                        else //si orienté vers la gauche
                        {
                          ptquad(nbquad*4,0)=(i+1)*dx;
                          ptquad(nbquad*4+1,0)=(i)*dx;
                          ptquad(nbquad*4+2,0)=(i)*dx;
                          ptquad(nbquad*4+3,0)=(i+1)*dx;
                        }
                        ptquad(nbquad*4,1)=(j)*dz;
                        ptquad(nbquad*4,2)=0.0;
                        ptquad(nbquad*4+1,1)=(j)*dz;
                        ptquad(nbquad*4+1,2)=0.0;
                        ptquad(nbquad*4+2,1)=(j+1)*dz-_interface(k,2);
                        ptquad(nbquad*4+2,2)=0.0;
                        ptquad(nbquad*4+3,1)=(j+1)*dz-_interface(k,4);
                        ptquad(nbquad*4+3,2)=0.0;

                    }
                    else
                    {
                        typinterf(i,j)=-4; //quadrillatère vers la droite

                        _interface(k,1)=dx*(p+ny/(2*nxx));
                        _interface(k,2)=0;
                        _interface(k,3)=dx*(p-ny/(2*nxx));
                        _interface(k,4)=dz*1;

                        //on rentre les coordonnées des sommets
                        if (nx>0)
                        {
                          ptquad(nbquad*4,0)=i*dx;
                          ptquad(nbquad*4+1,0)=i*dx+_interface(k,1);
                          ptquad(nbquad*4+2,0)=i*dx+_interface(k,3);
                          ptquad(nbquad*4+3,0)=i*dx;
                        }
                        else //si orienté vers la gauche
                        {
                          ptquad(nbquad*4,0)=(i+1)*dx;
                          ptquad(nbquad*4+1,0)=(i+1)*dx-_interface(k,1);
                          ptquad(nbquad*4+2,0)=(i+1)*dx-_interface(k,3);
                          ptquad(nbquad*4+3,0)=(i+1)*dx;
                        }
                        ptquad(nbquad*4,1)=(j+1)*dz;
                        ptquad(nbquad*4,2)=0.0;
                        ptquad(nbquad*4+1,1)=(j+1)*dz;
                        ptquad(nbquad*4+1,2)=0.0;
                        ptquad(nbquad*4+2,1)=(j)*dz;
                        ptquad(nbquad*4+2,2)=0.0;
                        ptquad(nbquad*4+3,1)=(j)*dz;
                        ptquad(nbquad*4+3,2)=0.0;

                        //case opposée
                        if (nx>0)
                        {
                          ptquad(nbquad*4,0)=(i+1)*dx;
                          ptquad(nbquad*4+1,0)=i*dx+_interface(k,1);
                          ptquad(nbquad*4+2,0)=i*dx+_interface(k,3);
                          ptquad(nbquad*4+3,0)=(i+1)*dx;
                        }
                        else //si orienté vers la gauche
                        {
                          ptquad(nbquad*4,0)=(i)*dx;
                          ptquad(nbquad*4+1,0)=(i+1)*dx-_interface(k,1);
                          ptquad(nbquad*4+2,0)=(i+1)*dx-_interface(k,3);
                          ptquad(nbquad*4+3,0)=(i)*dx;
                        }
                        ptquad(nbquad*4,1)=(j+1)*dz;
                        ptquad(nbquad*4,2)=0.0;
                        ptquad(nbquad*4+1,1)=(j+1)*dz;
                        ptquad(nbquad*4+1,2)=0.0;
                        ptquad(nbquad*4+2,1)=(j)*dz;
                        ptquad(nbquad*4+2,2)=0.0;
                        ptquad(nbquad*4+3,1)=(j)*dz;
                        ptquad(nbquad*4+3,2)=0.0;
                    }
                    quadvalcase(nbquad)=1;
                    nbquad++;
                    quadvalcase(nbquad)=0;
                    nbquad++;
                }
                if (nx<0) //si orienté vers la gauche
                {
                    if (typinterf(i,j)==-4)
                    {
                        typinterf(i,j)-=10;
                    }
                    else
                    {
                        typinterf(i,j)+=10;  //vers la gauche
                    }

                    _interface(k,2)=dx-_interface(k,2);
                    _interface(k,4)=dz-_interface(k,4);
                }
            }
            else
            {
              typinterf(i,j)=0;
              ptquad(nbquad*4,0)=(i+1)*dx;
              ptquad(nbquad*4,1)=(j+1)*dz;
              ptquad(nbquad*4,2)=0.0;
              ptquad(nbquad*4+1,0)=(i+1)*dx;
              ptquad(nbquad*4+1,1)=(j)*dz;
              ptquad(nbquad*4+1,2)=0.0;
              ptquad(nbquad*4+2,0)=(i)*dx;
              ptquad(nbquad*4+2,1)=(j)*dz;
              ptquad(nbquad*4+2,2)=0.0;
              ptquad(nbquad*4+3,0)=(i)*dx;
              ptquad(nbquad*4+3,1)=(j+1)*dz;
              ptquad(nbquad*4+3,2)=0.0;

              quadvalcase(nbquad)=0;
              nbquad++;
            }
          }
        }
        pttri.resize(nbtri*3,3);
        ptquad.resize(nbquad*4,3);
        ptpenta.resize(nbpenta*5,3);
}




// Sauvegarde la solution
void plic::SaveSol( int n)
{
	string name_file = "Results/solution_" + std::to_string(n) + ".vtk";

  int nb_vert = nbtri*3+nbquad*4+nbpenta*5;  //nombre de points

  //assert((sol.size() == _triangles.size()) && "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  solution << "# vtk DataFile Version 3.0 " << endl;
  solution << "2D Unstructured Grid" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET UNSTRUCTURED_GRID" << endl;

  solution << "POINTS " << nb_vert << " float " << endl;   //ajouter des points sur l'interface en fonction du type d'interface
  /*
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
  */
  for (int i=0;i<nbtri*3;++i)
  {
    solution << pttri(i,0)<<" "<<pttri(i,1)<<" "<<pttri(i,2)<<endl;
  }
  for (int i=0;i<nbquad*4;++i)
  {
    solution << ptquad(i,0)<<" "<<ptquad(i,1)<<" "<<ptquad(i,2)<<endl;
  }
  for (int i=0;i<nbpenta*4;++i)
  {
    solution << ptpenta(i,0)<<" "<<ptpenta(i,1)<<" "<<ptpenta(i,2)<<endl;
  }
  solution << endl;


  solution << "CELLS " << nbtri+nbquad+nbpenta << " " << (nbtri*4+nbquad*5+nbpenta*6) << endl; //deuxieme terme kesako ?
  for (int i;i<nbtri;++i)
  {
    solution << 3 << " "<< i*3 <<" "<<i*3 +1<<" "<<i*3 +2<<endl;
  }
  for (int i;i<nbquad;++i)
  {
    solution << 4 << " "<< nbtri+i*4 <<" "<<nbtri+i*4 +1<<" "<<nbtri+i*4 +2<<" "<<nbtri+i*4 +3<<endl;
  }
  for (int i;i<nbpenta;++i)
  {
    solution << 5 << " "<< nbtri+nbquad+i*5 <<" "<<nbtri+nbquad+i*5 +1<<" "<<nbtri+nbquad+i*5 +2<<" "<<nbtri+nbquad+i*5 +3<<" "<<nbtri+nbquad+i*5 +4<<endl;
  }

/*
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
  */
  solution << endl;

  solution << "CELL_TYPES " << nbtri+nbquad+nbpenta << endl;
  for (int i;i<nbtri;++i)
  {
    solution << 5<<endl;
  }
  for (int i;i<nbquad;++i)
  {
    solution << 9<<endl;
  }
  for (int i;i<nbpenta;++i)
  {
    solution << 7<<endl;
  }

  /*
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
  */
  solution << endl;

  solution << "CELL_DATA " << nbtri+nbquad+nbpenta << endl;
  solution << "SCALARS sol float 1" << endl;
  solution << "LOOKUP_TABLE default" << endl;
	//double sum=0;
  for (int i;i<nbtri;i++)
  {
    solution << trivalcase(i)<<endl;
  }
  for (int i;i<nbquad;i++)
  {
    solution << quadvalcase(i)<<endl;
  }
  for (int i;i<nbpenta;i++)
  {
    solution << pentvalcase(i)<<endl;
  }
  /*
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
  */
  solution << endl;

	//cout<<sqrt(sum)<<endl;
	solution.close();

}
