#include "recul.h"

using namespace std;
using namespace Eigen;

recul::recul(read_data &_data)
:_read_data(_data)//, _diff(0)
{
  //_read_data=&_data;
  //_diff=0;
  _C_solide=_read_data.Get_C_Solide();
  cout << _C_solide << endl;
  _dtmax=_read_data.Get_dt();
  _dt=_dtmax;
  _dx=_read_data.Get_dx();
  _dz=_read_data.Get_dz();
  _nx = _C_solide.cols();
  _nz = _C_solide.rows();
  _ninterf.resize(_nz,_nx);
}


recul::~recul()
{}

void recul::recul_surface()
{
  //récuperer _ninterf _interface _vitesse avec des arguments
  //_ninterf=_plic->Get_ninterf();
  _interface=_plic->Get_interface();
  _vitesse=_diff->GetVitesse();
  int affichage;
  affichage=0;

  //cout << "_C_solide" << endl << _C_solide << endl;
  if (affichage==1) {
    cout << "_ninterf" << endl << _ninterf << endl;
    cout << "_interface" << endl << _interface << endl;
    cout << "_vitesse" << endl << _vitesse << endl;
  }

  double maxvr;
  maxvr=_vitesse.maxCoeff();
  if ((maxvr*_dtmax<_dx) & (maxvr*_dtmax<_dz)) {
    _dt=_dtmax;
  } else {
    _dt=min(_dx,_dz)/maxvr;
  }

  //boucle sur les surfaces
  for(int i=0; i<_nz; i++)
  {
    for(int j=0; j<_nx; j++)
    {
      //cout << "i " << i << " j " << j << endl;
      double k;

      k = _ninterf(i,j);

      if (k>0) {
        //cout << "k " << k << endl;
        //calcul de l'angle alpha
        double xa, za, xb, zb, t_alpha, alpha;

        //repere local
        if (_interface(k-1,0)<_interface(k-1,2)) {
          xa = _interface(k-1,0);
          za = _interface(k-1,1);
          xb = _interface(k-1,2);
          zb = _interface(k-1,3);
        } else {
          xb = _interface(k-1,0);
          zb = _interface(k-1,1);
          xa = _interface(k-1,2);
          za = _interface(k-1,3);
        }

        /*xa = -j*_dx + _interface(0,k-1);
        za = (i+1)*_dz - _interface(1,k-1);
        xb = -j*_dx + _interface(2,k-1);
        zb = (i+1)*_dz - _interface(3,k-1);*/
        t_alpha = (zb-za)/(xb-xa);
        alpha = atan(t_alpha);

        //calcul des coordonnées des points C et D
        double xc, zc, xd, zd, vr, vrdt;

        vr = _vitesse(k-1);
        vrdt = vr*_dt;
        xc = xa + vrdt*sin(alpha);
        zc = za - vrdt*cos(alpha);
        xd = xb + vrdt*sin(alpha);
        zd = zb - vrdt*cos(alpha);

        /*MatrixXd coord;
        coord.resize(4,2);
        coord(0,0)=xa;
        coord(0,1)=za;
        coord(1,0)=xb;
        coord(1,1)=zb;
        coord(2,0)=xc;
        coord(2,1)=zc;
        coord(3,0)=xd;
        coord(3,1)=zd;*/

        double l;
        l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
        alpha=abs(alpha);

        //identification du cas et modification du tableau des concentrations en solide
        if (xc<0) {
          if (xd<0) {
            if (zc<0) {//yd forcément négatif
              //recul3(i, j, alpha, vrdt, coord);
              if (affichage==1) {
                cout << "cas3" << endl;
              }
              double Stot, S1, S2, S3, S4;
              Stot=vrdt*l;
              S1 = za*za*tan(alpha)/2;
              S4 = za*xb/2;
              S3 = xb*xb/(2*tan(alpha));
              S2 = Stot-(S1+S3+S4);
              _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
              _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
              if (i+1<_nz) {
                _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S2/(_dx*_dz);
                _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
              }
            } else {//yc>0
              if (zd<0) {
                //recul2(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas2" << endl;
                }
                double Stot, S1, S2, S3, S4;
                Stot=vrdt*l;
                S3 = xb*xb/(2*tan(alpha));
                S1 = -xc*zc-zc*zc/(2*tan(alpha))+xc*xc/(2*tan(alpha));
                S2 = xd*zd-xd*xd/(2*tan(alpha))+zd*zd/(2*tan(alpha));
                S4 = Stot-(S1+S2+S3);
                _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
                _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
                if (i+1<_nz) {
                  _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S2/(_dx*_dz);
                  _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
                }
              } else {//yd>0
                //recul1(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas1" << endl;
                }
                double Stot, S1, S2, l1, l2;
                Stot=vrdt*l;
                l1 = -xc*sin(alpha);
                l2 = -xd*sin(alpha);
                S1 = l*(l1+l2)/2;
                S2 = Stot-S1;
                _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
                _C_solide(i,j)=_C_solide(i,j)-(S2)/(_dx*_dz);
              }
            }
          } else {//0<xd<dx
            if (zc<0) {//yd forcémentnégatif
              //recul9(i, j, alpha, vrdt, coord);
              if (affichage==1) {
                cout << "cas9" << endl;
              }
              double Stot, S1, S2, S3, S4;
              Stot=vrdt*l;
              S1 = za*za*tan(alpha)/2;
              S2 = xc*zc-zc*zc*tan(alpha)/2+xc*xc*tan(alpha)/2;
              S3 = -xd*zd-xd*xd*tan(alpha)/2+zd*zd*tan(alpha)/2;
              S4 = Stot-(S1+S2+S3);

              _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
              _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
              if (i+1<_nz) {
                _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S2/(_dx*_dz);
                _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
              }
            } else {//yc>0
              if (zd<0) {
                if (zc-xc*(zd-zc)/(xd-xc)>0) {
                  //recul7(i, j, alpha, vrdt, coord);
                  if (affichage==1) {
                    cout << "cas7" << endl;
                  }
                  double S1, S2;

                  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
                  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
                  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
                  if (i+1<_nz) {
                    _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
                  }
                  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2)/(_dx*_dz);
                } else {
                  //recul8(i, j, alpha, vrdt, coord);
                  if (affichage==1) {
                    cout << "cas8" << endl;
                  }
                  double S1,S2,S3;
                  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
                  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
                  S3=(-xc*tan(alpha)-zc)*(-xc*tan(alpha)-zc)/(2*tan(alpha));

                  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-(S1-S3)/(_dx*_dz);
                  if (i+1<_nz) {
                    _C_solide(i+1,j)=_C_solide(i+1,j)-(S2-S3)/(_dx*_dz);
                    _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S3/(_dx*_dz);
                  }
                  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2+S3)/(_dx*_dz);
                }
              } else {//yd>0
                //recul4(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas4" << endl;
                }
                double S1;

                S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
                _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
                _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1)/(_dx*_dz);
              }
            }
          }
        } else if (xc<_dx) {
          if (xd<_dx) {
            if (zc<0) {
              if (zd<0) {
                //recul11(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas11" << endl;
                }
                double Stot, S1, S2;
                Stot=vrdt*l;
                S1=-zc*l;
                S2=Stot-S1;

                _C_solide(i,j)=_C_solide(i,j)-(S2)/(_dx*_dz);
                if (i+1<_nz) {
                  _C_solide(i+1,j)=_C_solide(i+1,j)-S1/(_dx*_dz);
                }
              } else {
                //recul6(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas6" << endl;
                }
                double S2;

                S2=zc*zc*(1/tan(alpha)+tan(alpha))/2;
                if (i+1<_nz) {
                  _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
                }
                _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S2)/(_dx*_dz);
              }
            } else {
              if (zd<0) {
                //recul10(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas10" << endl;
                }
                double S2;

                S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
                if (i+1<_nz) {
                  _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
                }
                _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S2)/(_dx*_dz);
              } else {
                //recul5(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas5" << endl;
                }
                _C_solide(i,j)=_C_solide(i,j)-l*vrdt/(_dx*_dz);
              }
            }
          } else {
            if (zc<0) {
              if (zd<0) {
                //recul17(i, j, alpha, vrdt, coord);
                if (affichage==1) {
                  cout << "cas17" << endl;
                }
                double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
                xb1=_dx-xa;
                zb1=za;
                xa1=_dx-xb;
                za1=zb;
                xd1=_dx-xc;
                zd1=zc;
                xc1=_dx-xd;
                zc1=zd;
                xa=xa1;
                za=za1;
                xb=xb1;
                zb=zb1;
                xc=xc1;
                zc=zc1;
                xd=xd1;
                zd=zd1;

                double Stot, S1, S2, S3, S4;
                Stot=vrdt*l;
                S1 = za*za*tan(alpha)/2;
                S2 = xc*zc-zc*zc*tan(alpha)/2+xc*xc*tan(alpha)/2;
                S3 = -xd*zd-xd*xd*tan(alpha)/2+zd*zd*tan(alpha)/2;
                S4 = Stot-(S1+S2+S3);

                _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
                _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
                if (i+1<_nz) {
                  _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S2/(_dx*_dz);
                  _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
                }
              } else {
                if (zc+(_dx-xc)*(zd-zc)/(xd-xc)>0) {
                  //recul13(i, j, alpha, vrdt, coord);
                  if (affichage==1) {
                    cout << "cas13" << endl;
                  }
                  double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
                  xb1=_dx-xa;
                  zb1=za;
                  xa1=_dx-xb;
                  za1=zb;
                  xd1=_dx-xc;
                  zd1=zc;
                  xc1=_dx-xd;
                  zc1=zd;
                  xa=xa1;
                  za=za1;
                  xb=xb1;
                  zb=zb1;
                  xc=xc1;
                  zc=zc1;
                  xd=xd1;
                  zd=zd1;

                  double S1, S2;
                  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
                  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
                  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
                  if (i+1<_nz) {
                    _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
                  }
                  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2)/(_dx*_dz);
                } else {
                  //recul14(i, j, alpha, vrdt, coord);
                  if (affichage==1) {
                    cout << "cas14" << endl;
                  }
                  double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
                  xb1=_dx-xa;
                  zb1=za;
                  xa1=_dx-xb;
                  za1=zb;
                  xd1=_dx-xc;
                  zd1=zc;
                  xc1=_dx-xd;
                  zc1=zd;
                  xa=xa1;
                  za=za1;
                  xb=xb1;
                  zb=zb1;
                  xc=xc1;
                  zc=zc1;
                  xd=xd1;
                  zd=zd1;

                  double S1, S2, S3;
                  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
                  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
                  S3=(-xc*tan(alpha)-zc)*(-xc*tan(alpha)-zc)/(2*tan(alpha));

                  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-(S1-S3)/(_dx*_dz);
                  if (i+1<_nz) {
                    _C_solide(i+1,j)=_C_solide(i+1,j)-(S2-S3)/(_dx*_dz);
                    _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S3/(_dx*_dz);
                  }
                  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2+S3)/(_dx*_dz);
                }
              }
            } else {
              //recul12(i, j, alpha, vrdt, coord);
              if (affichage==1) {
                cout << "cas12" << endl;
              }
              double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
              xb1=_dx-xa;
              zb1=za;
              xa1=_dx-xb;
              za1=zb;
              xd1=_dx-xc;
              zd1=zc;
              xc1=_dx-xd;
              zc1=zd;
              xa=xa1;
              za=za1;
              xb=xb1;
              zb=zb1;
              xc=xc1;
              zc=zc1;
              xd=xd1;
              zd=zd1;

              double S1;
              if (abs(alpha)>0.000000000001) {
                S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
              } else {
                S1=0;
              }
              _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
              _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1)/(_dx*_dz);
            }
          }
        } else {
          if (zc<0) {
            if (zd<0) {
              //recul18(i, j, alpha, vrdt, coord);
              if (affichage==1) {
                cout << "cas18" << endl;
              }
              double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
              xb1=_dx-xa;
              zb1=za;
              xa1=_dx-xb;
              za1=zb;
              xd1=_dx-xc;
              zd1=zc;
              xc1=_dx-xd;
              zc1=zd;
              xa=xa1;
              za=za1;
              xb=xb1;
              zb=zb1;
              xc=xc1;
              zc=zc1;
              xd=xd1;
              zd=zd1;

              double Stot, S1, S2, S3, S4;
              Stot=vrdt*l;
              S1 = za*za*tan(alpha)/2;
              S4 = za*xb/2;
              S3 = xb*xb/(2*tan(alpha));
              S2 = Stot-(S1+S3+S4);
              _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
              _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
              if (i+1<_nz) {
                _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S2/(_dx*_dz);
                _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
              }
            } else {
              //recul16(i, j, alpha, vrdt, coord);
              if (affichage==1) {
                cout << "cas16" << endl;
              }
              double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
              xb1=_dx-xa;
              zb1=za;
              xa1=_dx-xb;
              za1=zb;
              xd1=_dx-xc;
              zd1=zc;
              xc1=_dx-xd;
              zc1=zd;
              xa=xa1;
              za=za1;
              xb=xb1;
              zb=zb1;
              xc=xc1;
              zc=zc1;
              xd=xd1;
              zd=zd1;

              double Stot, S1, S2, S3, S4;
              Stot=vrdt*l;
              S3 = xb*xb/(2*tan(alpha));
              S1 = -xc*zc-zc*zc/(2*tan(alpha))+xc*xc/(2*tan(alpha));
              S2 = xd*zd-xd*xd/(2*tan(alpha))+zd*zd/(2*tan(alpha));
              S4 = Stot-(S1+S2+S3);
              _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
              _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
              if (i+1<_nz) {
                _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S2/(_dx*_dz);
                _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
              }
            }
          } else {
            //recul15(i, j, alpha, vrdt, coord);
            if (affichage==1) {
              cout << "cas15" << endl;
            }
            double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
            xb1=_dx-xa;
            zb1=za;
            xa1=_dx-xb;
            za1=zb;
            xd1=_dx-xc;
            zd1=zc;
            xc1=_dx-xd;
            zc1=zd;
            xa=xa1;
            za=za1;
            xb=xb1;
            zb=zb1;
            xc=xc1;
            zc=zc1;
            xd=xd1;
            zd=zd1;

            double Stot, S1, S2, l1, l2;
            Stot=vrdt*l;
            l1 = -xc*sin(alpha);
            l2 = -xd*sin(alpha);
            S1 = l*(l1+l2)/2;
            S2 = Stot-S1;
            _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
            _C_solide(i,j)=_C_solide(i,j)-(S2)/(_dx*_dz);
          }
        }

      }
    }
  }

  cpositive();

  //écriture dans un fichier
  ofstream mon_flux;
  string name_file("C_solide.dat");
  mon_flux.open(name_file, ios::out);
  if(mon_flux)
  {
    mon_flux << _C_solide << endl;
    /*for (int i = 0; i < _nz; i++) {
      for (int j = 0; j < _nx; j++) {
        mon_flux <<  << endl;
      }
    }*/
  }
  else
  {
    cout << "ERREUR: Impossible d’ouvrir le fichier." << endl;
  }
  mon_flux.close();

}

/*
void recul::recul1(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double Stot, S1, S2, l1, l2;
  Stot=vrdt*l;
  l1 = -xc*sin(alpha);
  l2 = -xd*sin(alpha);
  S1 = l*(l1+l2)/2;
  S2 = Stot-S1;
  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S2)/(_dx*_dz);

}


void recul::recul2(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S3 = xb*xb/(2*tan(alpha));
  S1 = -xc*zc-zc*zc/(2*tan(alpha))+xc*xc/(2*tan(alpha));
  S2 = xd*zd-xd*xd/(2*tan(alpha))+zd*zd/(2*tan(alpha));
  S4 = Stot-(S1+S2+S3);
  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S2/(_dx*_dz);
    _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
  }

}


void recul::recul3(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S1 = za*za*tan(alpha)/2;
  S4 = za*xb/2;
  S3 = xb*xb/(2*tan(alpha));
  S2 = Stot-(S1+S3+S4);
  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S2/(_dx*_dz);
    _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
  }

}

//C à gauche cas 4
void recul::recul4(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double S1;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1)/(_dx*_dz);

}

//C et D dans Cij cas 5
void recul::recul5(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  _C_solide(i,j)=_C_solide(i,j)-l*vrdt/(_dx*_dz);

}


void recul::recul6(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double S2;

  S2=zc*zc*(1/tan(alpha)+tan(alpha))/2;
  if (i+1<_nz) {
    _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
  }
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S2)/(_dx*_dz);

}

//C à gauche et D en dessous cas 7
void recul::recul7(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double S1, S2;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
  }
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2)/(_dx*_dz);

}

//C à gauche et D en dessous, depasse dans le coin cas 8
void recul::recul8(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double S1, S2, S3;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
  S3=(-xc*tan(alpha)-zc)*(-xc*tan(alpha)-zc)/(2*tan(alpha));

  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-(S1-S3)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,j)=_C_solide(i+1,j)-(S2-S3)/(_dx*_dz);
    _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S3/(_dx*_dz);
  }
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2+S3)/(_dx*_dz);

}

void recul::recul9(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S1 = za*za*tan(alpha)/2;
  S2 = xc*zc-zc*zc*tan(alpha)/2+xc*xc*tan(alpha)/2;
  S3 = -xd*zd-xd*xd*tan(alpha)/2+zd*zd*tan(alpha)/2;
  S4 = Stot-(S1+S2+S3);

  _C_solide(i,(j-1+_nx)%_nx)=_C_solide(i,(j-1+_nx)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,(j-1+_nx)%_nx)=_C_solide(i+1,(j-1+_nx)%_nx)-S2/(_dx*_dz);
    _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
  }

}

//D en dessous cas 10
void recul::recul10(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double S2;

  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
  if (i+1<_nz) {
    _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
  }
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S2)/(_dx*_dz);

}

void recul::recul11(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  xa=coord(0,0);
  za=coord(0,1);
  xb=coord(1,0);
  zb=coord(1,1);
  xc=coord(2,0);
  zc=coord(2,1);
  xd=coord(3,0);
  zd=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S1=-zc*l;
  S2=Stot-S1;

  _C_solide(i,j)=_C_solide(i,j)-(S2)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,j)=_C_solide(i+1,j)-S1/(_dx*_dz);
  }
}

void recul::recul12(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  //coordonnées symétriques
  xb=_dx-coord(0,0);
  zb=coord(0,1);
  xa=_dx-coord(1,0);
  za=coord(1,1);
  xd=_dx-coord(2,0);
  zd=coord(2,1);
  xc=_dx-coord(3,0);
  zc=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  //symétrie cas 4
  double S1;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1)/(_dx*_dz);

}

void recul::recul13(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  //coordonnées symétriques
  xb=_dx-coord(0,0);
  zb=coord(0,1);
  xa=_dx-coord(1,0);
  za=coord(1,1);
  xd=_dx-coord(2,0);
  zd=coord(2,1);
  xc=_dx-coord(3,0);
  zc=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  //symétrie cas 7
  double S1, S2;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,j)=_C_solide(i+1,j)-S2/(_dx*_dz);
  }
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2)/(_dx*_dz);

}

void recul::recul14(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  //coordonnées symétriques
  xb=_dx-coord(0,0);
  zb=coord(0,1);
  xa=_dx-coord(1,0);
  za=coord(1,1);
  xd=_dx-coord(2,0);
  zd=coord(2,1);
  xc=_dx-coord(3,0);
  zc=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  //symétrie cas 8
  double S1, S2, S3;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  S2=zd*zd*(1/tan(alpha)+tan(alpha))/2;
  S3=(-xc*tan(alpha)-zc)*(-xc*tan(alpha)-zc)/(2*tan(alpha));

  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-(S1-S3)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,j)=_C_solide(i+1,j)-(S2-S3)/(_dx*_dz);
    _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S3/(_dx*_dz);
  }
  _C_solide(i,j)=_C_solide(i,j)-(l*vrdt-S1-S2+S3)/(_dx*_dz);

}

void recul::recul15(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  //coordonnées symétriques
  xb=_dx-coord(0,0);
  zb=coord(0,1);
  xa=_dx-coord(1,0);
  za=coord(1,1);
  xd=_dx-coord(2,0);
  zd=coord(2,1);
  xc=_dx-coord(3,0);
  zc=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  //symétrie cas 1
  double Stot, S1, S2, l1, l2;
  Stot=vrdt*l;
  l1 = -xc*sin(alpha);
  l2 = -xd*sin(alpha);
  S1 = l*(l1+l2)/2;
  S2 = Stot-S1;
  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S2)/(_dx*_dz);

}

void recul::recul16(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  //coordonnées symétriques
  xb=_dx-coord(0,0);
  zb=coord(0,1);
  xa=_dx-coord(1,0);
  za=coord(1,1);
  xd=_dx-coord(2,0);
  zd=coord(2,1);
  xc=_dx-coord(3,0);
  zc=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  //symétrie cas 2
  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S3 = xb*xb/(2*tan(alpha));
  S1 = -xc*zc-zc*zc/(2*tan(alpha))+xc*xc/(2*tan(alpha));
  S2 = xd*zd-xd*xd/(2*tan(alpha))+zd*zd/(2*tan(alpha));
  S4 = Stot-(S1+S2+S3);
  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S2/(_dx*_dz);
    _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
  }

}

void recul::recul17(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  //coordonnées symétriques
  xb=_dx-coord(0,0);
  zb=coord(0,1);
  xa=_dx-coord(1,0);
  za=coord(1,1);
  xd=_dx-coord(2,0);
  zd=coord(2,1);
  xc=_dx-coord(3,0);
  zc=coord(3,1);

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  //symétrie cas 9
  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S1 = za*za*tan(alpha)/2;
  S2 = xc*zc-zc*zc*tan(alpha)/2+xc*xc*tan(alpha)/2;
  S3 = -xd*zd-xd*xd*tan(alpha)/2+zd*zd*tan(alpha)/2;
  S4 = Stot-(S1+S2+S3);

  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S2/(_dx*_dz);
    _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
  }

}

void recul::recul18(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,za,xb,zb,xc,zc,xd,zd,l;

  //coordonnées symétriques
  xb=_dx-coord(0,0);
  zb=coord(0,1);
  xa=_dx-coord(1,0);
  za=coord(1,1);
  xd=_dx-coord(2,0);
  zd=coord(2,1);
  xc=_dx-coord(3,0);
  zc=coord(3,1);

  double xa1,za1,xb1,zb1,xc1,zc1,xd1,zd1;
  xb1=_dx-xa;
  zb1=za;
  xa1=_dx-xb;
  za1=zb;
  xd1=_dx-xc;
  zd1=zc;
  xc1=-dx-xd;
  zc1=zd;
  xa=xa1;
  za=za1;
  xb=xb1;
  zb=zb1;
  xc=xc1;
  zc=zc1;
  xd=xd1;
  zd=zd1;

  l=sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za));
  alpha=abs(alpha);

  //symétrie cas 3
  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S1 = za*za*tan(alpha)/2;
  S4 = za*xb/2;
  S3 = xb*xb/(2*tan(alpha));
  S2 = Stot-(S1+S3+S4);
  _C_solide(i,(j+1)%_nx)=_C_solide(i,(j+1)%_nx)-S1/(_dx*_dz);
  _C_solide(i,j)=_C_solide(i,j)-(S4)/(_dx*_dz);
  if (i+1<_nz) {
    _C_solide(i+1,(j+1)%_nx)=_C_solide(i+1,(j+1)%_nx)-S2/(_dx*_dz);
    _C_solide(i+1,j)=_C_solide(i+1,j)-S3/(_dx*_dz);
  }

}
*/

void recul::cpositive()
{
  _nbinterface=0;
  for (int i=0; i<_nz; i++) {
    for (int j = 0; j <_nx; j++) {
      if (_C_solide(i,j)<=0.00000000001) {
        _C_solide(i,j)=0;
        _ninterf(i,j)=0;
      } else if (_C_solide(i,j)>=1) {
        _ninterf(i,j)=-1;
      } else {
        _nbinterface=_nbinterface+1;
        _ninterf(i,j)=_nbinterface;
      }
    }
  }
}
