#include "recul3D.h"

using namespace std;
using namespace Eigen;

recul3D::recul3D(double dt, double dx, double dy, double dz, int nx, int ny, int nz, vector<vector<vector<double> > > C_solide)
{
  //_read_data=&_data;
  //_diff=0;
  _dtmax=dt;
  _dt=_dtmax;
  _dx=dx;
  _dy=dy;
  _dz=dz;
  _nx=nx;
  _ny=ny;
  _nz=nz;
  _C_solide=C_solide;
}


recul3D::~recul3D()
{}



  VectorXd recul3D::eqplan(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc)
  {
    MatrixXd A;
    VectorXd vect;
    A.resize(4,4);
    vect.resize(4);
    A << xa,ya,za,1,  xb,yb,zb,1,  xc,yc,zc,1,  1,1,1,1;
    vect << 0, 0, 0, 1;
    VectorXd sol = A.colPivHouseholderQr().solve(vect);

    return sol;
  }

  double recul3D::volume_pyramide(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd)//MatrixXd coord)
  {
    VectorXd sol;
    sol.resize(4);
    sol=eqplan(xa,ya,za,xb,yb,zb,xc,yc,zc);
    cout << sol << endl;
    double a,b,c,d,e, distance;
    a=sol(0);
    b=sol(1);
    c=sol(2);
    d=sol(3);
    e=-(a*xd+b*yd+c*zd+d)/(a*a+b*b+c*c);
    distance=abs(e)*sqrt(a*a+b*b+c*c);

    double volume, surface;
    surface=surface_triangle(xa,ya,za,xb,yb,zb,xc,yc,zc);
    cout << surface << endl;
    volume=surface*distance/3;

    return volume;

  }

  double recul3D::surface_triangle(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc)
  {
    double alpha, surface, la, lb, lc;
    la=sqrt((xc-xb)*(xc-xb)+(yc-yb)*(yc-yb)+(zc-zb)*(zc-zb));
    lb=sqrt((xc-xa)*(xc-xa)+(yc-ya)*(yc-ya)+(zc-za)*(zc-za));
    lc=sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb)+(za-zb)*(za-zb));
    alpha=acos((lb*lb+lc*lc-la*la)/(2*lb*lc));
    surface=lb*lc*sin(alpha)/2;

    return surface;
  }

  //cas 7 coin solide tous reste dans E
  void recul3D::recul3D_1(MatrixXd repere, MatrixXd coord, double vrdt)
  {
    double il, jl, kl;
    //attention repère reduit à D, E, G, H, M, O, Q, R
    il=repere(1,0);//attention pas la bonne ligne
    jl=repere(1,1);
    kl=repere(1,2);

    double xa1,ya1,za1,xb1,yb1,zb1,xc1,yc1,zc1;
    double xa2,ya2,za2,xb2,yb2,zb2,xc2,yc2,zc2;
    xa1=coord(0,0);
    ya1=coord(0,1);
    za1=coord(0,2);
    xb1=coord(1,0);
    yb1=coord(1,1);
    zb1=coord(1,2);
    xc1=coord(2,0);
    yc1=coord(2,1);
    zc1=coord(2,2);
    xa2=coord(3,0);
    ya2=coord(3,1);
    za2=coord(3,2);
    xb2=coord(4,0);
    yb2=coord(4,1);
    zb2=coord(4,2);
    xc2=coord(5,0);
    yc2=coord(5,1);
    zc2=coord(5,2);

    double surf, voltot;
    surf=surface_triangle(xa1,ya1,za1,xb1,yb1,zb1,xc1,yc1,zc1);
    voltot=surf*vrdt;

    


    _C_solide[il][jl][kl]-=voltot/(_dx*_dy*_dz);
  }

  MatrixXd recul3D::repereglobal(int i, int j, int k)
  {
    MatrixXd repere;
    repere.resize(18,3);
    repere.row(0) << i-1,j-1,k;//A
    repere.row(1) << i,j-1,k;//B
    repere.row(2) << i+1,j-1,k;//C
    repere.row(3) << i-1,j,k;//D
    repere.row(4) << i,j,k;//E
    repere.row(5) << i+1,j,k;//F
    repere.row(6) << i-1,j+1,k;//G
    repere.row(7) << i,j+1,k;//H
    repere.row(8) << i+1,j+1,k;//I
    repere.row(9) << i-1,j-1,k+1;//J
    repere.row(10) << i,j-1,k+1;//K
    repere.row(11) << i+1,j-1,k+1;//L
    repere.row(12) << i-1,j,k+1;//M
    repere.row(13) << i,j,k+1;//O
    repere.row(14) << i+1,j,k+1;//P
    repere.row(15) << i-1,j+1,k+1;//Q
    repere.row(16) << i,j+1,k+1;//R
    repere.row(17) << i+1,j+1,k+1;//S

    return repere;
  }

  MatrixXd recul3D::rotationz(MatrixXd repere_prec)
  {
    MatrixXd repere_suiv;
    repere_suiv.resize(18,3);
    repere_suiv.row(0)=repere_prec.row(2);
    repere_suiv.row(1)=repere_prec.row(5);
    repere_suiv.row(2)=repere_prec.row(8);
    repere_suiv.row(3)=repere_prec.row(1);
    repere_suiv.row(4)=repere_prec.row(4);
    repere_suiv.row(5)=repere_prec.row(7);
    repere_suiv.row(6)=repere_prec.row(0);
    repere_suiv.row(7)=repere_prec.row(3);
    repere_suiv.row(8)=repere_prec.row(6);
    repere_suiv.row(9)=repere_prec.row(11);
    repere_suiv.row(10)=repere_prec.row(14);
    repere_suiv.row(11)=repere_prec.row(17);
    repere_suiv.row(12)=repere_prec.row(10);
    repere_suiv.row(13)=repere_prec.row(13);
    repere_suiv.row(14)=repere_prec.row(16);
    repere_suiv.row(15)=repere_prec.row(9);
    repere_suiv.row(16)=repere_prec.row(12);
    repere_suiv.row(17)=repere_prec.row(15);

    return repere_suiv;
  }

  MatrixXd recul3D::reductionrepere(MatrixXd repere_prec)
  {
    MatrixXd repere_suiv;
    repere_suiv.resize(8,3);
    repere_suiv.row(0)=repere_prec.row(3);//D
    repere_suiv.row(1)=repere_prec.row(4);//E
    repere_suiv.row(2)=repere_prec.row(6);//G
    repere_suiv.row(3)=repere_prec.row(7);//H
    repere_suiv.row(4)=repere_prec.row(12);//M
    repere_suiv.row(5)=repere_prec.row(13);//O
    repere_suiv.row(6)=repere_prec.row(15);//Q
    repere_suiv.row(7)=repere_prec.row(16);//R

    return repere_suiv;
  }

  MatrixXd recul3D::rotationcoin(MatrixXd repere_prec)
  {
    MatrixXd repere_suiv;
    repere_suiv.resize(8,3);
    repere_suiv.row(0)=repere_prec.row(5);//O->D
    repere_suiv.row(1)=repere_prec.row(1);//E->E
    repere_suiv.row(2)=repere_prec.row(4);//M->G
    repere_suiv.row(3)=repere_prec.row(0);//D->H
    repere_suiv.row(4)=repere_prec.row(7);//R->M
    repere_suiv.row(5)=repere_prec.row(3);//H->O
    repere_suiv.row(6)=repere_prec.row(6);//Q->Q
    repere_suiv.row(7)=repere_prec.row(2);//G->R

    return repere_suiv;
  }
