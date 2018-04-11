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



  VectorXd recul3D::eqplan(vector<double> pta, vector<double> ptb, vector<double> ptc)
  {
    MatrixXd A;
    VectorXd vect;
    A.resize(4,4);
    vect.resize(4);
    A << pta[0],pta[1],pta[2],1,  ptb[0],ptb[1],ptb[2],1,  ptc[0],ptc[1],ptc[2],1,  1,1,1,1;
    vect << 0, 0, 0, 1;
    VectorXd sol = A.colPivHouseholderQr().solve(vect);

    return sol;
  }

  double recul3D::volume_pyramide(vector<double> pta, vector<double> ptb, vector<double> ptc, vector<double> ptd)//MatrixXd coord)
  {
    VectorXd sol;
    sol.resize(4);
    sol=eqplan(pta,ptb,ptc);
    cout << sol << endl;
    double a,b,c,d,e, distance;
    a=sol(0);
    b=sol(1);
    c=sol(2);
    d=sol(3);
    e=-(a*ptd[0]+b*ptd[1]+c*ptd[2]+d)/(a*a+b*b+c*c);
    distance=abs(e)*sqrt(a*a+b*b+c*c);

    double volume, surface;
    surface=surface_triangle(pta,ptb,ptc);
    cout << surface << endl;
    volume=surface*distance/3;

    return volume;

  }

  double recul3D::surface_triangle(vector<double> pta, vector<double> ptb, vector<double> ptc)
  {
    double alpha, surface, la, lb, lc;
    la=sqrt((ptc[0]-ptb[0])*(ptc[0]-ptb[0])+(ptc[1]-ptb[1])*(ptc[1]-ptb[1])+(ptc[2]-ptb[2])*(ptc[2]-ptb[2]));
    lb=sqrt((ptc[0]-pta[0])*(ptc[0]-pta[0])+(ptc[1]-pta[1])*(ptc[1]-pta[1])+(ptc[2]-pta[2])*(ptc[2]-pta[2]));
    lc=sqrt((pta[0]-ptb[0])*(pta[0]-ptb[0])+(pta[1]-ptb[1])*(pta[1]-ptb[1])+(pta[2]-ptb[2])*(pta[2]-ptb[2]));
    alpha=acos((lb*lb+lc*lc-la*la)/(2*lb*lc));
    surface=lb*lc*sin(alpha)/2;

    return surface;
  }

  //cas 7 coin solide tous reste dans E
  void recul3D::recul3D_1(MatrixXd& repere, vector<vector<double>>& coord, double vrdt)
  {
    double il, jl, kl;
    //attention repère reduit à D, E, G, H, M, O, Q, R
    il=repere(1,0);//attention pas la bonne ligne
    jl=repere(1,1);
    kl=repere(1,2);

    //détermination des points
    vector<double> pta1,ptb1,ptc1,pta2,ptb2,ptc2;
    if ((coord[0][0]<coord[1][0]) && (coord[0][0]<coord[2][0])) {
      pta1=coord[0];
      pta2=coord[3];
      if (coord[1][1]<coord[2][1]) {
        ptb1=coord[1];
        ptb2=coord[4];
        ptc1=coord[2];
        ptc2=coord[5];
      } else {
        ptb1=coord[2];
        ptb2=coord[5];
        ptc1=coord[1];
        ptc2=coord[4];
      }
    } else if (coord[1][0]<coord[2][0]) {
      pta1=coord[1];
      pta2=coord[4];
      if (coord[0][1]<coord[2][1]) {
        ptb1=coord[0];
        ptb2=coord[3];
        ptc1=coord[2];
        ptc2=coord[5];
      } else {
        ptb1=coord[2];
        ptb2=coord[5];
        ptc1=coord[0];
        ptc2=coord[3];
      }
    } else {
      pta1=coord[2];
      pta2=coord[5];
      if (coord[0][1]<coord[1][1]) {
        ptb1=coord[0];
        ptb2=coord[3];
        ptc1=coord[1];
        ptc2=coord[4];
      } else {
        ptb1=coord[1];
        ptb2=coord[4];
        ptc1=coord[0];
        ptc2=coord[3];
      }
    }


    double surf, voltot;
    surf=surface_triangle(pta1,ptb1,ptc1);
    voltot=surf*vrdt;

    //compter le nombre de points en dehors de la surface
    int nbpt;
    nbpt=0;
    if (pta2[0]<0) {
      nbpt+=1;
    }
    if (ptb2[1]<0) {
      nbpt+=1;
    }
    if (ptc2[2]<0) {
      nbpt+=1;
    }

    //détermination du cas
    if (nbpt==0) {
      _C_solide[il][jl][kl]-=voltot/(_dx*_dy*_dz);
    } else if (nbpt==1) {
      if (ptb2[1]<0) {
        recul3D::rotationcoin(repere, coord);
        recul3D::rotationcoin(repere, coord);
        vector<double> point;
        point=pta1;
        pta1=ptb1;
        ptb1=ptc1;
        ptc1=point;
        point=pta2;
        pta2=ptb2;
        ptb2=ptc2;
        ptc2=point;
      } else if (ptc2[2]<0) {
        recul3D::rotationcoin(repere, coord);
        vector<double> point;
        point=pta1;
        pta1=ptc1;
        ptc1=ptb1;
        ptb1=point;
        point=pta2;
        pta2=ptc2;
        ptc2=ptb2;
        ptb2=point;
      }
      /* code *///cas 1 coin sortant
      vector<double> ptd, pte, ptf;
      ptd.resize(3);
      pte.resize(3);
      ptf.resize(3);
      //D intersection de A1A2 x=0
      double coeff_a;
      coeff_a=-pta2[0]/(pta1[0]-pta2[0]);
      ptd[0]=0;
      ptd[1]=pta2[1]+coeff_a*(pta1[1]-pta2[1]);
      ptd[2]=pta2[2]+coeff_a*(pta1[2]-pta2[2]);
      //E intersection A2B2 x=0
      coeff_a=-ptb2[0]/(pta2[0]-ptb2[0]);
      pte[0]=0;
      pte[1]=ptb2[1]+coeff_a*(pta2[1]-ptb2[1]);
      pte[2]=ptb2[2]+coeff_a*(pta2[2]-ptb2[2]);
      //F intersection A2C2 x=0
      coeff_a=-ptc2[0]/(pta2[0]-ptc2[0]);
      ptf[0]=0;
      ptf[1]=ptc2[1]+coeff_a*(pta2[1]-ptc2[1]);
      ptf[2]=ptc2[2]+coeff_a*(pta2[2]-ptc2[2]);

      //pyramide A2DEF dans cube D(ligne 0 dans repère)
      double vol_d;
      vol_d = volume_pyramide(pta2, ptd, pte, ptf);
      _C_solide[repere(0,0)][repere(0,1)][repere(0,2)]-=vol_d/(_dx*_dy*_dz);
      //reste dans cube E
      _C_solide[repere(1,0)][repere(1,1)][repere(1,2)]-=(voltot-vol_d)/(_dx*_dy*_dz);

    } else if (nbpt==2) {
      if (pta2[0]>0) {
        recul3D::rotationcoin(repere, coord);
        recul3D::rotationcoin(repere, coord);
        vector<double> point;
        point=pta1;
        pta1=ptb1;
        ptb1=ptc1;
        ptc1=point;
        point=pta2;
        pta2=ptb2;
        ptb2=ptc2;
        ptc2=point;
      } else if (ptb2[2]>0) {
        recul3D::rotationcoin(repere, coord);
        vector<double> point;
        point=pta1;
        pta1=ptc1;
        ptc1=ptb1;
        ptb1=point;
        point=pta2;
        pta2=ptc2;
        ptc2=ptb2;
        ptb2=point;
      }
      /* code *///cas 2 coins sortants

      vector<double> ptd1, pte1, ptf1;
      ptd1.resize(3);
      pte1.resize(3);
      ptf1.resize(3);
      double coeff_a;
      //D1 intersection de B1B2 y=0
      coeff_a=-ptb2[1]/(ptb1[1]-ptb2[1]);
      ptd1[0]=ptb2[0]+coeff_a*(ptb1[0]-ptb2[0]);
      ptd1[1]=0;
      ptd1[2]=ptb2[2]+coeff_a*(ptb1[2]-ptb2[2]);
      //E1 intersection A2B2 y=0
      coeff_a=-ptb2[1]/(pta2[1]-ptb2[1]);
      pte1[0]=ptb2[0]+coeff_a*(pta2[0]-ptb2[0]);
      pte1[1]=0;
      pte1[2]=ptb2[2]+coeff_a*(pta2[2]-ptb2[2]);
      //F1 intersection B2C2 x=0
      coeff_a=-ptc2[0]/(ptb2[0]-ptc2[0]);
      ptf1[0]=0;
      ptf1[1]=ptc2[1]+coeff_a*(ptb2[1]-ptc2[1]);
      ptf1[2]=ptc2[2]+coeff_a*(ptb2[2]-ptc2[2]);

      //pyramide B2D1E1F1 dans cube H(ligne 3 dans repère)
      double vol_h;
      vol_h = volume_pyramide(ptb2, ptd1, pte1, ptf1);
      _C_solide[repere(3,0)][repere(3,1)][repere(3,2)]-=vol_h/(_dx*_dy*_dz);

      vector<double> ptd, pte, ptf;
      ptd.resize(3);
      pte.resize(3);
      ptf.resize(3);
      //D intersection de A1A2 x=0
      coeff_a=-pta2[0]/(pta1[0]-pta2[0]);
      ptd[0]=0;
      ptd[1]=pta2[1]+coeff_a*(pta1[1]-pta2[1]);
      ptd[2]=pta2[2]+coeff_a*(pta1[2]-pta2[2]);
      //E intersection A2B2 x=0
      coeff_a=-ptb2[0]/(pta1[0]-ptb2[0]);
      pte[0]=0;
      pte[1]=ptb2[1]+coeff_a*(pta1[1]-ptb2[1]);
      pte[2]=ptb2[2]+coeff_a*(pta1[2]-ptb2[2]);
      //F intersection A2C2 x=0
      coeff_a=-ptc2[0]/(pta1[0]-ptc2[0]);
      ptf[0]=0;
      ptf[1]=ptc2[1]+coeff_a*(pta1[1]-ptc2[1]);
      ptf[2]=ptc2[2]+coeff_a*(pta1[2]-ptc2[2]);

      //pyramide A2DEF dans cube D(ligne 0 dans repère)
      double vol_d;
      vol_d = volume_pyramide(pta2, ptd, pte, ptf);
      _C_solide[repere(0,0)][repere(0,1)][repere(0,2)]-=vol_d/(_dx*_dy*_dz);

      //reste dans cube E
      _C_solide[repere(1,0)][repere(1,1)][repere(1,2)]-=(voltot-vol_d-vol_h)/(_dx*_dy*_dz);

    } else {
      /* code *///cas 3 coins sortants -> séparer les cas
      //compter le nb d'arrêtes sortantes
      int nbar;
      nbar=0;
      //tests arrêtes sortantes

      if (nbar==0) {
        /* code *///cas O arrête sortante

      } else if (nbar==1) {
        /* code *///cas 1 arrête sortante
        //attention rotation
      } else if (nbar==2) {
        /* code *///cas 2 arrête sortante
        //attention rotation
      } else {
        /* code *///cas 3 arrête sortante -> séparer les cas
      }
    }

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

  void recul3D::rotationz(MatrixXd& repere_prec, vector<vector<double>>& coord)
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

    repere_prec=repere_suiv;

    int n;
    n=coord.size();
    for (int i = 0; i < n; i++) {
      vector<double> pta, ptb;
      pta=coord[i];
      ptb.resize(3);
      ptb[0]=_dx/2+_dy/2-pta[1];
      ptb[1]=_dy/2-_dx/2+pta[0];
      ptb[2]=pta[2];
      coord[i]=ptb;
    }

  }

  void recul3D::reductionrepere(MatrixXd& repere_prec)
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

    repere_prec=repere_suiv;
  }

  void recul3D::rotationcoin(MatrixXd& repere_prec, vector<vector<double>>& coord)
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

    repere_prec=repere_suiv;

    int n;
    n=coord.size();
    for (int i = 0; i < n; i++) {
      vector<double> pta, ptb;
      pta=coord[i];
      ptb.resize(3);
      ptb[0]=pta[2];
      ptb[1]=pta[0];
      ptb[2]=pta[1];
      coord[i]=ptb;
    }
  }
