#include "Dense"

// classe maillage
class Maillage
{
  protected:

  public:
    // constructeur par défaut
    Maillage();
    // destructeur
    ~Maillage();
};

// classe fille maillage 2D cartésien
class Cartesien2D : public Maillage
{
  protected:
    // nom du maillage
    std::string _name;
    // taille du maillage et raffinement
    int _nbX, _nbZ;
    double _deltaX, _deltaZ, _Lz, _Lx;
    // vecteurs coordonnées
    Eigen::VectorXd _coordX, _coordZ;
    // vecteur indices
    Eigen::VectorXd _indices;

    // C'est Thomas tu peux rajouter Nx et Nz dans les attributs et mettre une méthode GetNx / GetNz
    // Okep c'est fait dude

  public:
    // constructeur par défaut
    Cartesien2D();
    // constructeur 1
    Cartesien2D(std::string name, double deltaX, double deltaZ, double Lx, double Lz);
    // destructeur
    ~Cartesien2D();

    //fonctions d'accès aux variables de l'objet
    Eigen::VectorXd GetCoordX(){return _coordX;};
    Eigen::VectorXd GetCoordZ(){return _coordZ;};
    Eigen::VectorXd GetIndices(){return _indices;};
    int GetNx(){return _nbX;};
    int GetNz(){return _nbZ;};

};

// classe fille maillage 3D cartésien
class Cartesien3D : public Maillage
{
  protected:
    // nom du maillage
    std::string _name;
    // taille du maillage et raffinement
    int _nbX, _nbY, _nbZ;
    double _deltaX, _deltaY, _deltaZ, _Lz, _Lx, _Ly;
    // vecteurs coordonnées
    Eigen::VectorXd _coordX, _coordZ, _coordY;
    // vecteur indices
    Eigen::VectorXd _indices;

  public:
    // constructeur par défaut
    Cartesien3D();
    // constructeur 1
    Cartesien3D(std::string name, double deltaX, double deltaY, double deltaZ, double Lx, double Lz, double Ly);
    // destructeur
    ~Cartesien3D();

    //fonctions d'accès aux variables de l'objet
    Eigen::VectorXd GetCoordX(){return _coordX;};
    Eigen::VectorXd GetCoordY(){return _coordY;};
    Eigen::VectorXd GetCoordZ(){return _coordZ;};
    Eigen::VectorXd GetIndices(){return _indices;};
    int GetNx(){return _nbX;};
    int GetNy(){return _nbY;};
    int GetNz(){return _nbZ;};

};
