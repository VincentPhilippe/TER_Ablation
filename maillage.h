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
class Cartesien : public Maillage
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
    Cartesien();
    // constructeur 1
    Cartesien(std::string name, double deltaX, double deltaZ, double _Lx, double _Lz);
    // destructeur
    ~Cartesien();

    //fonctions d'accès aux variables de l'objet
    Eigen::VectorXd GetCoordX(){return _coordX;};
    Eigen::VectorXd GetCoordZ(){return _coordZ;};
    Eigen::VectorXd GetIndices(){return _indices;};
    int GetNx(){return _nbX;};
    int GetNz(){return _nbZ;};

};
