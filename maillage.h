#include "Dense"

class Maillage2DCarre
{
  protected:
    // nom du maillage
    std::string _name;
    // taille du maillage et raffinement
    int _deltaX, _deltaZ, _Lz, _Lx;
    // vecteurs coordonnées
    Eigen::VectorXd _coordX, _coordZ;

  public:
    // constructeur par défaut
    void Maillage2DCarre();
    // constructeur 1
    void Maillage2DCarre(std::string name, double deltaX, double deltaZ, double _Lx, double _Lz);
    // destructeur
    void ~Maillage2DCarre();

    //fonctions d'accès aux variables de l'objet
    VectorXd GetCoordX(){return _coordX};
    VectorXd GetCoordZ(){return _coordZ};
}
