#ifndef OB_DISTGEOM_H
#define OB_DISTGEOM_H

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

#include <Eigen/Core>

// rdkit/Code/GraphMol/DistGeomHelpers/Embedder.cpp
// https://github.com/rdkit/rdkit/blob/8524652a86f8c050239bde3bf61f4dff642e4c7a/Code/GraphMol/DistGeomHelpers/Embedder.cpp#L924
class DistanceGeometry {
  public:
    void Setup(OpenBabel::OBMol mol);
    void GetGeometry(OpenBabel::OBMol);
  private:
    OpenBabel::OBMol _mol;
    Eigen::MatrixXf bounds;
    void InitBoundsMat(double defaultMin = 0.0, double defaultMax = 1000.0);
    void SetUpperBounds(int i, int j, double val);
    void SetLowerBounds(int i, int j, double val);
    void SetTopolBounds();
    void Set12Bounds();
    void Set13Bounds();
    void Set14Bounds();
    void Set15Bounds();
    void SetLowerBoundsVDW();
};

#endif // OB_DISTGEOM_H
