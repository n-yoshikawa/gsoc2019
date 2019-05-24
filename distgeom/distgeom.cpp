#include "distgeom.h"

void DistanceGeometry::Setup(OpenBabel::OBMol mol) {
  _mol = mol;
  InitBoundsMat();
  SetTopolBounds();
}

void DistanceGeometry::GetGeometry(OpenBabel::OBMol mol) {
}

// rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp
void DistanceGeometry::InitBoundsMat(double defaultMin, double defaultMax) {
  bounds = Eigen::MatrixXf::Zero(_mol.NumAtoms(), _mol.NumAtoms());
  int nrow = bounds.rows();
  for (int i = 1; i < nrow; ++i) {
    for (int j = 0; j < i; ++j) {
      SetUpperBounds(i, j, defaultMax);
      SetLowerBounds(i, j, defaultMin);
    }
  }
}

// rdkit/Code/DistGeom/BoundsMatrix.h
void DistanceGeometry::SetUpperBounds(int i, int j, double val) {
  if (i < j) {
    bounds(i, j) = val;
  } else {
    bounds(j, i) = val;
  }
}

void DistanceGeometry::SetLowerBounds(int i, int j, double val) {
  if (i < j) {
    bounds(j, i) = val;
  } else {
    bounds(i, j) = val;
  }
}

// rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp
void DistanceGeometry::SetTopolBounds(void) {
  unsigned int nb = _mol.NumBonds();
  unsigned int na = _mol.NumAtoms();
  Set12Bounds();
  Set13Bounds();
  Set14Bounds();
  Set15Bounds();
  SetLowerBoundsVDW();
}

void DistanceGeometry::Set12Bounds() { }
void DistanceGeometry::Set13Bounds() { }
void DistanceGeometry::Set14Bounds() { }
void DistanceGeometry::Set15Bounds() { }
void DistanceGeometry::SetLowerBoundsVDW() {}


