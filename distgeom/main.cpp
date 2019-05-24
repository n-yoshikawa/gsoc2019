#include <iostream>
#include <string>
#include <sstream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>

#include "distgeom.h"


int main(int argc, char** argv){
  std::string smiles = "FCl(=O)=O";
  std::stringstream ss(smiles);
  OpenBabel::OBConversion conv(&ss, &std::cout);
  conv.SetInAndOutFormats("smi", "sdf");
  OpenBabel::OBMol mol;
  conv.Read(&mol);

  DistanceGeometry dg;
  dg.Setup(mol);
  dg.GetGeometry(mol);

  conv.Write(&mol);
}
