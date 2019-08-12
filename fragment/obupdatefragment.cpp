/**********************************************************************
  obupdatefragment - Update coordinate database of ring fragments

  Copyright (C) 2019 Geoffrey R. Hutchison

  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/builder.h>
#include <openbabel/parsmart.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <openbabel/elements.h>

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <queue>
#include <map>

using namespace std;
using namespace OpenBabel;

int main(int argc, char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  ios::sync_with_stdio(false);

  if (argc <= 3) {
    cerr << "Usage: " << argv[0] << " <existing fragment DB> <New DB File name> <File to add>" << endl;
    exit(EXIT_FAILURE);
  }


  // Read though the existing fragment DB
  ifstream ifs(argv[1]);
  if(!ifs) {
    cerr << "Falied to open " << argv[1] << endl;
    exit(EXIT_FAILURE);
  }

  string tempFileName = "temp.txt";
  ofstream ofs(tempFileName);
  if(!ofs) {
    cerr << "Falied to open " << tempFileName << endl;
    exit(EXIT_FAILURE);
  }

  map<string, int> SMILES2pos;          // store position in temporary file
  vector<pair<int, string>> NumAtoms;   // store number of atoms of fragment

  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;

  char buffer[BUFF_SIZE];
  vector<string> vs;
  int pos = ifs.tellg();
  while (ifs.getline(buffer, BUFF_SIZE)) {
    ofs << buffer << '\n';
    if (buffer[0] != '#') {
      tokenize(vs, buffer);

      if (vs.size() == 1) { // SMARTS pattern
        stringstream ss(vs[0]);
        conv.SetInStream(&ss);
        conv.Read(&mol);
        SMILES2pos[vs[0]] = pos;
        unsigned int n = mol.NumAtoms();
        NumAtoms.push_back({n, vs[0]});
      }
    }
    pos = ifs.tellg();
  }
  ifs.close();

  // Read new databases
  conv.SetOptions("O", conv.OUTOPTIONS);
  OBFormat *inFormat, *canFormat;
  OBAtom *atom;
  OBBond *bond;
  map<string, vector<OBMol> > fragment_list;
  map<string, int> fragment_count;
  vector<pair<int, string> > fragment_size;

  canFormat = conv.FindFormat("can"); // Canonical SMILES format
  conv.SetOutFormat(canFormat);

  for (int fileIdx = 3; fileIdx < argc; fileIdx++) {
    cerr << " Reading file " << argv[fileIdx] << endl;

    inFormat = conv.FormatFromExt(argv[fileIdx]);
    if(inFormat==NULL || !conv.SetInFormat(inFormat))
    {
      cerr << " Cannot read file format for " << argv[fileIdx] << endl;
      continue; // try next file
    }

    ifs.open(argv[fileIdx]);

    if (!ifs)
    {
      cerr << "Cannot read input file: " << argv[fileIdx] << endl;
      continue;
    }

    while(ifs.peek() != EOF && ifs.good())
    {
      conv.Read(&mol, &ifs);
      mol.DeleteHydrogens();

      unsigned int size = mol.NumAtoms();
      OBBitVec atomsToCopy(size+1);
      for (unsigned int i = 1; i <= size; ++i) {
        atom = mol.GetAtom(i);
        atomsToCopy.SetBitOn(atom->GetIdx());
      }

      size = mol.NumBonds();
      OBBitVec bondsToExclude(size);
      for (unsigned int i = 0; i < size; ++i) {
        bond = mol.GetBond(i);
        if (bond->IsRotor()) {
          bondsToExclude.SetBitOn(bond->GetIdx());
        }
      }

      OBMol mol_copy;
      mol.CopySubstructure(mol_copy, &atomsToCopy, &bondsToExclude);
      vector<OBMol> fragments = mol_copy.Separate(); // Copies each disconnected fragment as a separate
      for (unsigned int i = 0; i < fragments.size(); ++i) {
        if (fragments[i].NumAtoms() < 5) // too small to care
          continue;

        string smiles = conv.WriteString(&fragments[i], true);

        if (SMILES2pos.count(smiles) == 1) { // fragment alreadly exists
          continue;                          // TODO: need to decide how to decide multiple structure fragment
        }

        OBPairData *pd = dynamic_cast<OBPairData*>(fragments[i].GetData("SMILES Atom Order"));
        istringstream iss(pd->GetValue());
        vector<unsigned int> canonical_order;
        canonical_order.clear();
        copy(istream_iterator<unsigned int>(iss),
            istream_iterator<unsigned int>(),
            back_inserter<vector<unsigned int> >(canonical_order));

        unsigned int order;
        OBAtom *atom;

        fragments[i].Center(); // Translate to the center of all coordinates
        fragments[i].ToInertialFrame(); // Translate all conformers to the inertial frame-of-reference.

        SMILES2pos[smiles] = ofs.tellp();
        unsigned int n = fragments[i].NumAtoms();
        NumAtoms.push_back({n, smiles});

        // Write out an XYZ-style file with the CANSMI as the title
        ofs << smiles << '\n'; // endl causes a flush

        for (unsigned int index = 0; index < canonical_order.size(); ++index) {
          order = canonical_order[index];
          atom = fragments[i].GetAtom(order);

          snprintf(buffer, BUFF_SIZE, "%d %9.3f %9.3f %9.3f\n",
              atom->GetAtomicNum(),
              atom->x(), atom->y(), atom->z());
          ofs << buffer;
        } 
      }
    } // while reading molecules (in this file)
    ifs.close();
    ifs.clear();
  } // // while reading files

  ofs.close();

  // sort fragments by decreasing order of number of atoms
  sort(NumAtoms.rbegin(), NumAtoms.rend());

  ifs.open(tempFileName);
  if(!ifs) {
    cerr << "Failed to open " << tempFileName << endl;
    exit(EXIT_FAILURE);
  }

  ofs.open(argv[2]);
  if(!ofs) {
    cerr << "Failed to open " << argv[2] << endl;
    exit(EXIT_FAILURE);
  }

  // Write to the final output file
  for(auto x : NumAtoms) {
    string smiles = x.second;
    ifs.clear();
    ifs.seekg(SMILES2pos[smiles]);

    char buffer[BUFF_SIZE];
    vector<string> vs;
    bool isFirst = true;
    while (ifs.getline(buffer, BUFF_SIZE)) {
      tokenize(vs, buffer);
      if (vs.size() == 1) { // SMARTS pattern
        if (isFirst) isFirst = false;
        else break;
      }
      ofs << buffer << '\n';
    }
  }
}
