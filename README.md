move
====

Quick and dirty program for moving atoms in PDB format

Compile using g++ or compiler of choice:

  g++ move.cpp -o move.exe

Modify file 'input' with indexes (starting from 1) to move
and coordinates to move by.

Then run:
  ./move.exe INITIAL_PDB FINAL_PDB

Only works with PDB files.
Only moves ATOM and HETATM type elements.
Not tested for extreme cases.
No periodic boundary stuff.
