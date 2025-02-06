#!/bin/sh

cp ./Doxyfile ./Doxyfile.back
sed -i 's/m=m2cpp.pl/m="perl .\/m2cpp.pl"/' ./Doxyfile
doxygen ./Doxyfile
mv ./Doxyfile.back ./Doxyfile
