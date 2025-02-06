#!/bin/sh

cp ./docs/Doxygen/Doxyfile ./docs/Doxygen/Doxyfile.back
sed -i 's/m=m2cpp.pl/m="perl .\/m2cpp.pl"/' ./Doxyfile
doxygen ./docs/Doxygen/Doxyfile
mv ./docs/Doxygen/Doxyfile.back ./docs/Doxygen/Doxyfile
