#!/bin/bash
                                                                                
INSTALL_PATH=$STREAM_INSTALL_DIR
                                                                                
echo INSTALL_PATH = $INSTALL_PATH
                                                                                
echo Making Directories
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/lib
mkdir -p $INSTALL_PATH/bin
mkdir -p $INSTALL_PATH/database/properties

echo Installing STREAM Library Files in $INSTALL_PATH/lib
mv gridMotion/gridMotion_m.so $INSTALL_PATH/lib
mv turbo/turbo_m.so $INSTALL_PATH/lib
mv tetra/extensions_tetra_m.so $INSTALL_PATH/lib
mv FSIEXPLICIT/FSIEXPLICIT_m.so $INSTALL_PATH/lib

echo Installing STREAM Executables in $INSTALL_PATH/bin
if [ -e streamUns ] ; then
  mv streamUns $INSTALL_PATH/bin/$STREAM_EXEC
fi
if [ -e piso/streamUns ] ; then
  mv piso/streamUns $INSTALL_PATH/bin/$STREAM_EXEC-piso
fi

# If the install directory is the compile directory, do not copy database files.
if [ $PWD != $INSTALL_PATH/src ]
then
  echo Installing STREAM Property Files in $INSTALL_PATH/database/properties
  cp -p ../database/properties/*.tran $INSTALL_PATH/database/properties
fi

chmod -R a+rX $INSTALL_PATH
