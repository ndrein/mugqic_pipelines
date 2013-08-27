

###################
################### BWA
###################
VERSION="2.3.14"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/igvtools/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download
wget http://www.broadinstitute.org/igv/projects/downloads/igvtools_$VERSION.zip
unzip igvtools_$VERSION.zip
# Compile
mv IGVTools igvtools-$VERSION
cd igvtools-$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - IGVtools \"
}
module-whatis \"MUGQIC - IGVtools  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/igvtools/igvtools-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/igvtools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/igvtools/


