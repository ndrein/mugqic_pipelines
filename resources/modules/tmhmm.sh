#!/bin/bash

#
# TMHMM
#

## NOTE (FL): 
# - Tools from cbs.dtu.dk under some sort of danish license and need to be downloaded manually.
# - Assumes gnuplot, gs and are avail on the system
# - Empty output problem: http://sourceforge.net/p/trinotate/mailman/message/31228028/
#
# You need an executable of the program decodeanhmm that runs under
# Unix. The program may already be in bin/decodeanhmm.
# 
# The scripts require perl 5.x
# 
# For plotting gnuplot is needed (making postscript plots).
# 
# When generating html output the postscript plots are converted to
# gif, and for this you need the programs ghostscript (gs) and ppmtogif.
# 
# After unpacking the directory you should
# 
# 1. Insert the correct path for perl 5.x in the first line of the scripts
#    bin/tmhmm and bin/tmhmmformat.pl (if not /usr/local/bin/perl).
# 2. Make sure you have an executable version of decodeanhmm in the bin
#    directory.
# 3. Include the directory containing tmhmm in your path.
# 4. Read the TMHMM2.0.guide.html.
# 5. Run the program.

SOFTWARE=tmhmm
VERSION=2.0c

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
ARCHIVE=$SOFTWARE-$VERSION.Linux.tar.gz
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget http://www.dropbox.com/s/2lgbq3ci3zvvm2d/$ARCHIVE -O $ARCHIVE
fi
tar zxvf $ARCHIVE

SOFTWARE_DIR=$SOFTWARE-$VERSION
cd $SOFTWARE_DIR

# ppm2gif, could not find this on the web. Link to ppm2tiff may work instead
ln -s `which ppm2tiff` bin/ppm2gif

# Update Perl script shebang
sed -i "s,#!/usr/.*/perl.*,#!/usr/bin/env perl," bin/tmhmm bin/tmhmmformat.pl

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX .
mv -i $SOFTWARE_DIR $INSTALL_DIR
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
then
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive/
fi

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by removing '_INSTALL_HOME' in $INSTALL_HOME and lowercasing the result
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME/_INSTALL_HOME/} | tr '[:upper:]' '[:lower:]'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD