#! /bin/csh -f 

cd $OBJROOT/lnd/obj

if (-f $CASEBUILD/clmconf/Filepath) then
   cp $CASEBUILD/clmconf/Filepath ./tmp_filepath 
else
   echo "clm.buildexe.csh ERROR - missing $CASEBUILD/clmconf/Filepath"
   exit -1
endif
if (-f Filepath) then
  cmp -s tmp_filepath Filepath || mv -f tmp_filepath Filepath 
else
  mv -f tmp_filepath Filepath 
endif

#------------------------------------------------------------------------------
# giac gcam and glm build
#------------------------------------------------------------------------------

set iacobj = $OBJROOT/iac/obj
if (! -d $iacobj) mkdir -p $iacobj

if ("$CLM_IAC_MODE" == "giac") then
#------------------------------------------------------------------------------
# Build GCAM
cd $iacobj
set gcamtop = $CODEROOT/lnd/clm/src/iac/giac/gcam
cp -p -r $gcamtop/cvs .
cp -p $INCROOT/iac_fields_mod* cvs/objects/ccsmcpl/source/ 
cd cvs/objects/build/linux
mkdir ./objs
gmake gcamccsm
cp libgcamlib.a $LIBROOT/libiac.a
#cp libgcamlib.a $LIBROOT/libgcam.a
cp ../../ccsmcpl/source/gcam_comp_mod.mod $INCROOT/
#------------------------------------------------------------------------------
#Build GLM
cd $iacobj
set glmtop = $CODEROOT/lnd/clm/src/iac/giac/glm
cp -p -r $glmtop .
cd glm
gmake clean
cp -p $INCROOT/iac_fields_mod* .
gmake glmlib
ar ru $LIBROOT/libiac.a ccsm_glm_interface.o  glm_comp_mod.o  glm.future.iesm.o
cp glm_comp_mod.mod $INCROOT/
#------------------------------------------------------------------------------
endif

#------------------------------------------------------------------------------
# iac build
#------------------------------------------------------------------------------
cd $iacobj
set iacpath = $CODEROOT/lnd/clm/src/iac
set iacfiles = "unknown"
if ("$CLM_IAC_MODE" == "siac") then
  set iacfiles = "$iacpath/siac"
endif
if ("$CLM_IAC_MODE" == "diac") then
  set iacfiles = "$iacpath/diac"
endif
if ("$CLM_IAC_MODE" == "giac") then
  set iacfiles = "$iacpath/giac/coupling"
endif

if ("$iacfiles" == "unknown") then
  echo "CLM_IAC_MODE problem in clm.buildexe.csh"
  exit -1
endif

cat >! Filepath <<EOF
$CASEROOT/SourceMods/src.${CLM_IAC_MODE}
$iacfiles
EOF

$GMAKE complib -j $GMAKE_J MODEL=${CLM_IAC_MODE} COMPLIB=$LIBROOT/libiac.a -f $CASETOOLS/Makefile || exit 2
cp iac_comp_mod.mod $INCROOT/
#------------------------------------------------------------------------------

cd $OBJROOT/lnd/obj
set clmdefs = "`cat $CASEBUILD/clmconf/CESM_cppdefs`"
$GMAKE complib -j $GMAKE_J MODEL=clm COMPLIB=$LIBROOT/liblnd.a USER_CPPDEFS="$clmdefs" -f $CASETOOLS/Makefile || exit 2

wait



