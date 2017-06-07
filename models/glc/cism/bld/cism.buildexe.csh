#! /bin/csh -f 

cd $OBJROOT/glc/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.cism
$CODEROOT/glc/cism/drivers/cpl_share
$CODEROOT/glc/cism/drivers/cpl_$comp
$CODEROOT/glc/cism/source_glc
$CODEROOT/glc/cism/source_glimmer-cism
$CODEROOT/glc/cism/source_slap
$CODEROOT/glc/cism/mpi
EOF

gmake complib -j $GMAKE_J MODEL=cism COMPLIB=$LIBROOT/libglc.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

