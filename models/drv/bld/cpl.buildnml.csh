#! /bin/csh -f

set cwd = `pwd`
cd ${CASEROOT}
source ./Tools/ccsm_getenv || exit -1
cd $cwd

if !(-d $CASEBUILD/cplconf) mkdir $CASEBUILD/cplconf
rm $CASEBUILD/cplconf/* >& /dev/null

cd $CASEBUILD/cplconf || exit -1

if (-e $CASEROOT/user_nl_cpl) then
  $CCSMROOT/scripts/ccsm_utils/Tools/user_nlcreate -user_nl_file \
    $CASEROOT/user_nl_cpl -namelist_name cpl_inparm >! $CASEBUILD/cplconf/cesm_namelist 
endif

$CCSMROOT/models/drv/bld/build-namelist \
    -infile $CASEBUILD/cplconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $CCSMROOT/scripts \
    -grid $GRID -rof_grid $ROF_GRID -atm_grid $ATM_GRID -lnd_grid $LND_GRID -ocn_grid $OCN_GRID || exit -2

if (-d ${RUNDIR}) then
   cp $CASEBUILD/cplconf/drv_in       ${RUNDIR}/. || exit -3
   cp $CASEBUILD/cplconf/seq_maps.rc  ${RUNDIR}/. || exit -4
   cp $CASEBUILD/cplconf/*modelio*    ${RUNDIR}/. || exit -5
endif



