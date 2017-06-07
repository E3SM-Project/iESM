#! /bin/csh -f

#-----------------------------------------------------------------------
## Run script for running CAM in single column mode on a pc platform.
## This runs a low resolution
## T5 Eulerian Spectral case in serial mode.
#-----------------------------------------------------------------------
##

## Do our best to get sufficient stack memory
limit stacksize unlimited

## netCDF stuff
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = $HOME/cam_trunk

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA    $HOME/inputdata

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: startup, continue, or branch.
## $stop_n is the number of timesteps to integrate (units depends on stop_option value)
set dyn          = "eul"
set ocn          = "dom"
set case         = scamrun
if ( $dyn != "eul" ) then
   echo "SCAM is supported using eularian dynamics only" && exit 1
endif
set runtype      = startup

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = ~/runs/
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld
set res          = "8x16"

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## build exec
if ( ! -x $blddir/cam )then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure -test -res $res -ocn $ocn -dyn $dyn -debug -scam || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    make -j4 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist for scam.  You must specify scamlat,scamlon,iopfile.
## The starting date/time should correspond to the start of the IOP file
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist  -namelist "&camexp \
scmlat=36.6 \
scmlon=262.5  \
start_ymd=19950718 \
start_tod=19800 \
stop_n=840 \
stop_option='nsteps' \
iopfile='$CSMDATA/atm/cam/scam/iop/arm0795v1.2.nc' \
mfilt=1400 \
nhtfrq=1/" -ignore_ic_date || echo "build-namelist failed" && exit 1

## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
mv $blddir/*in .
echo "running SCAM in $rundir"
$blddir/cam                 || echo "CAM run failed" && exit 1

exit 0
