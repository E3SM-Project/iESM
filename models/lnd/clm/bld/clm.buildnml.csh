#! /bin/csh -fv
set echo

if !(-d $CASEBUILD/clmconf) mkdir -p $CASEBUILD/clmconf

#--------------------------------------------------------------------
# Invoke clm configure - output will go in CASEBUILD/clmconf
#--------------------------------------------------------------------

set config_opts=" "
if ($LND_GRID == "reg" && $GRID != "CLM_USRDAT" ) then
   set config_opts=" -sitespf_pt $GRID"
endif
if ("$CCSM_COMPSET" =~ P* || "$CCSM_COMPSET" =~ R* ) then
   set config_opts=" -sitespf_pt $LND_GRID"
endif

cd $CASEBUILD/clmconf  
$CODEROOT/lnd/clm/bld/configure  $config_opts -comp_intf $COMP_INTERFACE \
    $CLM_CONFIG_OPTS -usr_src $CASEROOT/SourceMods/src.clm || exit -1 

#--------------------------------------------------------------------
# Create clm.buildnml.csh
#--------------------------------------------------------------------

if ($RUN_TYPE == startup ) then
   if ($CLM_FORCE_COLDSTART == on) then
     set START_TYPE = "cold"
   else
     set START_TYPE = "default"
   endif
else
   if ($RUN_TYPE == hybrid ) then
     set START_TYPE = "startup"
   else
     set START_TYPE = $RUN_TYPE
   endif
endif

set RESOLUTION = $LND_GRID
set clmusr     = ""
if ($LND_GRID == reg ) then
   if ( $GRID == CLM_USRDAT ) then
      set RESOLUTION = $CLM_USRDAT_NAME
      set clmusr     = " -clm_usr_name $CLM_USRDAT_NAME"
   else
      set RESOLUTION = $GRID
   endif
endif

set default_lnd_in_filename = "lnd_in"

set inst_counter = 1
while ($inst_counter <= $NINST_LND)

if ($NINST_LND > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set lnd_in_filename = ${default_lnd_in_filename}${inst_string}

setenv INST_STRING $inst_string

cd $CASEBUILD/clmconf  

if (-e $CASEBUILD/clm.input_data_list) rm $CASEBUILD/clm.input_data_list

if (-e $CASEROOT/user_nl_clm${inst_string}) then
  $UTILROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_clm${inst_string} \
    -namelist_name clm_inparm >! $CASEBUILD/clmconf/cesm_namelist  || exit -2
endif

set glc_opts = ""
if ("$COMP_GLC" != "sglc" )then
   set glc_opts = "-glc_grid $GLC_GRID -glc_smb .$GLC_SMB. "
endif

set usecase = " "
if ($CLM_NML_USE_CASE != "UNSET") set usecase = "-use_case $CLM_NML_USE_CASE"

set clm_startfile = " "
if ( $RUN_TYPE == "hybrid" || $RUN_TYPE == "branch" ) then
   set clm_startfile = "-clm_startfile ${RUN_REFCASE}.clm2.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
endif

$CODEROOT/lnd/clm/bld/build-namelist -infile $CASEBUILD/clmconf/cesm_namelist \
    -csmdata $DIN_LOC_ROOT  \
    -inputdata $CASEBUILD/clm.input_data_list \
    -namelist "&clm_inparm $CLM_NAMELIST_OPTS /" $usecase $glc_opts \
    -res $RESOLUTION $clmusr -clm_start_type $START_TYPE $clm_startfile \
    -l_ncpl $LND_NCPL -lnd_frac "${LND_DOMAIN_PATH}/${LND_DOMAIN_FILE}" \
    -glc_nec $GLC_NEC -co2_ppmv $CCSM_CO2_PPMV -co2_type $CLM_CO2_TYPE \
    -config $CASEBUILD/clmconf/config_cache.xml $CLM_BLDNML_OPTS || exit -3
    
if (-d ${RUNDIR}) then
  cp $CASEBUILD/clmconf/lnd_in ${RUNDIR}/$lnd_in_filename || exit -2
  # Only copy drv_flds_in namelist file if one doesn't already exist
  if ( ! -f "${RUNDIR}/drv_flds_in" ) cp $CASEBUILD/clmconf/drv_flds_in ${RUNDIR}/. >& /dev/null
endif

@ inst_counter = $inst_counter + 1

end


#===========================================
# iac stuff
#===========================================


if ($LND_GRID =~ 360x720* ) then
  set map_landuse  = "lnd/clm2/mappingdata/maps/halfdeg/map_0.5x0.5_MODIS_to_0360x0720_nomask_aave_da_c120329.nc"
  set map_fmax     = "lnd/clm2/mappingdata/maps/halfdeg/map_0.5x0.5_USGS_to_0360x0720_nomask_aave_da_c120329.nc"
  set map_flanwat  = "lnd/clm2/mappingdata/maps/halfdeg/map_0.5x0.5_AVHRR_to_0360x0720_nomask_aave_da_c120329.nc"
  set map_fvocef   = "lnd/clm2/mappingdata/maps/halfdeg/map_0.5x0.5_AVHRR_to_0360x0720_nomask_aave_da_c120329.nc"
  set map_fsoitex  = "lnd/clm2/mappingdata/maps/halfdeg/map_5x5min_IGBP-GSDP_to_0360x0720_nomask_aave_da_c120329.nc"
  set map_firrig   = " "
  set map_fglctopo = "lnd/clm2/mappingdata/maps/halfdeg/map_10x10min_nomask_to_0360x0720_nomask_aave_da_c120329.nc"
  set map_flndtopo = "lnd/clm2/mappingdata/maps/halfdeg/map_10x10min_nomask_to_0360x0720_nomask_aave_da_c120329.nc"
  set clmC_bfn_dir = "iac/giac/iac2gcam/iESM_exp1.1"
  set clm2gcam_mapfile = "iac/giac/iac2gcam/CCSM05d_2_GCAM_Lut.nc"
#####  a new base_clmfile needs to be generated based on 2000-2004
  set iac_base_clmfile = "iac/giac/iac2gcam/clm_for_gcam_15yr_means/year_1990-2004_exp1.1/iESM_exp1.1.clm2.h1.1990-2004_mean_exp1.1.nc"
#####  set mksurfout_init = "surfdata.pftdyn_0.5x0.5_C_S2_simyr2005-2100_Expt1.1.nc"; a new one of these for 2005 has to be generated
  set mksurfout_init = "surfdata.pftdyn_360x720_iESM112_simyear1850_c04022014.nc"
endif

if ($LND_GRID =~ 0.9x1.25* ) then
  set map_landuse  = "lnd/clm2/mappingdata/maps/0.9x1.25/map_0.5x0.5_landuse_to_0.9x1.25_aave_da_110307.nc"
  set map_fmax     = "lnd/clm2/mappingdata/maps/0.9x1.25/map_0.5x0.5_landuse_to_0.9x1.25_aave_da_110307.nc"
  set map_flanwat  = "lnd/clm2/mappingdata/maps/0.9x1.25/map_0.5x0.5_lanwat_to_0.9x1.25_aave_da_110307.nc"
  set map_fvocef   = "lnd/clm2/mappingdata/maps/0.9x1.25/map_0.5x0.5_vocef_to_0.9x1.25_aave_da_110307.nc"
  set map_fsoitex  = "lnd/clm2/mappingdata/maps/0.9x1.25/map_5minx5min_soitex_to_0.9x1.25_aave_da_110307.nc"
  set map_firrig   = "lnd/clm2/mappingdata/maps/0.9x1.25/map_5minx5min_irrig_simyr2000_to_0.9x1.25_aave_da_110307.nc"
  set map_fglctopo = "lnd/clm2/mappingdata/maps/0.9x1.25/map_10minx10min_topo_to_0.9x1.25_aave_da_110307.nc"
  set map_flndtopo = "lnd/clm2/mappingdata/maps/0.9x1.25/map_10minx10min_topo_to_0.9x1.25_aave_da_110307.nc"
  set clmC_bfn_dir = "iac/giac/iac2gcam/iESM_exp1.2"
  set clm2gcam_mapfile = "iac/giac/iac2gcam/CCSM_2_GCAM_lut.nc"
  set iac_base_clmfile = "iac/giac/iac2gcam/b.e11.B20TRBPRP.f09_g16.iESM_exp12_ctrl.001.clm2.h1_2000_2004_mean.nc"
  # this is the file for starting in 1850
  #set mksurfout_init = "surfdata.pftdyn_192x288_iESM112_simyear1850_c04232014.nc"
  # this is the file for starting in 2005
  set mksurfout_init = "surfdata.pftdyn_192x288_iESM112_simyear2005_FC12Hist_c05292014.nc"
endif


if ("$CLM_IAC_MODE" == "giac") then

set exedir = $RUNDIR; cd $exedir

# copy in input files
set gcamtop = $CODEROOT/lnd/clm/src/iac/giac/gcam
cp $gcamtop/exe/configuration_exp1_clmcarbon.xml ./configuration.xml
cp $gcamtop/exe/log_conf.xml .
cp $gcamtop/exe/carbon_tax_4p5.xml .
# cp $gcamtop/cvs/objects/model_data/EERE/erbinput2.csv .
cp $gcamtop/cvs/objects/model_data/base/core_model_input.xml .
cp $gcamtop/cvs/objects/model_data/base/iesm_cdensadj.xml .
cp $gcamtop/cvs/objects/magicc/inputs/input_gases.emk .
cp $gcamtop/cvs/objects/model_data/base/solver_config.xml .
cp $gcamtop/cvs/objects/magicc/inputs/co2hist_c.in .
cp $gcamtop/cvs/objects/magicc/inputs/maguser_c.cfg .
cp $gcamtop/cvs/objects/magicc/inputs/maggas_c.cfg .
cp $gcamtop/cvs/objects/magicc/inputs/magice_c.cfg .
cp $gcamtop/cvs/objects/magicc/inputs/magmod_c.cfg .
cp $gcamtop/cvs/objects/magicc/inputs/magrun_c.cfg .
cp $gcamtop/cvs/objects/magicc/inputs/magxtra_c.cfg .
cp $gcamtop/cvs/objects/magicc/inputs/qhalos_c.in .
cp $gcamtop/cvs/objects/magicc/inputs/qextra_c.in .
cp $gcamtop/cvs/objects/magicc/inputs/Co2input_c.dat .
cp $gcamtop/cvs/objects/magicc/inputs/BCOCHist_c.csv .

set glmtop = $CODEROOT/lnd/clm/src/iac/giac/glm

#create glm namelist files 

cat >! glm.fut.conf <<EOF
[old control]
hist_option       = 2
trun              = tone
cat3              = inputs/other/
cat5              = inputs/hyde_3.0/half_deg_grids/
cat6              = inputs/future/RCP_MiniCam/
foutput_dir        = output/

[control]
top level glm dir = $DIN_LOC_ROOT/iac/giac/glm/
case name         = test_chdata
runtype           = initial   # initial or restart
start year        = 2005
stop year         = 2100
output_text         = 0         # write out text state and lu
output_netcdf       = 1         # write out netcdf state and lu
input_text_files    = 0         # read in txt wood harvest files
input_netcdf_files  = 1         # read in netcdf wood harvest files
use_urban           = 0       # 1[Y] or 0[N]
res option          = 2       # 1 (1 degree) or 2 (0.5 degree) only option 2 is valid
number of countries = 192
future rate option  = 0       # HIST(0) | GCAM(3) | AIM(4) | IMAGE(5) | MESSAGE(6)
future scenario     = 0       # HIST(0) | GCAM(3) | AIM(4) | IMAGE(5) | MESSAGE(6)
num_regions         = 14      # HIST=0,GCAM=14,AIM=24,IMAGE=24,MESSAGE=24
gridded_woodharvest = 2       # HIST=0,GCAM=0,AIM=1,IMAGE=0,MESSAGE=1,GCAM_AEZ=2
logging_option      = 1       #  0=wh=zero,1=standard wh data, 4="nodata"
zdis_option         = 1       # option for algorithm for spatial allocation of wood harvest
                              # only option 1 is supported in this version

#priority for clearing and wood harvest
smart_flow_option        = 1       # primary(1) or secondary(2)

#agricultural residence option
#minimum flows only or 
#shifting cultivation within the locations defined by our SC map
adjust_smart_flow_option = 5  # minimum(1) or shifting(5)

# choose whether clearing for agriculture is counted 
# towards meeting wood harvest demand Y(2) or not N(1)
converted_forest_land_option  = 1 

#historical land-use dataset = HYDE         # HYDE or No-Data
nodata_option = 5                           # 5=HYDE3 or 6=No-Data

maxz                         = 21       # Maximum z before we get tired, and spread remaining harvest over all forested cells with z >= this value
best_case                    = 1        # 1= best case, 0= other case
best_case_min_flows_t5       = 1
best_case_min_flows_t4       = 0
total_harvest_switch         = 1
secondary_harvest_switch     = 1
virgin_harvest_switch        = 1
force_harvest_switch         = 1
cpavg                        = 1
tb2bb                        = 1.0     # total biomass to bole biomass ratio, value of 2.0 assumes the WH numbers are bole biomass only 
                                       # and we need to cut twice as much biomass  */

phbio_filename               =  inputs/other/phbio.average.7states.txt
phbio_length                 =  50    # filelength of probability of harvest given biomass

output_updated_states        = 1
output_updated_states2       = 0
output_updated_states3       = 0
output_lu                    = 1
country_primeflow_print      = 1

#glm.future options
region_test_tmp              = 0
region_test_gcode            = 11
region_test_index            = 6
conterminous                  = 0



[output_files]
static_regions_file=static_regions_file.txt
static_vba_file=static_vba_file.txt

[output directory]
output dir = /lustre/data/jet/GLM/output

[hyde_datasets]

hyde_crop_path = $DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gcrop_1500-2005.nc
hyde_othr_path =$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gothr_1500-2005.nc
hyde_past_path =$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gpast_1500-2005.nc
hyde_watr_path =$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gwatr.1500-1501.nc
hyde_icew_path =$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gicew.1500-1501.nc

[hyde_datasets_nodata]

crop_nodata_path =$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/nodata/gcrop_1500-1510.nc
other_nodata_path =$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/nodata/gothr_1500-1510.nc
past_nodata_path  =$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/nodata/gpast_1500-1510.nc

[future_datasets]

future_crop_constructed_states = inputs/GCAM_Expt1_constructed_states_2005-2100/gcrop.2005-2100.nc
future_secd_constructed_states = inputs/GCAM_Expt1_constructed_states_2005-2100/gsecd.2005-2100.nc
future_othr_constructed_states = inputs/GCAM_Expt1_constructed_states_2005-2100/gothr.2005-2100.nc
future_past_constructed_states = inputs/GCAM_Expt1_constructed_states_2005-2100/gpast.2005-2100.nc
future_watr_constructed_states = inputs/hyde_3.0/half_deg_grids/gwatr.1500-1501.nc
future_icew_constructed_states = inputs/hyde_3.0/half_deg_grids/gicew.1500-1501.nc
updated_initial_state = $DIN_LOC_ROOT/iac/giac/glm/inputs/initial/1700-2005/initial_state_Exp1.2_1700-2005.nc

[woodharvest_datasets]
woodharvest_file =$DIN_LOC_ROOT/iac/giac/glm/inputs/other/wood_harvest/woodharvest_1500-2005.nc
woodharvest_nodata_file =$DIN_LOC_ROOT/iac/giac/glm/inputs/woodharvest_nodata_1500-2005.nc
rcp_woodharvest_file =$DIN_LOC_ROOT/iac/giac/glm/inputs/other/wood_harvest/RCP/woodharvest_minicam_rcp_2005_2100.nc
rcp_woodharvest_nodata_file =$DIN_LOC_ROOT/iac/giac/glm/inputs/other/wood_harvest/RCP/woodharvest_minicam_rcp_2005_2100.nc
rcp_woodharvest_aez_file =$DIN_LOC_ROOT/iac/giac/glm/inputs/other/wood_harvest/RCP/woodharvest_rcp45_aez_nd_2005_2100.nc
rcp_woodharvest_nodata_aez_file =$DIN_LOC_ROOT/iac/giac/glm/inputs/other/wood_harvest/RCP/woodharvest_rcp45_aez_nd_2005_2100.nc

[other_datasets]

wh_region_codes_file   = inputs/other/wood_harvest/codes_halfdeg_minicam.txt
wh_cont2region_codes_file = inputs/other/wood_harvest/continent_codes_minicam_test.txt
whcodes2glm_map  = inputs/other/wood_harvest/codes2glm_gcam_turkey_in_ee.txt
aez_region_grid_file = inputs/other/wood_harvest/AEZ_region_grid.txt
aez_region_zone_file = inputs/other/wood_harvest/AEZ_zone_grid.txt

cellinfo_file  = inputs/other/cellarea/cellarea_halfdeg.txt
ccodes_file    = inputs/other/ccodes/ccodes.txt.sort2wh
ccodes_map     = inputs/other/ccodes/ccodes_half_deg.txt
cnames_file    = inputs/other/ccodes/cnames.txt.sort2wh
regnames_file  = inputs/other/wood_harvest/names_minicam.txt
contcodes_file = inputs/other/ccodes/continent.codes.txt.sort2wh
shiftcult_map  = inputs/other/shift_cult/shiftcult_map_halfdeg.txt
regcodes_map     = inputs/other/wood_harvest/regcodes_halfdeg.txt
gcodes_cont_map  = inputs/other/ccodes/gcodes_continent_half_deg_DUMMY.asc
miami_biomass_file_vba = inputs/other/miami_biomass_v3/miami_halfdeg_conform.txt
miami_biomass_file_vnppa = inputs/other/miami_npp/miami.half_deg.in_conform

[debug options]
smart_flow_bug_print         = 1
state_print                  = 0
state_bug_print              = 1
flow_bug_print               = 1

EOF

#1 change start stop years in glm config file
sed -i 's/start year.*/start year =   2005/' glm.fut.conf
sed -i 's/stop year.*/stop year  =   2100/' glm.fut.conf
sed -i 's/runtype.*/runtype    =    initial/' glm.fut.conf

if (${CONTINUE_RUN} == TRUE) then

#2 change runtype in glm config file
   sed -i 's/runtype.*/runtype    =    restart/' glm.fut.conf

#3 modify configuration.xml with restart-period value from rpointer.file
   sed -i '/restart-period/d' configuration.xml
   set period = `sed -n '2p' rpointer.gcam`
   sed -i '/climateOutputInterval/i \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ <Value name="restart-period">'"$period"'</Value>' configuration.xml

#4 modify configure with restart filename from rpointer.gcam
   sed -i '/xmlInputFileName/d' configuration.xml
   set gcamfile = `sed -n '1p' rpointer.gcam`
   sed -i '/xmlOutputFileName/i \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ <Value name="xmlInputFileName">'"$gcamfile"'</Value>' configuration.xml

else
#2 change runtype in glm config file
   sed -i 's/runtype.*/runtype    =    initial/' glm.fut.conf

#3 modify configuration.xml with restart-period value from rpointer.file
   sed -i '/restart-period/d' configuration.xml

#4 modify configure with restart filename from rpointer.gcam
   sed -i '/xmlInputFileName/d' configuration.xml
   set gcamfile = "core_model_input.xml"
   sed -i '/xmlOutputFileName/i \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ <Value name="xmlInputFileName">'"$gcamfile"'</Value>' configuration.xml
endif

# for glm2iac
set giactop = $DIN_LOC_ROOT/iac/giac
set mksurfin = mksurf_landuse_iESM_720x360.nc
set mksurfout = surfdata_iESM
set co2fluxout = ./co2flux_iESM
cp -f $giactop/glm2iac/pftregridparams.txt .
cp -f $giactop/glm2iac/surfdata_360x720_potveg.nc .
cp -f $giactop/glm2iac/mksrf_landuse_ingrid.nc ./$mksurfin
# these are the files for starting in 1850
#cp -f $giactop/glm2iac/iESM_Ref_CropPast1850_c04022014.nc ./iESM_Ref_CropPast.nc
#cp -f $giactop/glm2iac/surfdata_360x720_mcrop1850_c04022014.nc ./surfdata_360x720_mcrop.nc
# these are the files for starting in 2005
cp -f $giactop/glm2iac/iESM_Ref_CropPast2005_FC12Hist_c05292014.nc ./iESM_Ref_CropPast.nc
cp -f $giactop/glm2iac/surfdata_360x720_mcrop2005_FC12Hist_c05292014.nc ./surfdata_360x720_mcrop.nc
# these are needed as the reference for model years <2000
cp -f $giactop/glm2iac/surfdata_360x720_mcrop2000_c03062014.nc ./surfdata_360x720_mcrop2000.nc
cp -f $giactop/glm2iac/iESM_Ref_CropPast2000_c03282014.nc ./iESM_Ref_CropPast2000.nc

if (${CONTINUE_RUN} == FALSE) then
    cp -f $giactop/atm/cam/ggas/co2flux_fossil_1980-2005-monthly_0.9x1.25_c20100204.nc ${co2fluxout}_dyn.nc
    cp -f $giactop/glm2iac/${mksurfout_init} ${mksurfout}_dyn.nc
    cp -f ./iESM_Ref_CropPast.nc ./iESM_Dyn_CropPast.nc
    cp -f ./surfdata_360x720_mcrop.nc ./surfdata_360x720_mcrop_dyn.nc
endif
chmod 666 ./iESM_Dyn_CropPast.nc  ./surfdata_360x720_mcrop_dyn.nc ./iESM_Ref_CropPast.nc ./surfdata_360x720_mcrop.nc surfdata_360x720_potveg.nc  ${mksurfout}_dyn.nc ./iESM_Ref_CropPast2000.nc ./surfdata_360x720_mcrop2000.nc ${co2fluxout}_dyn.nc

# The following was the old base year 2000 reference file we are now switching to 1850 reference
# cp -f $giactop/glm2iac/surfdata_360x720_mcrop2000.nc ./surfdata_360x720_mcrop.nc
# cp -f $giactop/glm2iac/surfdata_360x720_mcrop1850_c08202013.nc ./surfdata_360x720_mcrop.nc
# use a new 1850 start file to facilitate the dynamic pl code file use
# The following was the old base year 2000 reference file we are now switching to 1850 reference
# cp -f $giactop/glm2iac/iESM_Expt1_C_S2_CropPast_Ref.nc ./iESM_Ref_CropPast.nc
# cp -f $giactop/glm2iac/iESM_Ref_CropPast1850_c08132013.nc ./iESM_Ref_CropPast.nc
# Use a new 1850 base file that has a standard format and only one record so it can be copied to create the dynamic pl code file
# mksrf_fdynuse      = './mksrf_landuse.nc'

set initrun = '.true.'
if ($CONTINUE_RUN == 'TRUE') set initrun = '.false.'

cat >! iac_in <<EOF
&iacnml
 clm_iac_carbon_scaling = .false.
 co2flux_coupling = .false.
 fast_oneway_iac_coupling = .false.
 npp_hr_on = .true.
 sneakermode = .false.
 initial_run = $initrun
 clmC_bfn_dir = '$DIN_LOC_ROOT/${clmC_bfn_dir}'
 clm2gcam_mapfile = '$DIN_LOC_ROOT/${clm2gcam_mapfile}'
 iac_base_clmfile = '$DIN_LOC_ROOT/${iac_base_clmfile}'
 gcam2glm_basecrop = '$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gcrop_1500-2005.nc'
 gcam2glm_basepast =  '$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gpast_1500-2005.nc'
 gcam2glm_baseothr =  '$DIN_LOC_ROOT/iac/giac/glm/inputs/hyde_3.0/half_deg_grids/gothr_1500-2005.nc'
 gcam2glm_basebiomass =  '$DIN_LOC_ROOT/iac/giac/glm/inputs/miami_biomass_conform_0.5x0.5_map.nc'
 gcam2glm_aezmap =  '$DIN_LOC_ROOT/iac/giac/glm/inputs/aez_18_reg_14_0.5x0.5_map.nc'
 gcam2emisfile_co2base2000='$DIN_LOC_ROOT/iac/giac/gcam2emisfile/gridcar_2000.nc'
 gcam2emisfile_co2shipbase2000='$DIN_LOC_ROOT/iac/giac/gcam2emisfile/IPCC_emissions_CO2_ships_2000_0.5x0.5_v1_20_04_2009.nc'
 gcam2emisfile_grid720x360='$DIN_LOC_ROOT/iac/giac/gcam2emisfile/gridcar_2000.720x360.gridarea.nc'
 gcam2emisfile_grid288x192='$DIN_LOC_ROOT/iac/giac/gcam2emisfile/co2flux_fossil_2005-monthly_0.9x1.25_gridarea_c20100204.nc'
 gcam2emisfile_lut720x360map='$DIN_LOC_ROOT/iac/giac/gcam2emisfile/LUT_cry2grid.nc'
 gcam2emisfile_downscaleinfo='$DIN_LOC_ROOT/iac/giac/gcam2emisfile/GCAM_Downscale_0.5deg.interp5.nc'
 gcam2emisfile_rcp45allsteps='$DIN_LOC_ROOT/iac/giac/gcam2emisfile/GHG_GCAM_CO2_all_steps_rcp4.5_Stab_pralit.interp5.new.nc'
 clm_nx = $LND_NX
 clm_ny = $LND_NY
/

&mksurfnml
 mksrf_fgrid    = '$DIN_LOC_ROOT/${map_landuse}'
 mksrf_fsoitex  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_soitex.10level.c010119.nc'
 mksrf_forganic = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_organic.10level.0.5deg.081112.nc'
 mksrf_flanwat  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_lanwat.050425.nc'
 mksrf_fmax     = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_fmax.070406.nc'
 mksrf_fglacier = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_glacier.060929.nc'
 mksrf_furban   = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_urban_3den_0.5x0.5_simyr2000.c090223_v1.nc'
 mksrf_fvegtyp  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_landuse_rc1850_c090630.nc'
 mksrf_fsoicol  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_soilcol_global_c090324.nc'
 mksrf_flai     = '$DIN_LOC_ROOT/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_lai_global_c090506.nc'
 mksrf_fglctopo = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_topo.10min.c080912.nc'
 mksrf_flndtopo = '$DIN_LOC_ROOT/lnd/clm2/rawdata/topodata_10min_USGS_071205.nc'
 mksrf_fvocef   = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_vocef.c060502.nc'
 mksrf_firrig   = ' '
 mksrf_fdynuse  = '$mksurfin'
 map_fpft       = '$DIN_LOC_ROOT/${map_landuse}'
 map_fglacier   = '$DIN_LOC_ROOT/${map_landuse}'
 map_fsoicol    = '$DIN_LOC_ROOT/${map_landuse}'
 map_furban     = '$DIN_LOC_ROOT/${map_landuse}'
 map_fmax       = '$DIN_LOC_ROOT/${map_landuse}'
 map_forganic   = '$DIN_LOC_ROOT/${map_landuse}'
 map_fglcmec_g2g= '$DIN_LOC_ROOT/${map_landuse}'
 map_flai       = '$DIN_LOC_ROOT/${map_landuse}'
 map_fharvest   = '$DIN_LOC_ROOT/${map_landuse}'
 map_flanwat    = '$DIN_LOC_ROOT/${map_flanwat}'
 map_fvocef     = '$DIN_LOC_ROOT/${map_fvocef}'
 map_fsoitex    = '$DIN_LOC_ROOT/${map_fsoitex}'
 map_firrig     = '$DIN_LOC_ROOT/${map_firrig}'
 map_fglctopo   = '$DIN_LOC_ROOT/${map_fglctopo}'
 map_flndtopo   = '$DIN_LOC_ROOT/${map_flndtopo}'
 map_fglcmec_t2g= '$DIN_LOC_ROOT/lnd/clm2/mappingdata/maps/topo/map_10minx10min_topo_to_0.5x0.5_landuse_aave_da_110228.nc'
 outnc_double   = .true.
 outnc_dims     = 2
 fsurdat        = '${mksurfout}.nc'
 fsurlog        = '${mksurfout}.log'
 fdyndat        = '${mksurfout}_dyn.nc '

/
EOF

# remove the temporary dbxml file to prevent hang
rm -f __db.00*

endif  

if ("$CLM_IAC_MODE" == "diac") then

set exedir = $RUNDIR; cd $exedir

# for glm2iac
set giactop = $DIN_LOC_ROOT/iac/giac
set mksurfin = mksurf_landuse_iESM_720x360.nc
set mksurfout = surfdata_iESM
cp -f $giactop/glm2iac/pftregridparams.txt .
cp -f $giactop/glm2iac/surfdata_360x720_potveg.nc .
cp -f $giactop/glm2iac/mksrf_landuse_ingrid.nc ./$mksurfin
# these are the files for starting in 1850
#cp -f $giactop/glm2iac/iESM_Ref_CropPast1850_c04022014.nc ./iESM_Ref_CropPast.nc
#cp -f $giactop/glm2iac/surfdata_360x720_mcrop1850_c04022014.nc ./surfdata_360x720_mcrop.nc
# these are the files for starting in 2005
cp -f $giactop/glm2iac/iESM_Ref_CropPast2005_FC12Hist_c05292014.nc ./iESM_Ref_CropPast.nc
cp -f $giactop/glm2iac/surfdata_360x720_mcrop2005_FC12Hist_c05292014.nc ./surfdata_360x720_mcrop.nc
# these are needed as the reference for model years <2000
cp -f $giactop/glm2iac/surfdata_360x720_mcrop2000_c03062014.nc ./surfdata_360x720_mcrop2000.nc
cp -f $giactop/glm2iac/iESM_Ref_CropPast2000_c03282014.nc ./iESM_Ref_CropPast2000.nc


if (${CONTINUE_RUN} == FALSE) then
    cp -f $giactop/glm2iac/${mksurfout_init} ${mksurfout}_dyn.nc
    cp -f ./iESM_Ref_CropPast.nc ./iESM_Dyn_CropPast.nc
    cp -f ./surfdata_360x720_mcrop.nc ./surfdata_360x720_mcrop_dyn.nc
endif
chmod 666 ./iESM_Dyn_CropPast.nc  ./surfdata_360x720_mcrop_dyn.nc ./iESM_Ref_CropPast.nc ./surfdata_360x720_mcrop.nc surfdata_360x720_potveg.nc ${mksurfout}_dyn.nc ./iESM_Ref_CropPast2000.nc ./surfdata_360x720_mcrop2000.nc 

# mksrf_fdynuse      = './mksrf_landuse.nc'
#cp -f $giactop/glm2iac/surfdata_360x720_mcrop1850_c08202013.nc ./surfdata_360x720_mcrop.nc
# use a new 1850 start file to facilitate the dynamic pl code file use
# now copy this initial file to the dynamic file
#cp -f $giactop/glm2iac/iESM_Ref_CropPast1850_c08132013.nc ./iESM_Ref_CropPast.nc
# Use a new 1850 base file that has a standard format and only one record so it can be copied to create the dynamic pl code file
# Now copy this initial file to the dynamic file

set initrun = '.true.'
if ($CONTINUE_RUN == 'TRUE') set initrun = '.false.'

cat >! iac_in <<EOF
&iacnml
 glm2iac_glmofile = '$DIN_LOC_ROOT/iac/giac/glm2iac/glmo_1700-2005_c08132013.nc'
/


&mksurfnml
 mksrf_fgrid    = '$DIN_LOC_ROOT/${map_landuse}'
 mksrf_fsoitex  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_soitex.10level.c010119.nc'
 mksrf_forganic = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_organic.10level.0.5deg.081112.nc'
 mksrf_flanwat  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_lanwat.050425.nc'
 mksrf_fmax     = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_fmax.070406.nc'
 mksrf_fglacier = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_glacier.060929.nc'
 mksrf_furban   = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_urban_3den_0.5x0.5_simyr2000.c090223_v1.nc'
 mksrf_fvegtyp  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_landuse_rc1850_c090630.nc'
 mksrf_fsoicol  = '$DIN_LOC_ROOT/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_soilcol_global_c090324.nc'
 mksrf_flai     = '$DIN_LOC_ROOT/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_lai_global_c090506.nc'
 mksrf_fglctopo = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_topo.10min.c080912.nc'
 mksrf_flndtopo = '$DIN_LOC_ROOT/lnd/clm2/rawdata/topodata_10min_USGS_071205.nc'
 mksrf_fvocef   = '$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_vocef.c060502.nc'
 mksrf_firrig   = ' '
 mksrf_fdynuse  = '$mksurfin'
 map_fpft       = '$DIN_LOC_ROOT/${map_landuse}'
 map_fglacier   = '$DIN_LOC_ROOT/${map_landuse}'
 map_fsoicol    = '$DIN_LOC_ROOT/${map_landuse}'
 map_furban     = '$DIN_LOC_ROOT/${map_landuse}'
 map_fmax       = '$DIN_LOC_ROOT/${map_landuse}'
 map_forganic   = '$DIN_LOC_ROOT/${map_landuse}'
 map_fglcmec_g2g= '$DIN_LOC_ROOT/${map_landuse}'
 map_flai       = '$DIN_LOC_ROOT/${map_landuse}'
 map_fharvest   = '$DIN_LOC_ROOT/${map_landuse}'
 map_flanwat    = '$DIN_LOC_ROOT/${map_flanwat}'
 map_fvocef     = '$DIN_LOC_ROOT/${map_fvocef}'
 map_fsoitex    = '$DIN_LOC_ROOT/${map_fsoitex}'
 map_firrig     = '$DIN_LOC_ROOT/${map_firrig}'
 map_fglctopo   = '$DIN_LOC_ROOT/${map_fglctopo}'
 map_flndtopo   = '$DIN_LOC_ROOT/${map_flndtopo}'
 map_fglcmec_t2g= '$DIN_LOC_ROOT/lnd/clm2/mappingdata/maps/topo/map_10minx10min_topo_to_0.5x0.5_landuse_aave_da_110228.nc'
 outnc_double   = .true.
 outnc_dims     = 2
 fsurdat        = '${mksurfout}.nc'
 fsurlog        = '${mksurfout}.log'
 fdyndat        = '${mksurfout}_dyn.nc '

/
EOF

# remove the temporary dbxml file to prevent hang
rm -f __db.00*

endif  

