# User’s Guide

The Integrated Earth System Model (iESM) was designed to use the CESM infrastructure to create cases and run the model. Therefore, the steps used to run iESM are identical to the steps used to run CESM1.2. The model generates CESM history and restart files much like those of standalone CESM. The GCAM options are invoked by changing various options within CESM. GCAM results are written to a separate file within the CESM output directories. Information on running iESM and using GCAM within iESM is included in this guide.

## Running iESM

Setting up and running iESM uses a set of four commands:
``` 	
./create_newcase
./cesm_setup
./CASE.build
./CASE.submit 
```

### Create a case
To create a case, you need to specify the name of the case, the machine, the component set, and the resolution. For example, the command below would create an emissions-forced RCP4.5 called “CASE” on Titan with a resolution of roughly 1 degree.

``` 
./create_newcase -case CASE -mach titan -compset BRCP45BPRP -res f09_g16 
```

iESM has been run using RCP4.5 and RCP8.5 component sets, including both concentration- and emissions-driven runs. 

### Configuring GCAM
To use GCAM within iESM, you will need to “turn it on”. The following command will do this:

``` 
./xmlchange -file env_run.xml -id CLM_IAC_MODE -val 'giac' 
```

Out of the box, GCAM is configured to run an RCP4.5. If you want to try a different GCAM pathway, you will need to adjust the GCAM configuration file (configuration.xml). The default is to use the `configuration_exp1_clmcarbon.xml` as the configuration. To change this, you will need to either overwrite that file, or use a `sed` command to change the configuration file used. For example, the following command will use `my_configuration.xml` as the GCAM configuration.

``` 
sed -i 's/configuration_exp1_clmcarbon.xml/my_configuration.xml/g' Buildconf/clm.buildnml.csh
``` 

For information on how to adjust GCAM’s configuration, please see the GCAM repository (www.github.com/JGCRI/gcam-core). 

### Enabling feedbacks between GCAM and CESM
GCAM and CESM can exchange information about land use, land cover, CO<sub>2</sub> emissions, and land productivity every five years. The following subsections explain how to enable these features. By default, they are turned off. (Note: the commands listed in the subsections may need small syntax changes depending on the shell used. These commands work for bash).

#### Land Use, Land Cover Exchange
To use GCAM’s land use, land cover change to drive CESM, the following change is needed: 

```
cat > user_nl_clm << EOF                                                                                                                         
!----------------------------------------------------------------------------------                                                              
! Users should add all user specific namelist changes below in the form of                                                                       
! namelist_var = new_namelist_value                                                                                                              
!                                                                                                                                                
! Include namelist variables for drv_flds_in ONLY if -megan and/or -drydep options                                                               
! are set in the CLM_NAMELIST_OPTS env variable.                                                                                                 
!                                                                                                                                                
! EXCEPTIONS:                                                                                                                                    
! Set co2_ppmv           with CCSM_CO2_PPMV                      option                                                                          
! Set dtime              with L_NCPL                             option                                                                          
! Set fatmlndfrc         with LND_DOMAIN_PATH/LND_DOMAIN_FILE    options                                                                         
! Set finidat            with RUN_REFCASE/RUN_REFDATE/RUN_REFTOD options for hybrid or branch cases                                              
! Set glc_grid           with GLC_GRID                           option                                                                          
! Set glc_smb            with GLC_SMB                            option                                                                          
! Set maxpatch_glcmec    with GLC_NEC                            option                                                                          
!----------------------------------------------------------------------------------                                                              
  fsurdat     =  '/pic/projects/climate/csmdata/iac/giac/glm2iac/surfdata_0.9x1.25_CESM12_simyear1850_c131015.nc'                                
  fpftdyn     = './surfdata_iESM_dyn.nc'                                                                                                         
  hist_fincl2 = 'TOTVEGC','FROOTC','LIVECROOTC','DEADCROOTC','TOTSOMC','TOTLITC','CWDC','NPP','AGNPP','BGNPP','HR'                               
  hist_dov2xy = .true., .false.                                                                                                                  
  hist_mfilt  = 1,1                                                                                                                              
  hist_nhtfrq = 0,0                                                                                                                              
EOF  
```           

You will need to update the path to the `fsurdat` file depending on the machine you use.

### Land Productivity Exchange
To use changes in the productivity of land from CESM in GCAM, the following command is needed: 
```
sed -i 's/clm_iac_carbon_scaling = .false./clm_iac_carbon_scaling = .true./g' Buildconf/clm.buildnml.csh
```

#### CO<sub>2</sub> Emissions exchange
To pass CO<sub>2</sub> emissions from GCAM to CESM, the following two changes are needed. 

``` 
sed -i 's/co2flux_coupling = .false./co2flux_coupling = .true./g' Buildconf/clm.buildnml.csh
```

```                                                                                                             cat > user_nl_cam <<EOF                                                                                                                          
! Users should add all user specific namelist changes below in the form of                                                                       
! namelist_var = new_namelist_value                                                                                                              
  co2flux_fuel_file = './co2flux_iESM_dyn.nc'                                                                                                    
EOF                                                                           
```

### Machines
We have successfully run iESM on Cori, Edison, Titan, Mira, and Constance. Machine files for the first four machines are included in the code. If you want to configure iESM for a new machine, we recommend using those files as examples. There are other machine files in the repository, but they will need some extra attention to make sure GCAM works properly. In particular, you’ll want to ensure the GCAM libraries, compilers, and paths are appropriately defined for your new machine.
