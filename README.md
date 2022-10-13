# O Riche August 4th 2022

This is a stable fork built from Duo's original NEMO4 integration branch:
commit: 39f4cb337542ce9fbbf364b32c9179d45d7d341f
date:   May 24th 2022

Changes will be marked with git tags. The current commit will be tagged as
cannemo4_or_dvlp_v0.0

Summary: 

- for now there are 4 BGC tracers (DIC/Alkalinity/O2/NO3) activated and
  DIC and O2 are being affected by surface solubity and gas exchange.
    
- for now some of PISCES code, e.g. sms_pisces.F90/p4zche.F90/p4zflx.F90
  has been modified/moved.
        
- the code has been moved to TOP with the original idea to have CanOE/CMOC
  shared code running there but in insight and previously discussed it
  will interfer with run using CFC or C14 components of NEMO4/PISCES.
  This is likely to be migrated to a BGCM component that contains CanOE
  and CMOC and will host shared code like the one above in the BGCM level
  and called within CanOE and CMOC code.
      
- code moved resides in new TOP files, i.e. sms_top.F90, trcche.F90, and
  trcflx.F90.
      
- PISCES array indices have been migrated to TOP and are shared across
  the code for now. Index values are assigned like previously in PISCES
  in trcini.F90/trc_init.
    
- diagnostics code native from NEMO4/PISCES has been activated and use
  to test the next feature
    
- new code added to be able to read external sources. The code resides
  within TOP as trcsrc.F90; there 3 subroutines, i.e. trc_src_init, 
  trc_src3d, and trc_src2d, the former reads namelist_top_cfg to get 
  info on the external source files and the 2 others load the variables
  in arrays src3d_dta and src3d_dta respectively. Index system similar
  to trb/trn/tra is used and based on the order of the sources in the
  namelist file.

The rest of the text below refers to many old features of the CanNEMO standalone 
built from NEMO 3.4

=========================================================================

## CanNEMO

The development of NEMO at CCCma, including ocean physics, biogeochemistry and sea-ice as well as associated scripts and code for running the model. 
CanNEMO is a component of [CanESM](https://gitlab.com/cccma/canesm).

## Layout

There are four directories at the top level of the repo:

 - bin : scripts for compiling code, and setting up runs on the ECCC HPCs.

 - lib : the `xnemo` module used for building CCCma jobstrings.

 - rtd : Run time diagnostic code.

 - nemo: The "actual" NEMO source code. Under nemo/CONFIG there are
         several directories beginning with "CCC_", and these house
         the official CCCma configurations of the model.


## License

CanESM is distributed under the [Open Government License - Canada version 2.0](https://open.canada.ca/en/open-government-licence-canada).
The NEMO code is government by the CeCILL FREE SOFTWARE LICENSE AGREEMENT, which permits distribution of modified code.

## Support disclaimer

This code is made available on an as-is basis. lt has been tested only on the computing facilities
within Environment and Climate Change Canada (ECCC). There is no guarantee that it will run on
other platforms or if it does, that it will run correctly. No support of any kind will be made available
to help users to run the model on their own system. README documents linked below describe the development process.
