O Riche, 2014 December 13th

Trimming of most of PISCES code from NEMO-CMOC code
===================================================

Summary of the changes
----------------------

### CCC NEMO v0.1 code with CMOC configuration

The present code is a trimming of the code of run NQD (my user-specific 3-letter run name) which is known as __ccc-nemo-v0.1__ with **cmoc configuration** under the current nomenclature and version control system.

For the sake of clarity, I would like to amplify that the cmoc configuration is the implementation of the cmoc carbon model *[Zahariev et al 2008]* into the NEMO/PISCES code (stable version 3.4.1).

It uses the PISCES forcing. An iron mask (monthly, as far as I recall) diagnoses the iron limitation on phytoplankton growth.

The carbon chemistry is calculated by the native PISCES code (see *p4zchem.F90* and *p4zlys.F90*).

The balance of alkalinity, DIC, and N fixation/denitrification is enforced in the water column at each grid point (see *p4zrem.F90*). 

Some modules have been already removed, e.g. mesoplankton module, *p4zmeso.F90*.

This code expects a new namelist containing the cmoc parameters, i.e. *namelist_cmoc*.

All changes have been marked with the prefix "<OR CMOC ..." and a date and some details about the changes.


### The New Code

The present code is the continuation of the effort of accelerating the execution of the BGC model under orca 2 configuration and a fortiori under higher-resolution configurations.

All lines of code relevant to PISCES useless to CMOC and extranumerary tracers (CaCo3, PO4, Si, DOC, PHY2, DSi, Fer, BFe, SFe, DFe, GSi, NFe, DCHL, NH4) have been commented out as before in the previous instance of the NEMO-CMOC code.

The removal of tracers requires to modify not only the modules of the code but all the modules required for intialization of these tracers and associated housekeeping necessary to the model, e.g. par_pisces.F90 that contains the index associated with each tracer.

This is also true in the case of parameters which are stored in namelists. Finally, this is also the case of arrays and their memory allocation.


The Changes
-----------------

### The Namelists

The following namelists and relevant sections have been modified:

* *namelist_top*  (&namtrc for the tracers and their indices; &namtrc_dta for the initialization by input data files; &namtrc_trd for the diagnostic of tracers via their indices)

* *namelist_pisces* (&nampisbio to remove ferat3 Fe/C ratio; &nampissed to remove dust, n deposition, and iron sedimentation, in effect most of this section has been trimmed): have to restore dust as required by Jim (done in a later version of the NEMO-CMOC code, but contains an argon tracer for the time being). Some sections like &nampiskrs and &nampisdia, kriest parameterization and diagnostics respectively, could be trimmed to make the namelists clearer to browse and modify.

* *orca2_lim_pisces_iodef.xml* (definition section, removed field Irondep; output files definition section, added 1m_grid_T fields "qsr" and "qsn" for diagnostic of PAR, they could be removed in a later version; in 1m_diad_T section, potientially remove sDIC, sAlk, BUPOC, BUCALC, LNut,LNFe, LNligh, PAR, Mumax, all of these are diagnostics for debugging purposes; the corresponding lines in the code should be removed as well)

### The PISCES code

The following modules have been removed: p4zint.F90, p4zmeso.F90, p4zlim.F90.

The following modules have been modified (trimming summarized earlier):

* *p4zbio.F90*: removed "use p4zlim" in the header
   subroutine "p4z_bio": removed preprocessor blocks (beginning and end of the subroutine) associated with key_kriest, and call to module p4z_lim
* *p4zlys.F90*: removed lines associated with calcite tracer, with the calculation of parameter *zdispot*, 3 lines updating the tracer alkalinity (*jptal*), calcite (*jpcal*), and DIC (*jpdic*).

* *p4zmicro.F90*: removed "use p4zlim", "use p4zsink", "use p4zint" from the header

* *p4zmort.F90*: removed "use p4zsink" from the header, removed p4z_diat, merged p4z_mort and p4z_nano into p4z_mort (this implies that any reference to p4z_nano and p4z_diat have been removed from other modules too)

* *p4zprod.F90*: removed "use p4zlim"

* *p4zopt.F90*: removed diatom chlorophyll attenuation of light and fixed bug reported on mailing list

* *p4zrem.F90*: removed "use p4zopt", "use p4zprod", "use p4zint" from the header, removed most of the PISCES parameters and associated arrays.
   This module contains the CMOC N2 fixation and denitrification.

* *p4zsed.F90*: removed "use p4zopt", "use p4zlim", "use p4zrem", "use p4zint".
   This module contains the dust, river, and iron sedimentation initialization and forcing.

* *p4zsink.F90*: trimmed the PISCES code as done in the other modules.
   This module contains the particles aggregation and the Kriest parameterization, the later has not been commented out because it depends on pre-processor key key_kriest

* *par_pisces.F90*: reduced the number of tracers and associated a variable name to the appropriate tracer index.

* *sms_pisces.F90*: removed pisces parameters and arrays.

* *trcini_pisces.F90*: removed "use p4zlim"

* *trcnam_pisces.F90*: removed pisces parameters

* *trcsms_pisces.F90*: removed "use p4zint", removed calls to trc_sms_pisces_mass_conserv, p4z_int, and trc_sms_pisces_dmp (relaxation of PO4 and Si tracers), removed trc_sms_pisces_mass_conserv subroutine

* *trcrst_pisces.F90*: removed refences to Silica parameters

* *trcwri_pisces.F90*: removed references to PO4 and NH4, removed associated conversion block

