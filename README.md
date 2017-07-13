# CanNEMO

The development of NEMO at CCCma, including ocean physics, biogeochemistry and sea-ice as well as associated scripts and code for running the model. 

The `official` version of the code is on the Science Network gitlab server [https://gitlab.science.gc.ca/ncs001/CanNEMO], which includes issue tracking
and a wiki. This is only available from within the ECCC network, and hence stategic mirrors also exist, but are not formal.

## Layout

There are four directories at the top level of the repo:

 - bin : scripts for compiling code, and setting up runs on the ECCC HPCs.

 - lib : the `xnemo` module used for building CCCma jobstrings.

 - rtd : Run time diagnostic code.

 - nemo: The "actual" NEMO source code. Under nemo/CONFIG there are
         several directories beginning with "CCC_", and these house
         the official CCCma configurations of the model.

## Branches, workflow and contributions

Development follows the gitflow-like workflow.

- The head of the `master` branch always reflects the latest frozen verion of the code. Master
  is updated approximately once every 18 months, as defined under the CCCma development cycle.

- All work is done on `feature` branches.  

- Features are merged together on `develop`, which represents a working version of the code.

More information can be found on the CCCma [twiki](http://wiki.cccma.ec.gc.ca/cgi-bin/twiki/view/Main/NemoVersionControl)


## Contact

For questions: neil.swart@canada.ca

-----------------------------------
# Development of the CMIP6 diagnostics

## Options to output the CMIP6 diagnostics or the default CCCma fields

For the CMIP6 diagnostics: 
    set nemo_ar6_diag = on in make_orca_job (default).

For the default CCCma fields:
    1) set nemo_ar6_diag = off in make_orca_job;
    2) set output_level="0" in iodef.xml file;
    3) remove ' "add_key="key_trdtra key_diaar5" ' from compile file.

## Supported configurations

    1. CCC_ORCA1_LIM
    2. CCC_ORCA1_LIM_CMOC
    3. CCC_ORCA1_LIM_CANOE
    4. CCC_ORCA1_OFF_CMOC
    5. CCC_ORCA1_OFF_CANOE
    6. CCC_ORCA025_LIM

    Note that the CMIP6 diagnostics are only available for configs 1-3.

## Contact

    For problems and questions: duo.yang@canada.ca

D. YANG, 11/JUL/2017
-----------------------------------


