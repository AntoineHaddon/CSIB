# CanNEMO

The development of NEMO at CCCma, including ocean physics, biogeochemistry and sea-ice as well as associated scripts and code for running the model. 

The `official` version of the code is on the ECCC gitlab server [https://eccc-gitlab.science.gc.ca/ncs001/CanNEMO], which includes issue tracking
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
