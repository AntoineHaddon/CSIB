## CanNEMO

The development of NEMO at CCCma, including ocean physics, biogeochemistry and sea-ice as well as associated scripts and code for running the model. 
CanNEMO is a component of [CanESM](https://gitlab.com/cccma/canesm).

## Layout

There are four directories at the top level of the repo:

 - {bin/,lib/} : scripts for compiling code, and setting up runs on the ECCC HPCs.

 - rtd/ : Run time diagnostic code.

 - diag/: CMIP6 diagnostics calculations.

 - nemo/: The "actual" NEMO source code. Under nemo/CONFIG there are
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
