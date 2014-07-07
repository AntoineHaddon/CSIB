Default PISCES configuration
============================

_O.Riche July 07/09/2014_ This configuration is based on the NPB run, a slightly upgraded verions of the so-called NKZ run. There are a few differences with the original PISCES configuration. I added an ideal age tracer (a few lines mainly in the p4zbio.F90). The relaxation of nitrate, silicic acid and alkalinity has been turned off. Finally, this version uses the old hack of two top namelists to get the first year initialized with observation and the next years initialized with the restart file.
