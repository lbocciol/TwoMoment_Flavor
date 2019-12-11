# TwoMoment_Flavor
Decoherence following Richers (2019)

With this code I reproduced FIG 3, 18, 1st column of 23 from Richers (2019).

- Decoherence.f90 is the main driver. I'm not sure if it works with nNodes > 1, it might not.
- InitializationModule.f90 initializes the profile from CHIMERA, or if you prefer
you can assing constant values and do a SingleZone calculation. 
The radiation is initialized to a FD, while the non diagonal elements are initialized
according to Richers (2019) (see IsotropicSQA)
- FlavorOpacitiesModule.f90: here is where the opacities from weaklib are converted to 
the form detailed in Richers (2019), and a brand new kernel (nu-nu scattering) is computed
according to Blaschke & Cirigliano (2016). The structure of the subroutines is the same as in 
the version of Nulib contained in IsotropicSQA.
- IntegrationModule.f90 does the implicit step. The timestep is adjusted by requiring that 
the difference between t + dt and t + 0.8dt is less than some tolerance.
There's the possibility to perform also an explicit RK6 integration, but I think it only
works for nu-nu interactions only (i.e. without EmAb)
- ReadProfileModuel.f90 reads the CHIMERA profile
- RadiationFieldsModule.f90 is here as an example of how the indices are defined. 
You need to modify thonado/Modules/Fields/RadiationFieldsModule.f90 in order for this 
code to work. Carefule because the order in which the indicied are defined
DOES matter when you do Packf and Flattenf in FlavorOpacitiesModule.
