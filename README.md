# MESA-Nu-Interactions
Adds reactions in MESA which simulate the charged current neutrino interactions on free nucleons in core-collapse supernovae

MESA-Nu-Interactions is a set of modified and added FORTRAN files which enable the user to simulate a time varying electron neutrino flux and the resultant charged current neutrino interactions on free nucleons in MESA simulated core-collapse supernovae. The rates for these reactions are based off the calculations by Qian & Woosley (1997). The code which performs the rates calculations is placed in run_star_extras so as to be easily accessible to the user.
To install this functionality the user should replace the original files in the mesa-r10398 directory with the modified and added files from MESA-Nu-Interactions/MESA_nu_int_files.
The following is a table of the files and the correct placement within $MESA_DIR:

File, Location within mesa-r10398

net_lib.f90, /net/public

net_approx21.f90, /net/private

net_approx21_plus_co56.f90, /net/private

net_approx21_plus_fe53 _fe55_co56.f90, /net/private

net_approx21_procs.inc, /net/private

net_eval.f90, /net/private

net_derivs.f90, /net/private

net_burn.f90, /net/private

net_burn_const_density.f90, /net/private

net_burn_const_p.f90, /net/private

makefile, /net/make

test_net_support.f, /net/test/src

sample_net.f, /net/test/src

test_net_do_one, /net/test/src

other_net_get.f90, /star/other

net.f90, /star/private

reactions.list, /data/rates_data

raw_rates.f90, /rates/private

rates_names.f90, /rates/private

rates_def.f90, /rates/public

test_output, /rates/test

run_star_extras.f90, In the work directory of a project where the functionality is desired
                                        
net_neu_reac.f90, /net/private

MESA can then be recompiled. If the compiler indicates MESA has been installed correctly then so has MESA-Nu-Interaction.
To use the added reactions in a core-collapse supernova simulation the user must include them in their reaction network. They have been named r_prot_wk_neut and r_neut_wk-minus_prot. However, by default the reactions have a null rate which does not effect the composition or energy generation. Only when the user has set the "use_other_net_get" inlist control parameter has been set to true, and the x_ctrls defining the properties of the neutrino flux have been set will the proper rates be calculated and used. This allows the user to continue to use MESA in the same way as pre-addition, but also turn on the neutrino interactions when desired entirely through the use of inlists. The charged nurrent neutrino interactions have also been added to the approx21 series of reaction networks. If the user does not wish to use them in these networks, the above still applies (i.e. if use_other_net_get = .false. then the rates used is effectively zero).

The MESA_Nu_Interactions/test_ccsn directory contains a MESA star work directory which can be used to test the implemented charged current neutrino interactions and familiarize the user.
