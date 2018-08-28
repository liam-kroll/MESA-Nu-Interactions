# ***********************************************************************
#
#  Copyright (C) 2010  Bill Paxton, edited 2018 Liam Kroll
#
#  MESA is free software; you can use it and/or modify
#  it under the combined terms and restrictions of the MESA MANIFESTO
#  and the GNU General Library Public License as published
#  by the Free Software Foundation; either version 2 of the License,
#  or (at your option) any later version.
#
#  You should have received a copy of the MESA MANIFESTO along with
#  this software; if not, it is available at the mesa website:
#  http://mesa.sourceforge.net/
#
#  MESA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU Library General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# ***********************************************************************
# This is the makefile for the mesa nuclear reaction networks library

MESA_DIR = ../..

#################################################################
#
# PREAMBLE

include $(MESA_DIR)/utils/makefile_header


#################################################################
#
# SOURCE FILES
# Changed below, added the net_neu_reac.f90 file

SRCS = \
   net_def.f90 \
   net_neu_reac.f90 \
   net_screen.f90 \
   net_derivs_support.f90 \
   net_derivs.f90 \
   net_approx21.f90 \
   net_approx21_plus_co56.f90 \
   net_approx21_plus_fe53_fe55_co56.f90 \
   net_burn_support.f90 \
   net_initialize.f90 \
   net_eval.f90 \
   net_burn.f90 \
   net_burn_const_density.f90 \
   net_burn_const_p.f90 \
   net_lib.f90

#################################################################
#
# TARGETS

LIB = libnet.a
LIB_OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
$(LIB) : $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

#################################################################
#
# COMPILATION RULES

# also in test makefile
APROX21_FXT_FLAG =
#APROX21_FXT_FLAG = -DAPROX21_FXT

#COMPILE = $(COMPILE_NO_OPT) $(FCfree) $(APROX21_FXT_FLAG)

#COMPILE = $(COMPILE_TO_TEST) $(FCfree) $(APROX21_FXT_FLAG)
COMPILE = $(COMPILE_TO_DEPLOY) $(FCfree) $(APROX21_FXT_FLAG)

COMPILE_CMD = $(COMPILE)

%.o : %.f
	$(COMPILE_CMD) $<
%.o : %.f90
	$(COMPILE_CMD) $<
%.mod : %.o
	@true
%.o : %.mod

clean:
	-@rm -f *.o *.mod *.modmic *genmod.f90 *.a .depend

#################################################################
#
# DEPENDENCIES

# Paths to search for source & include files

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

# Automatic dependency generation using makedepf90

nodeps : $(.DEFAULT_GOAL)

NODEPS = $(or $(filter nodeps,$(MAKECMDGOALS)),$(filter clean,$(MAKECMDGOALS)))

ifeq ($(NODEPS),)

  DEPEND_DIRS = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

  .depend :
	makedepf90 -m %m.mod -I$(DEPEND_DIRS) $(SRCS) > .depend
  -include .depend

endif
