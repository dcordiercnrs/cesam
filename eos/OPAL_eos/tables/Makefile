#--------------------------------------------------------------------
#
#     Makefile for interpolation program of OPAL EOS Tables
#
#--------------------------------------------------------------------
# D. Cordier, April 2003, Warsaw. Novembre 2009, Rennes.
#
#--------------------------------------------------------------------
# Compilation avec "g77" :
#FC  = g77
#
#--------------------------------------------------------------------
# Compilation avec "gfortran" :
FC  = gfortran
#
all: Z_interpol_dc clear
#
Z_interpol_dc: Z_interpol_dc.f long.o
	$(FC) $(OPT) -o Z_interpol_dc Z_interpol_dc.f long.o
#
long.o: long.f
	$(FC) -c long.f
#--------------------------------------------------------------------
clear:
	rm -f *~