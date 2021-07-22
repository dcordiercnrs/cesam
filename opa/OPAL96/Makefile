#-------------------------------------------------------
#
# "Makefile" du programme de tabulation de l'opacité
#
#-------------------------------------------------------
# D.C., mai 2002, novembre 2009.
#
#-------------------------------------------------------
all: fmttobin_gfortran
# Compilation avec "g77" :
#FC = g77
#OPT     = -O6 -fno-automatic -fugly-assumed -malign-double
#
#-------------------------------------------------------
# Compilation avec "gfortran" :
FC = gfortran
#
fmttobin_gfortran: fmttobin_gfortran.f
	$(FC) $(OPT) -o fmttobin_gfortran fmttobin_gfortran.f
#
clear1:
	rm -f *~ ; rm -f *.c
