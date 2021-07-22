#-------------------------------------------------------
#
# "Makefile" du programme de tabulation de l'opacité
#
#-------------------------------------------------------
# D.C., mai 2002, novembre 2009.
#-------------------------------------------------------
all: fmaker4.1
#-------------------------------------------------------
# Compilation avec "g77" :
#FC = g77
#OPT     = -O6 -fno-automatic -fugly-assumed -malign-double
#
#-------------------------------------------------------
# compilation avec "gfortran" :
#
FC = gfortran
#
fmaker4.1: fmaker4.1.f
	$(FC) $(OPT) -o fmaker4.1 fmaker4.1.f
#
clear1:
	rm -f *~ ; rm -f *.c
