#-------------------------------------------------------
#
# "Makefile" du programme de tabulation des réactions
#            nucléaires pour CESAM 3.2
#
#-------------------------------------------------------
# D.C., avril 2001.
#----------------------------------------
# Compilateur, options et bibliothèques :
#
# Compilation avec "g77" :
#FC = g77
#OPT     = -O6 -fno-automatic -fugly-assumed -malign-double
#
# Compilation avec "gfortran" (il n'y a pas besoin d'options) :
FC = gfortran
#
#-------------------------------------------------------
#
all: tab_bin tab_bin_NACRE
#
#-------------------------------------------------------
tab_bin: tab_bin.f
	$(FC) $(OPT) -o tab_bin tab_bin.f
#
#-------------------------------------------------------
tab_bin_NACRE: tab_bin_NACRE.f
	$(FC) $(OPT) -o tab_bin_NACRE tab_bin_NACRE.f
#
#
clear1:
	rm -f *~