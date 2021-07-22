#################################################################
#
#        --- Travail sur les céphéides juin 2010 ---
#
# #################################################################
# D. Cordier, adaptation au CLUSTER de Besançon, 30 septembre 2011
# Avec compilateur GNU :
#FC       = g77
#FC      = ifort
#OPT     = -O6 -fno-automatic -fugly-assumed -malign-double -o
#OPT     = -extend-source -o
#LIB     = ../lib/lib_CRES4.1.a -lpgplot -L/usr/X11R6/lib -lX11
#LIB     = ../lib/lib_CRES4.0_ifort.a -lpgplot -lX11
#
# Avec compilateur INTEL (ne marche pas pour le moment, 8 mai 2002)
#CC  = ifc
#OPT = -O3 -tpp6 -cm -w -WB -extend_source -c
#
# Compilation sur le cluster avec gfortran :
FC = gfortran
OPT = -w -fno-automatic -ffixed-line-length-none -o 
LIB = ../lib/lib_CRES4.0_gfortran.a -lpgplot -lX11
LIBNODESS = ../lib/lib_CRES4.0_gfortran_nodess.a 
#
#################################################################
#
all: exe_accre_2011_gfortran 
#exe_accre_2011_gfortran_nodess
#
#################################################################
exe_accre_2010: exe_opal_nacre_gain.f
	$(FC) $(OPT) exe_accre_2010 exe_opal_nacre_gain.f \
	$(LIB)
#
exe_accre_2011_gfortran: exe_opal_nacre_gain.f
	$(FC) $(OPT) exe_accre_2011_gfortran_cluster exe_opal_nacre_gain.f \
        $(LIB)
#
exe_accre_2011_gfortran_nodess: exe_opal_nacre_gain_NOdess.f
	$(FC) $(OPT) exe_accre_2011_gfortran_nodess exe_opal_nacre_gain_NOdess.f $(LIBNODESS)
#
