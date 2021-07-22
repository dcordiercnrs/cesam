#################################################################
#
#            Makefile des executables CESAM 
#
#################################################################
#
FC       = g77
OPT     = -O6 -fno-automatic -fugly-assumed -malign-double -o
LIB     = ../lib_CRES4.0.a -lpgplot -L/usr/X11R6/lib -lX11
#
#################################################################
#
all: exe_cres4.0
#################################################################
#
exe_cres4.0: exe_cres4.0.f
	$(FC) $(OPT) exe_cres4.0 exe_cres4.0.f $(LIB)
