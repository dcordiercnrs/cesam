#--------------------------------------------------------------------
#
#              Makefile for OPAL EOS calling routines
#
#--------------------------------------------------------------------
# D. Cordier, April 2003, Warsaw.
#
FC  = g77
OPT = -O6 -c -fno-automatic -fugly-assumed -malign-double
#
all: esac_dc.o etat_opal_dc.o clear
#
esac_dc.o: esac_dc.f
	$(FC) $(OPT) -c esac_dc.f
#
etat_opal_dc.o: etat_opal_dc.f
	$(FC) $(OPT) -c etat_opal_dc.f
#
#--------------------------------------------------------------------
clear:
	rm -f *~