####################################################
#
# Makefile de construction de l'executable de CRES
# ce Makefile doit etre lance par le Script "mex"
# ou utilis� seulement en tapant "make"
#
# D. Cordier, <daniel.cordier@ensc-rennes.fr,
#	Novembre 1999,
#       Mai 2000,
#	F�vrier 2001,
#       Octobre 2001,
#	Mai 2002,
#       Ao�t 2003, Varsovie.
#
####################################################
# Estimation du temps de compilation de l'ensemble :
#
#  - sur PIII 450 Mhz/384 Mo /g77 : env. 12 minutes
#
####################################################
#
# Avec compilateur GNU :
#FC       = g77
#FC = /Users/cordier/Downloads/usr/local/bin/g77 
#OPT      = -O6 -c -fno-automatic -fugly-assumed -malign-double
#
# Compilateur Intel "ifort", cluster d'UTINAM :
#FC = ifort
#OPT = -c -extend-source
#
# Compilateur gfortran install� sur le cluster d'UTINAM
FC  = gfortran
OPT = -c -w -O3 -fno-automatic -ffixed-line-length-none
# VERSION :
VER      = 4.0
#
HEAD = HEAD
# 
LIB_DIR = ../lib/
#
# Avec compilateur INTEL (ne marche pas pour le moment, 8 mai 2002)
#FC  = ifc
#OPT = -O3 -tpp6 -cm -w -WB -extend_source -c
#---------------------------------------------------
#
all: CESAM lib
#
lib:
	ar r lib_CRES4.0_gfortran.a *.o ; ranlib lib_CRES4.0_gfortran.a ; \
	mv lib_CRES4.0_gfortran.a $(LIB_DIR)
#
CESAM:  base_her_cte.o desfull_3.o etat_gong1_3.o int_roxburgh_m.o\
	no_des_3.o polyder.o setcnsc.o bilin.o dess_hr_DC.o etat_gong2_3.o\
	inteff.o noedif.o pp1dn.o setf4.o bilinc.o dess_ol_dig01.o\
	gain_masse.o etat_int_herm_3.o inteffc.o noeu_dis.o pp2d.o sgauss.o\
	bleqstc.o dgrad.o etat_opal.o inter_3.o noeud.o ppcno3a18.o shell.o\
	blk_eos.o difdiv.o evol_3.o intertot_3.o nrj_nuc.o ppcno3a_12_3.o\
	simq.o blstio.o diff_cte_3.o f4der.o intgauss.o nrj_nuc_Xcno.o\
	ppcno3a_12e.o slinf.o bspint.o diff_micro_3.o f4mhd.o k_5750.o nrj_nuc_cn.o\
	ppcno3a_9_3.o smatpr.o bval.o diff_z16.o f4n.o k_5750_3.o nuc_gong_3.o\
	quad.o sn1dn_dis.o bvald.o diffus_3.o ferdir.o kappa_cond_3.o numerik.o\
	radsub.o snoein.o chim_gram_3.o dmpeqs.o fich_osc_3.o leq.o omega.o\
	reac_c_3.o static_m_3.o clmnew.o dnunl_3.o gausdn.o lim_gong1_3.o omegac.o\
	read_osc_3.o store.o colatm_3.o edding_3.o gausdp_g.o lim_tau1_3.o\
	opa_gong_3.o storec.o colloc.o eq_diffus_3.o gmass.o lim_tau1_dc.o\
	opa_zams_he.o resout_cephee1.o strans.o colpnt.o eqatm_3.o herm_etat.o\
	lim_teff_3.o opal96_tops99.o rk_imps_3.o sum_n.o conv_gong_3.o eqstf.o\
	hheion.o lim_zc_3.o opal96_tops99_2000.o ro_new.o t6rinterp.o\
	conv_jdc.o eqstfc.o hmnion.o list_3.o opal96_tops99_2001.o sbsp1d.o taueff_3.o\
	conv_jmj_3.o eqstp.o hopf_3.o lit_nl_cephee1.o opal96_tops99_exSMC.o sbsp1dn.o\
	taux_reac.o coulombln.o eqstpc.o horder.o long.o opal96_tops99_fixe01.o sbsp_dis.o\
	taux_tab.o cres.o horner.o minv.o perte_ext.o sbval.o theffp.o ctes_3.o\
	hsra_3.o neuville.o perte_masse.o schu58.o thermo_3.o ctes_ceph.o\
	etat_ceff_3.o hviona.o neweta.o perte_tot.o schu58_n.o update_3.o des_cesam_3.o\
	etat_eff_3.o hvionac.o newsp1.o phder.o scopy.o zero.o des_m.o etat_eff_dc.o\
	iben_3.o newspn.o phderc.o sdot.o zeroc.o des_r.o etat_eff_ps_3.o int_roxburgh.o\
	newton.o pminmax.o setcns.o zone_3.o edding_TESTp20.o edding_TESTm20.o\
	conv_cm_new_3.o conv_cgm_new_3.o conv_cm_3.o nrj_nacre.o taux_nacre.o\
	etat_opal_dc.o esac_dc.o opal96_tops99_2003.o
#
#===================================================
#
base_her_cte.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	base_her_cte.f
	$(FC) $(OPT) base_her_cte.f
#
bilin.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	bilin.f
	$(FC) $(OPT) bilin.f
#
bilinc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	bilinc.f
	$(FC) $(OPT) bilinc.f
#
bleqstc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	bleqstc.f
	$(FC) $(OPT) bleqstc.f
#
blk_eos.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	blk_eos.f
	$(FC) $(OPT) blk_eos.f
#
blstio.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	blstio.f
	$(FC) $(OPT) blstio.f
#
bspint.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	bspint.f
	$(FC) $(OPT) bspint.f
#
bval.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	bval.f
	$(FC) $(OPT) bval.f
#
bvald.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	bvald.f
	$(FC) $(OPT) bvald.f
#
chim_gram_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	chim_gram_3.f
	$(FC) $(OPT) chim_gram_3.f
#
clmnew.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	clmnew.f
	$(FC) $(OPT) clmnew.f
#
colatm_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	colatm_3.f
	$(FC) $(OPT) colatm_3.f
#
colloc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	colloc.f
	$(FC) $(OPT) colloc.f
#
colpnt.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	colpnt.f
	$(FC) $(OPT) colpnt.f
#
conv_gong_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	conv_gong_3.f
	$(FC) $(OPT) conv_gong_3.f
#
conv_jdc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	conv_jdc.f
	$(FC) $(OPT) conv_jdc.f
#
conv_jmj_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	conv_jmj_3.f
	$(FC) $(OPT) conv_jmj_3.f
#
coulombln.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	coulombln.f
	$(FC) $(OPT)  coulombln.f
#
cres.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	cres.f
	$(FC) $(OPT) cres.f
#
ctes_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ctes_3.f
	$(FC) $(OPT) ctes_3.f
#
ctes_ceph.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ctes_ceph.f
	$(FC) $(OPT) ctes_ceph.f
#
des_cesam_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	des_cesam_3.f
	$(FC) $(OPT) des_cesam_3.f
#
des_dasgal1.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	des_dasgal1.f
	$(FC) $(OPT) des_dasgal1.f
#
des_m.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	des_m.f
	$(FC) $(OPT) des_m.f
#
des_r.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	des_r.f
	$(FC) $(OPT) des_r.f
#
desfull_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	desfull_3.f
	$(FC) $(OPT) desfull_3.f
#
dess_hr_DC.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	dess_hr_DC.f
	$(FC) $(OPT) dess_hr_DC.f
#
dess_ol_dig01.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	dess_ol_dig01.f
	$(FC) $(OPT) dess_ol_dig01.f
#
dgrad.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	dgrad.f
	$(FC) $(OPT) dgrad.f
#
difdiv.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	difdiv.f
	$(FC) $(OPT) difdiv.f
#
diff_cte_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	diff_cte_3.f
	$(FC) $(OPT) diff_cte_3.f
#
diff_micro_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	diff_micro_3.f
	$(FC) $(OPT) diff_micro_3.f
#
diff_z16.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	diff_z16.f
	$(FC) $(OPT) diff_z16.f
#
diffus_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	diffus_3.f
	$(FC) $(OPT) diffus_3.f
#
dmpeqs.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	dmpeqs.f
	$(FC) $(OPT) dmpeqs.f
#
dnunl_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	dnunl_3.f
	$(FC) $(OPT) dnunl_3.f
#
edding_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	edding_3.f
	$(FC) $(OPT) edding_3.f
#
eq_diffus_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	eq_diffus_3.f
	$(FC) $(OPT) eq_diffus_3.f
#
eqatm_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	eqatm_3.f
	$(FC) $(OPT) eqatm_3.f
#
eqstf.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	eqstf.f
	$(FC) $(OPT) eqstf.f
#
eqstfc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	eqstfc.f
	$(FC) $(OPT) eqstfc.f
#
eqstp.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	eqstp.f
	$(FC) $(OPT) eqstp.f
#
eqstpc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	eqstpc.f
	$(FC) $(OPT) eqstpc.f
# OPAL EOS perso.:
esac_dc.o: cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	esac_dc.f
	$(FC) $(OPT) esac_dc.f
#
etat_opal_dc.o: cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_opal_dc.f
	$(FC) $(OPT) etat_opal_dc.f
#
etat_ceff_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_ceff_3.f
	$(FC) $(OPT) etat_ceff_3.f
#
etat_eff_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_eff_3.f
	$(FC) $(OPT) etat_eff_3.f
#
etat_eff_dc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_eff_dc.f
	$(FC) $(OPT) etat_eff_dc.f
#
etat_eff_ps_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_eff_ps_3.f
	$(FC) $(OPT) etat_eff_ps_3.f
#
etat_gong1_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_gong1_3.f
	$(FC) $(OPT) etat_gong1_3.f
#
etat_gong2_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_gong2_3.f
	$(FC) $(OPT) etat_gong2_3.f
#
etat_int_herm_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_int_herm_3.f
	$(FC) $(OPT) etat_int_herm_3.f
#
etat_opal.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	etat_opal.f
	$(FC) $(OPT) etat_opal.f
#
evol_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	evol_3.f
	$(FC) $(OPT) evol_3.f
#
f4der.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	f4der.f
	$(FC) $(OPT) f4der.f
#
f4mhd.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	f4mhd.f
	$(FC) $(OPT) f4mhd.f
#
f4n.o : cesam_3.parametres atmosphere_3.common ctephy.common\
	evol_chim_3.common modele_3.common f4n.f
	$(FC) $(OPT) f4n.f
#
ferdir.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ferdir.f
	$(FC) $(OPT) ferdir.f
#
fich_osc_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	fich_osc_3.f
	$(FC) $(OPT) fich_osc_3.f
#
gain_masse.o: cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	gain_masse.f
	$(FC) $(OPT) gain_masse.f
#
gausdn.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	gausdn.f
	$(FC) $(OPT) gausdn.f
#
gausdp_g.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	gausdp_g.f
	$(FC) $(OPT) gausdp_g.f
#
gmass.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	gmass.f
	$(FC) $(OPT) gmass.f
#
herm_etat.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	herm_etat.f
	$(FC) $(OPT) herm_etat.f
#
hheion.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	hheion.f
	$(FC) $(OPT) hheion.f
#
hmnion.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	hmnion.f
	$(FC) $(OPT) hmnion.f
#
hopf_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	hopf_3.f
	$(FC) $(OPT) hopf_3.f
#
horder.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	horder.f
	$(FC) $(OPT) horder.f
#
horner.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	horner.f
	$(FC) $(OPT) horner.f
#
hsra_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	hsra_3.f
	$(FC) $(OPT) hsra_3.f
#
hviona.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	hviona.f
	$(FC) $(OPT) hviona.f
#
hvionac.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	hvionac.f
	$(FC) $(OPT) hvionac.f
#
iben_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	iben_3.f
	$(FC) $(OPT) iben_3.f
#
int_roxburgh.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	int_roxburgh.f
	$(FC) $(OPT) int_roxburgh.f
#
int_roxburgh_m.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	int_roxburgh_m.f
	$(FC) $(OPT) int_roxburgh_m.f
#
inteff.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	inteff.f
	$(FC) $(OPT) inteff.f
#
inteffc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	inteffc.f
	$(FC) $(OPT) inteffc.f
#
inter_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	inter_3.f
	$(FC) $(OPT) inter_3.f
#
intertot_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	intertot_3.f
	$(FC) $(OPT) intertot_3.f
#
intgauss.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	intgauss.f
	$(FC) $(OPT) intgauss.f
#
k_5750.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	k_5750.f
	$(FC) $(OPT) k_5750.f
#
k_5750_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	k_5750_3.f
	$(FC) $(OPT) k_5750_3.f
#
kappa_cond_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	kappa_cond_3.f
	$(FC) $(OPT) kappa_cond_3.f
#
leq.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	leq.f
	$(FC) $(OPT) leq.f
#
lim_gong1_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	lim_gong1_3.f
	$(FC) $(OPT) lim_gong1_3.f
#
lim_tau1_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	lim_tau1_3.f
	$(FC) $(OPT) lim_tau1_3.f
#
lim_tau1_dc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	lim_tau1_dc.f
	$(FC) $(OPT) lim_tau1_dc.f
#
lim_teff_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	lim_teff_3.f
	$(FC) $(OPT) lim_teff_3.f
#
lim_zc_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	lim_zc_3.f
	$(FC) $(OPT) lim_zc_3.f
#
lim_zc_3_save.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	lim_zc_3_save.f
	$(FC) $(OPT) lim_zc_3_save.f
#
list_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	list_3.f
	$(FC) $(OPT) list_3.f
#
lit_nl_cephee1.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	lit_nl_cephee1.f
	$(FC) $(OPT) lit_nl_cephee1.f
#
long.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	long.f
	$(FC) $(OPT) long.f
#
minv.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	minv.f
	$(FC) $(OPT) minv.f
#
neuville.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	neuville.f
	$(FC) $(OPT) neuville.f
#
neweta.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	neweta.f
	$(FC) $(OPT) neweta.f
#
newsp1.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	newsp1.f
	$(FC) $(OPT) newsp1.f
#
newspn.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	newspn.f
	$(FC) $(OPT) newspn.f
#
newton.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	newton.f
	$(FC) $(OPT) newton.f
#
no_des_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	no_des_3.f
	$(FC) $(OPT) no_des_3.f
#
noedif.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	noedif.f
	$(FC) $(OPT) noedif.f
#
noeu_dis.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	noeu_dis.f
	$(FC) $(OPT) noeu_dis.f
#
noeud.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	noeud.f
	$(FC) $(OPT) noeud.f
#
nrj_nuc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	nrj_nuc.f
	$(FC) $(OPT) nrj_nuc.f
#
nrj_nacre.o: cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	nrj_nacre.f
	$(FC) $(OPT) nrj_nacre.f
#
nrj_nuc_Xcno.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	nrj_nuc_Xcno.f
	$(FC) $(OPT) nrj_nuc_Xcno.f
#
nrj_nuc_cn.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	nrj_nuc_cn.f
	$(FC) $(OPT) nrj_nuc_cn.f
#
nuc_gong_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	nuc_gong_3.f
	$(FC) $(OPT) nuc_gong_3.f
#
numerik.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	numerik.f
	$(FC) $(OPT) numerik.f
#
omega.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	omega.f
	$(FC) $(OPT) omega.f
#
omegac.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	omegac.f
	$(FC) $(OPT) omegac.f
#
opa_gong_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opa_gong_3.f
	$(FC) $(OPT) opa_gong_3.f
#
opa_zams_he.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opa_zams_he.f
	$(FC) $(OPT) opa_zams_he.f
#
opal96_tops99.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opal96_tops99.f
	$(FC) $(OPT) opal96_tops99.f
#
opal96_tops99_2000.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opal96_tops99_2000.f
	$(FC) $(OPT) opal96_tops99_2000.f
#
opal96_tops99_2001.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opal96_tops99_2001.f
	$(FC) $(OPT) opal96_tops99_2001.f
#
opal96_tops99_2003.o: cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opal96_tops99_2003.f
	$(FC) $(OPT) opal96_tops99_2003.f
#
opal96_tops99_exSMC.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opal96_tops99_exSMC.f
	$(FC) $(OPT) opal96_tops99_exSMC.f
#
opal96_tops99_fixe01.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	opal96_tops99_fixe01.f
	$(FC) $(OPT) opal96_tops99_fixe01.f
#
perte_ext.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	perte_ext.f
	$(FC) $(OPT) perte_ext.f
#
perte_masse.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	perte_masse.f
	$(FC) $(OPT) perte_masse.f
#
perte_tot.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	perte_tot.f
	$(FC) $(OPT) perte_tot.f
#
phder.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	phder.f
	$(FC) $(OPT) phder.f
#
phderc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	phderc.f
	$(FC) $(OPT) phderc.f
#
pminmax.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	pminmax.f
	$(FC) $(OPT) pminmax.f
#
polyder.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	polyder.f
	$(FC) $(OPT) polyder.f
#
pp1dn.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	pp1dn.f
	$(FC) $(OPT) pp1dn.f
#
pp2d.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	pp2d.f
	$(FC) $(OPT) pp2d.f
#
ppcno3a18.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ppcno3a18.f
	$(FC) $(OPT) ppcno3a18.f
#
ppcno3a_12_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ppcno3a_12_3.f
	$(FC) $(OPT) ppcno3a_12_3.f
#
ppcno3a_12e.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ppcno3a_12e.f
	$(FC) $(OPT) ppcno3a_12e.f
#
ppcno3a_9_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ppcno3a_9_3.f
	$(FC) $(OPT) ppcno3a_9_3.f
#
quad.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	quad.f
	$(FC) $(OPT) quad.f
#
radsub.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	radsub.f
	$(FC) $(OPT) radsub.f
#
reac_c_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	reac_c_3.f
	$(FC) $(OPT) reac_c_3.f
#
read_osc_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	read_osc_3.f
	$(FC) $(OPT) read_osc_3.f
#
resout_cephee1.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	resout_cephee1.f
	$(FC) $(OPT) resout_cephee1.f
#
rk_imps_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	rk_imps_3.f
	$(FC) $(OPT) rk_imps_3.f
#
ro_new.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	ro_new.f
	$(FC) $(OPT) ro_new.f
#
sbsp1d.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sbsp1d.f
	$(FC) $(OPT) sbsp1d.f
#
sbsp1dn.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sbsp1dn.f
	$(FC) $(OPT) sbsp1dn.f
#
sbsp_dis.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sbsp_dis.f
	$(FC) $(OPT) sbsp_dis.f
#
sbval.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sbval.f
	$(FC) $(OPT) sbval.f
#
schu58.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	schu58.f
	$(FC) $(OPT) schu58.f
#
schu58_n.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	schu58_n.f
	$(FC) $(OPT) schu58_n.f
#
scopy.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	scopy.f
	$(FC) $(OPT) scopy.f
#
sdot.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sdot.f
	$(FC) $(OPT) sdot.f
#
setcns.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	setcns.f
	$(FC) $(OPT) setcns.f
#
setcnsc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	setcnsc.f
	$(FC) $(OPT) setcnsc.f
#
setf4.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	setf4.f
	$(FC) $(OPT) setf4.f
#
sgauss.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sgauss.f
	$(FC) $(OPT) sgauss.f
#
shell.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	shell.f
	$(FC) $(OPT) shell.f
#
simq.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	simq.f
	$(FC) $(OPT) simq.f
#
slinf.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	slinf.f
	$(FC) $(OPT) slinf.f
#
smatpr.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	smatpr.f
	$(FC) $(OPT) smatpr.f
#
sn1dn_dis.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sn1dn_dis.f
	$(FC) $(OPT) sn1dn_dis.f
#
snoein.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	snoein.f
	$(FC) $(OPT) snoein.f
#
static_m_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	static_m_3.f
	$(FC) $(OPT) static_m_3.f
#
store.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	store.f
	$(FC) $(OPT) store.f
#
storec.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	storec.f
	$(FC) $(OPT) storec.f
#
strans.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	strans.f
	$(FC) $(OPT) strans.f
#
sum_n.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	sum_n.f
	$(FC) $(OPT) sum_n.f
#
t6rinterp.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	t6rinterp.f
	$(FC) $(OPT) t6rinterp.f
#
taueff_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	taueff_3.f
	$(FC) $(OPT) taueff_3.f
#
taux_reac.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	taux_reac.f
	$(FC) $(OPT) taux_reac.f
#
taux_tab.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	taux_tab.f
	$(FC) $(OPT) taux_tab.f
#
taux_nacre.o:  cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	taux_nacre.f
	$(FC) $(OPT) taux_nacre.f
#
theffp.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	theffp.f
	$(FC) $(OPT) theffp.f
#
thermo_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	thermo_3.f
	$(FC) $(OPT) thermo_3.f
#
update_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	update_3.f
	$(FC) $(OPT) update_3.f
#
zero.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	zero.f
	$(FC) $(OPT) zero.f
#
zeroc.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	zeroc.f
	$(FC) $(OPT) zeroc.f
#
zone_3.o : cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	zone_3.f
	$(FC) $(OPT) zone_3.f
#==============================================
# Routine pour avoir un gradient de T � la Canuto :
conv_cm_new_3.o: conv_cm_new_3.f ctephy.common  modele_3.common
	$(FC) $(OPT) conv_cm_new_3.f
conv_cgm_new_3.o: conv_cgm_new_3.f ctephy.common  modele_3.common
	$(FC) $(OPT) conv_cgm_new_3.f
conv_cm_3.o: conv_cm_3.f
	$(FC) $(OPT) conv_cm_3.f
#==============================================
# Routines exp�rimentales :
edding_TESTp20.o: cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	edding_TESTp20.f
	$(FC) $(OPT) edding_TESTp20.f
#
edding_TESTm20.o: cesam_3.parametres atmosphere_3.common\
	ctephy.common evol_chim_3.common modele_3.common\
	edding_TESTm20.f
	$(FC) $(OPT) edding_TESTm20.f
#
# -----------------------------------------------------------------------------------------------------------
clean:
	rm -f *.o *.a
#
