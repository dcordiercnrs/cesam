c****************************************************************************
 
	subroutine cres(etat,opa,conv,nuc,lim_ext,tdetau,
     1	coeff_diff,perte_m,ctes,des_cesam)


c----------------------------------------------------------------------------
c Modifications apportees par D. Cordier  :

c       * ecritures dans le fichier '*.HR' : chercher le mot cle 
c         'lehr'

c       * variante dans le calcul du terme en 'TdS' : chercher le mot cle
c         'modtds'       (Avril 98)

c       * nouvelle conditions d'arret de l'evolution (sur Ycentrale et
c	  la luminosite) : chercher le mot cle 
c         'modifarret'   (Nov. 98)

c       * modifications concernant la possibilite de faire des enregistrements
c         de modeles intermediaires : chercher le mot cle 
c         'modifrep'     (Nov. 98)

c       * question concernant les "pauses" durant le lancement du code :
c         si la reponse est 'o' ---> il y a des "pauses" pendant le
c         lancement est on ne peut pas lancer le programme en batch si
c         le reponse est differentes de 'o' ---> pas de pauses on
c         peut utiliser un batch, voir aussi 'modele_3.common'
c         Chercher le mot cle :
c         'despau' (Fev. 99)

c       * enregistrement des info. concernant Roxburgh dans un fichier
c         'Roxy' (Fev. 99)

c       * enregistrement dans le fichier '*.HR' de la temperature (en log10)
c         au bord de zone convective avec et sans eventuel over/undershoot
c         Chercher le mot cle :
c         'tzc'

c       * création d'un nouveau fichier de sortie : '*.struc' avec la structure
c         de l'étoile enregistrée de façon plus claire que dans '*.osc'. Ce fichier
c         est directement utilisable avec un plotter.
c         D.C., spetembre 2003, Varsovie.

c----------------------------------------------------------------------------


c	calcul de modele de structure interne et d'evolution stellaire
 
c	methode collocation de de Boor p. 277 ou Galerkin
c	adapte a la structure interne
 
c	avec fonction de repartition
 
c	Auteurs : P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c                 D. Cordier, DASGAL, Observatoire de Paris - Meudon
c                 Version du 2 Novembre 1998
 
c	CESAM version 3
 
c	MODIFS
c	07 10 96 gravite effective dans vaissala
c	10 10 96 introduction de la rotation dans thermo
c	21 03 97 bug dans calcul de z_opa (signale par Sacha)
c	08 05 97 : introduction de l'option kipp pour calcul du TdS
 
c entrees:
c	noms des routines de physique
 
c---------------------------------------------------------------------
 
c routines externes de physique
c	etat: equation d'etat
c	opa: opacite
c	conv: convection
c	nuc: reactions thermonucleaires
c	lim_ext: condition limite externe
c	tdetau: loi t(tau)
c	coeff_diff: calcul du coefficient de diffusion
c	ctes: definition des constantes
c	coeff_diff: coefficients de diffusion
c	perte_m : perte de masse
 
c fonctions des principales routines externes numeriques
c	cesam_3 : gestion generale des calculs, des ecritures
c	des_cesam: dessin du modele durant l'execution
 
c fonctions des principales routines internes numeriques
c	resout_3 : initialisation, formation, resolution du probleme aux limites
c	evol_3 : gestion de l'evolution de la composition chimique
c	static_m_3 : formation des equations de l'equilibre quasi statique
c	lim_zc_3 : determination des limites ZR/ZC
c	update_3 : transfert t+dt ---> t ---> t-dt ---> t-dt_-dt__
c	repartit_3 : initialisation des repartitions
c	diffus_3, eq_diffus_3 : calcul de la diffusion
c	coll_atm_3, eq_atm_3 : restitution d'atmosphere
 
c------------------------------------------------------------------------
 
c	Les COMMONs ne contiennent que des quantites constantes
c	au cours de l'evolution (elles peuvent changer lors de l'initialisation)
 
c COMMON /ctephys/	initialise dans ctes
c	constantes physiques utilisees
 
c COMMON /fich_atm/	initialise cesam_3
c	sert a passer le nom du modele au fichier atmosphere
 
c COMMON /gong/	initialise dans cesam_3, lim_gong1_3
c	sert a passer les parametres beta et lambda pour les comparaisons
c	du SMC
 
c COMMON /modele_3/ initialise dans lit_nl_cephee2, cesam_3, resout_3
c	caracteristiques du modele
c	mtot : masse totale
c	alpha : longueur de melange
c	dtmax, dtmin : pas temporel maximal, minimal
c	dtlist : intervalle de temps entre deux listes completes
c	modele : identificateur du modele
c	x0, y0, z0 : abondance initiale de H, He4, elements lourds
c	agemax : age maximal
c	dt0 : pas temporel initial
c	precix : precision en abscisse
c	precit : precision integration temporelle
c	d_grav : variation max du TdS
c	vsal : parametre de convection
c	ovshts, ovshti : parametres d'overshoot
c	w_rot : vitesse angulaire initiale
c	ne,m_qs,m_ch,mpr : nb equations, ordre des splines pour
c	quasi-statique et chimie, mpr=m_qs+1
c	opa1, opa2, opa3, opa4 : noms de fichiers d'opacite
c	m_ch : ordre des splines d'interpolation pour la comp. chim.
c	m_qs+1=mpr: ordre des splines pour integration spatiale
c	ne: nombre d'equations
c	modele,opa1......opa8 : nom du modele et des fichiers d'opacite
c	ledoux=.true.: on utilise le critere de Ledoux
c	jpz=.true.: convection penetratrive de JPZ
c	der_num: derivees numeriques
c	discon=.true. : discontinuite de comp. chim. aux limites ZR/ZC
c	rot_solid=.true.: rotation solide
c	z_cte=.true. : Z est fixe
 
c COMMON /evol_chim_3/	initialise dans lit_nl_cephee1, resout_3, nuc
c	donnees pour les reactions thermonucleaires
c	nucleo(pnelem) : masses atomiques des elements
c	ab_min : abondances minimale
c	t_inf : T < t_inf plus de reactions thermonucleaires
c	mdot : taux de perte de masse
c	nbelem : nombre d'especes chimiques generalisees
c	nchim : nombre d'especes chimiques
c	nom_elem : noms des especes chimiques generalisees
c	nreac : nombre de reactions nucleaires
c	iw, iz, ihe4 : indice du moment angulaire, de Z, de He4
c	theta : ordre d'integration par RK
c	diffusion=.true. : on tient compte de la diffusion
 
c COMMON /premain/	initialise dans	iben_3
c		transporte la constante d'Iben pour la PMS
c	c_iben: constante d'Iben
 
c COMMON /atmosphere_3/	initialise dans lit_nl et dans lim_ext
c		elements pour le calcul de l'atmosphere
c	tau_min: ep. opt. a l'exterieur
c	tau_max: ep. opt. au fond
c	n_atm: nb. de couches pour l'atm.
c	n23: indice de la couche ou le rayon est defini
c	nea: nombre d'equations dans l'atmosphere
 
c	les COMMONs atmosphere_3, modele_3
c	ctephy, reac_nuc sont definis dans les fichiers
c	atmosphere_3.common, modele.common
c	ctephy.common, inputs_3.common, evol_chim_3.common
c	ils sont introduits par des INCLUDEs
 
c	les parametres sont definis dans le fichier cesam_3.parametres
c	ils sont introduits par un INCLUDE
 
c------------------------------------------------------------------
 
c dimensions maximales (en parametre) :
c	pn : nombre maximal de points de raccord
c	pne : nombre maximal d'equations d'ordre 1
c	pm_qs : ordre maximal des splines pour la collocation
c	pm_ch : ordre maximal des splines pour l'interpolation de la comp.chim.
c	pnzc : nombre maximal de zones convectives
c	pnelem : nombre maximal d'elements diffuses
c	pn_atm : nombre maximal de couches dans l'atmosphere
c	pnreac : nombre maximal de reactions nucleaires
 
c	le centre est a la couche 1, la limite externe en couche n
c	le nombre de couches est VARIABLE il est determine dans le SSP repartit
c	une couche est place a chaque limite ZR/ZC (a 5% pres)
 
c	increment temporel :
 
c		dt=delta t / 10**6 ans
 
c	variable independante : q l'indice (Reel) de couche
 
c	variables dependantes :
 
c		ksi=ln P pression cgs
c		eta=ln T temperature K
c		dzeta=(r/rsol)**2
c		lambda=(l/lsol)**2/3
c		mu=(m/mtot)**2/3-1
c		psi=dQ/dq=cte avec Q(mu,t)=fonction d'espacement
 
c		xchim(i),i=1,nbelem : abondances / mole en f(m/Mtot)
c		 + MA + Z  (MA : moment angulaire par unite de masse w r**2)
 
c	ordre et indices d'identification des elements:
c	H D He3 He4 Li7 Be7 C12 C13 N14 N15 O16 O17 MA Z
c	1 2  3   4   5   6   7   8   9  10   11  12 13 14	
c	avant le debut d'un pas temporel:
 
c	mc_t,mct_t,nc_t,knotc_t,chim_t:   comp. chim. a age-dt
c	mc,mct,nc,knotc,chim: comp. chim. a age
c	bp,q,qt,knot,p,t,r,l,m:  var. pples. a age
c	m23,m23t_t,m_ch,n,knot23_t,tdst: TdS a age t
 
c_______________________________________________________________________
 
c	Fichiers utilises
 
c	les fichiers de sortie ont le nom choisi par l'utilisateur nom_fich
c	suivit d'un attribut qui les distingue
c	l'attribut "_B" est specifique aux fichiers binaires
 
c	FOR002	listing de sortie, cree dans cesam ferme et garde dans
c		cesam, nom : nom_fich.lis
 
c	FOR003	lecture des donnees, utilise temporairement dans resout3
c		nom : cesam_3.don
 
c	FOR004  modele initial ou point de reprise (binaire ou en ASCII),
c		utilise temporairement dans resout_3
 
c	FOR011 lecture des tables d'opacite
c		utilise temporairement dans routines d'opacite
 
c	FOR011-->017 lecture des tables d'equation d'etat
c		utilise temporairement dans certains SSP d'equation d'etat
 
c	FOR020	tables de reactions nucleaires, ouvert et ferme dans reac_t_3
 
c	FOR024  modele initial d'age 0 homogene, ZAMS ou PMS,
c		ouvert et ferme dans cesam_3
c		nom : nom_fich_B.hom en binaire
 
c	FOR025  dernier modele, ouvert et ferme dans cesam_3
c		nom : nom_fich_B.dat en binaire
 
c	FOR026  point de reprise, ouvert et ferme dans cesam_3
c		nom : nom_fich_B.rep en binaire
 
c	FOR030  liste pour oscillation ouvert et ferme dans cesam_3
c		nom : nom_fich.osc en ASCII
 
c	FOR031	dernier modele calcule de l'atmosphere
c		ouvert et ferme dans lim_teff
c		nom : nom_fich_B.atm en binaire

c       FOR032 le fichier de structure de D.C.

c	FOR048  liste pour dessins ouvert et ferme dans cesam_3
c		nom : nom_fich.des en ASCII
c		des interventions sont necessaires dans cesam_3 pour sa
creation:
c		enlever les "c" de commentaires
 
c	FOR053	fichier pour tracer un diag. HR et les ZC ouvert et ferme dans
cesam_3
c	        nom : nom_fich.HR en ASCII
 
c	FOR060	tables d'EOS de LAOL
	
c-----------------------------------------------------------------------
 
c		P R O G R A M M E   P R I N C I P A L
 
c	declaration de dimensions, gestion du calcul
 
	implicit none
 
c parametres: des maxima pour :
c	nb. de points de liaison, ordre splines pour coll., ordre des equations,
c	nb. d'elements chimiques, nb. d'equations, nb. de couches dans atmosph.
c	nb. zones convectives, ordre des splines d'interp. de la comp. chim.
c	nb. reactions thermonucleaires
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
c	include 'cesam_3.version'
 
	integer n,i,k,j,un23,n_t,knot,i_pp,i_cno,i_3a,i_gr,ipms,knota,
     1	nc_t,knotc_t,nc,knotc,lim_t,jlim_t(2*pnzc),
     2	lx,nadd,knot_t,jlim(2*pnzc),lim,reprise,ipn,
     3	nc_,knotc_,knot23_t,n_,knot_
c	data ipms/3/
 
	integer*4 long
	external long

c---------
c modifrep
 
        integer n_fich_rep ! Modif. DC Nov. 98

	logical enre_bin

	real*8 dtrep

	character key_word*5
 
c-----------
c modifarret
 
        logical lcroissant ! Modif. DC Nov. 98
 
        real*8 log_lslsol, y_stop ! idem

c---------
c lehr
        integer indi, n_pp, n_cno, n_3a, n_grav_p,
     +          n_grav_m, change_compo, num_pp(10), num_cno(10),
     +          num_3a(10), num_grav_p(10), num_grav_m(10)
 
        logical lpp, lcno, l3a, lgrav_p, lgrav_m, natpp(10),
     +          natcno(10), nat3a(10), natgrav_p(10), natgrav_m(10),
     +          l_c
 
	real*8 pos_lpp(10), pos_lcno(10), pos_l3a(10),
     +  pos_lgrav_p(10), pos_lgrav_m(10),
     +  t_lpp(10), t_lcno(10), t_l3a(10), t_lgrav_p(10),t_lgrav_m(10),
     +  rho_lcno(10), x_lcno(10), y_lcno(10), Xcno_lcno(10),
     +  l_xchi(pnelem), bid1, bid2, bid3, bid4, bid5, bid6,
     +  bid7, bid8, bid9, bid10, bid11, bid12, bid13, bid14,
     +  bid15, bid16, bid17,
     +  lambda_c, z_local, z_lim, m_change_compo

	real*8 delrad, xxx, yyy

        real*8 log_lum_c, masse_c, rayon_c, log_phi_c, xsurf
c----------
c modtds
	real*8 dlum, dlum_s_dm
	logical tdstest

	common/tdstest/tdstest

c----------
c tzc
	real*8 logtzc(2*pnzc), logrozc(2*pnzc), logtov(2*pnzc)

	common /temp_limZCZR/ logtzc, logrozc, logtov

c----------
	real*8	xchim(pnelem),xchim1(pnelem),dxchim(pnelem),dcompg(pnelem*pn),
     1	p(pn),t(pn),r(pn),l(pn),m(pn),w(pn),z(pn),gradient(pn),ro(pn),
     2	grad_mj(pn),d_grad(pn),compchim(pnelem*pn),eint(pn),
     3	epsilon(5*pn),kap(pn),nh1(pn),bidt(pn),compg(pnelem*pn),
     4	nhe1(pn),nhe2(pn),beta(pn),delta(pn),hp(pn),gradad(pn),	
     4	gradrad(pn),gradconv(pn),alfa(pn),dcapdr(pn),dcapdt(pn),
     5	depsdr(pn),depsdt(pn),anupp(pn),anupep(pn),anub8(pn),
     6	anube7(pn),anun13(pn),anuo15(pn),dcompchim(pnelem*pn),
     7	vaissala(pn),dlpdxt(pn),dlpdxr(pn),cp(pn),xdot(pn),ydot(pn),
     8	rdot(pn),lambda(pn)
	
	real*8	drop,drot,dup,dut,d2p,d2ro,dgradp,dgradt,dgradl,dgradr,
     1	dgradm,epsilo(5),depsp,depst,dkapp,dkapt,drox,dux,
     3	dgradx,depsx(pnelem),dkapx,aradias3,dcomp(pnelem),pertem,
     4	lteff,teff,agep,y,deltap,deltat,deltax,dcpp,dcpt,dcpx,
     4	drott,drotp,drotx,dutt,dutp,dutx,ro_t,u_t,bid,	
     5	dm,e_tot,e_pp,e_cno,e_3a,e_gr,p_t,t_t,nuc_m(pnelem),
     6	f(pnelem),dfdx(pnelem),tdst(pn),m23t_t(pqt),
     7	fonc(0:5),poly(0:5),absc(0:5)	
	
	real*8	p_atm(pn_atm),t_atm(pn_atm),r_atm(pn_atm),m_atm(pn_atm),
     1	tau(pn_atm),ro_atm(pn_atm),k_atm(pn_atm),gradr_atm(pn_atm),
     2	grada_atm(pn_atm),grad_atm(pn_atm),chim1g(pnelem),
     3	nh1_atm(pn_atm),nhe1_atm(pn_atm),nhe2_atm(pn_atm),
     4	alfa_atm(pn_atm),delta_atm(pn_atm),cp_atm(pn_atm),
     5	vais_atm(pn_atm),gradc_atm(pn_atm),ldcapdr_a(pn_atm),
     6	ldcapdt_a(pn_atm),lnp_a(pn_atm),lnt_a(pn_atm),grad_mj_a(pn_atm),
     7	rstar,pext,text,dpdl,dpdr,dtdl,dtdr,mext,dml,dmr,
     8	epsa(5),u,lamb,logg,hpa,rp1,rp2,lp1,lp2,t_inf0,tds
	
	real*8	cte3,cte5,cte7,cte8
 
	real*8	bp(pbp),q(pn),qt(pqt),chim(pchim),mc(pnch),mct(pchimt),
     1	r2(pn),m23(pn),m_zc(2*pnzc),r_zc(2*pnzc),r_ov(2*pnzc),
     2	age,bp_t(pbp),q_t(pn),qt_t(pqt),chim_t(pchim),mc_t(pnch),
     3	mct_t(pchimt),r2_t(pn),m23_t(pn),dt,m_zc_t(2*pnzc),
     4	r_zc_t(2*pnzc),r_ov_t(2*pnzc),bp_(pbp),q_(pn),qt_(pqt),
     5	chim_(pchim),mc_(pnch),mct_(pchimt),r2_(pn),m23_(pn),
     6	dt_,mstar,mstar_t,chimg(pnelem),dchimg(pnelem),dw(pn)
	
	integer iconst,ivar,ivers	!pour ecrire fichier oscillations
	real*8 glob(15),var(20+pnelem,pn+pn_atm)
 
	real*8 log_teff,x_stop,t_stop		!pour test d'arret d'evolution
	character*4 arret

c------------
c Roxy

	real*8 diff, mroxSmsch

	common/roxburgh/diff, age, mroxSmsch
 
c------------
c modifarret
 
	common/arret/log_teff,log_lslsol,y_stop,x_stop,t_stop,arret,
     +               lcroissant ! Modif. DC Nov. 98

c____________
c modifrep

	common/enrep/dtrep,enre_bin,key_word ! idem

	real*8 beta_g,lambda_g
	common/gong/beta_g,lambda_g
 
	real*8 c_iben
	common/premain/c_iben


	character*1 oui
	character*2 c_qs,ci,cj,l2,d_qs
	character*3 cr,ch,l3
	character*4 cm,l4,ck
	character*9 today
	character*30 chaine
	character*50 nom_fich1,modeleb
	character*60 methode
	character*80 titre
 
	logical convec(pn),llist,statut,lconv(2*pnzc),lconv_t(2*pnzc),
     1	pms,pms_pr,zams,zams_pr,post_e,zams_e,passe,dessin,post,
     2	post_pr,cohe,cohe_e,cohe_pr,sort,ecrit
c	data passe/.false./
 
	external ctes,etat,opa,conv,nuc,lim_ext,tdetau,perte_m,
     1	coeff_diff,iben_3,des_cesam
	
	include 'cesam_3.version'
	data ipms/3/
	data passe/.false./
 
	save
 
c-----------------------------------------------------------------------
 
2000	format((1x,1p8d10.3))
 
	beta_g=1.	!pour GONG1 beta et lambda seront differents
	lambda_g=1.
 
c	calcul a effectuer, point de reprise
 
	write(6,*)' '
	write(6,*)'CESAM, version: ',version
	write(6,*)' '
199	print*,'Si on poursuit une evolution : taper 1 puis RETURN'
	write(6,*)' '
	print*,'Si on initialise un modele de ZAMS : taper 2 puis RETURN'
	write(6,*)' '
	print*,'Si on initialise un modele de PMS : taper 3 puis RETURN'
	
c	on entre un23 = 1, 2 ou 3 comme indiquer ci dessus
c	puis, suivant les cas on donne a un23 les valeurs suivantes:
	
c	un23 :  1 poursuite d'une evolution, modele repris en binaire
c		2 modele initial en binaire, ZAMS
c		3 modele initial en binaire, PMS
c	       -2 modele initial en ASCII, ZAMS
c	       -3 modele initial en ASCII, PMS
 
	read(5,*)un23
	

c despau

	rep_pause='o'
	write(6,*) '*****************************************'
	write(6,*) 'Accepte-t-on les pauses ? (o/n)'
	read(5,'(a)') rep_pause

	if(un23 .eq. 1)then
	 write(6,6)
6	 format(1x,'entrer le nom du fichier binaire du modele dont on ',
     1	'poursuit l''evolution')
	 read(5,'(a)')nom_fich1
	 write(6,*)'on utilise le modele ',nom_fich1
	 write(6,*)' '
 
c-------------
c modifrep
 
         n_fich_rep=0 ! Modif. DC Nov. 98
 
	 open(unit=4,form='unformatted',status='old',file=nom_fich1)
	 read(4)age,nbelem,nchim,mtot,alpha,modele,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	 close (unit=4)
 
	 write(6,*)'On poursuit l''evolution du modele: ',modele(:long(modele))
	 write(6,201)age
201	 format(1x,' d''age:',1pd10.3)
	 write(6,*)
	 pms=.false.
	 zams=.false.
	 post=.false.
	 cohe=.false.	!debut de la combustion de He	
	
	elseif(un23 .eq. 2)then
	 pms=.false.
	 zams=.true.
	 post=.false.
	 cohe=.false.
	 print*,'le modele initial de ZAMS est-il donne en binaire ? o/n'
	 read(5,'(a)')oui
	 if(oui .ne. 'n')then
	  write(6,7)
7	  format(1x,'entrer le nom du fichier binaire du modele pris pour ',
     1	'modele initial')	
	  read(5,'(a)')nom_fich1
	  write(6,*)'on utilise le modele ',nom_fich1
	  write(6,*)' '
	
	  open(unit=4,form='unformatted',status='old',file=nom_fich1)
	  read(4)age,nbelem,nchim,mtot,alpha,modele,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	  close (unit=4)
 
	  write(6,*)'On prend pour modele initial: ',modele(:long(modele))
	  write(6,*)' '	
	  un23=2	!fichier en binaire (intruction redondante)
	 else
	  write(6,9)
9	  format(1x,'entrer le nom du fichier ASCII du modele pris pour ',
     1	'modele initial')
	  read(5,'(a)')nom_fich1	!on lira ce fichier dans resout_3
	  open(unit=4,form='formatted',status='old',file=nom_fich1)
	  write(6,*)'CESAM utilise le modele ',nom_fich1
	  write(6,*)' '
	  un23=3		!fichier en ASCII
	 endif
	 age=0.d0
	 dt=0.
	
	elseif(un23 .eq. 3)then
	 pms=.true.
	 zams=.false.
	 post=.false.
	 cohe=.false.
	 print*,'le modele initial de PMS est-il donne en binaire ? o/n'
	 read(5,'(a)')oui
	 if(oui .eq. 'o')then
	  write(6,*)'CESAM utilise un modele initial de PMS en binaire'	
	  write(6,7)
	  read(5,'(a)')nom_fich1
	  write(6,*)'CESAM utilise le modele ',nom_fich1
	  write(6,*)' '
	  open(unit=4,form='unformatted',status='old',file=nom_fich1)
	  read(4)age,nbelem,nchim,mtot,alpha,modele,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	  close (unit=4)
	  un23=-2 !fichier en binaire et, pour la PMS, pas de modele homogene
	 else
	  write(6,*)'CESAM utilise un modele initial de PMS en ASCII'
	  write(6,*)' '
	  write(6,9)
	  read(5,'(a)')nom_fich1
	  open(unit=4,form='formatted',status='old',file=nom_fich1)
	  write(6,*)'CESAM utilise le modele ',nom_fich1
	  write(6,*)' '
	  un23=-3  !fichier en ASCII et, pour la PMS, pas de modele homogene
	 endif
	 age=0.d0
	 dt=0.	
	 write(6,*)'calcul de l''evolution pre-sequence principale:'
	 write(6,*)'entrer la constante de contraction C (2.d-2)'
	 read(5,*)c_iben
	 ipms=1
	
c	 ipms =	1, 2, 3 gere le calcul PMS,
c		1 et 2 : premier et second modele d'initialisation
c		3 : le modele ignore qu'il est issu de la pms
	
	else
	 print*,'ERREUR, entrer 1, 2 ou 3, vous avez tape',un23
	 goto199
	endif
 
c	identificateurs des fichiers des donnees, du point de reprise,
c       des resultats,
c	du modele d'age 0, des oscillations
 
	write(6,*)'entrer l''identificateur du modele'
	write(6,*)'exemple: soleil_evolue'
	read(5,'(a)')nom_fich2
	write(6,*)'identificateur des fichiers du modele : ',
     1	nom_fich2(:long(nom_fich2))
	write(6,*)' '
 
c	fichier pour le listing
 
	open(unit=2,form='formatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'.lis')
 
c	fichier pour le diagramme HR
 
	open(unit=53,form='formatted',status='unknown',access='append',
     1	file=nom_fich2(:long(nom_fich2))//'.HR')
	agep=-1.d30	!initialisation pour lister le premier modele
c-------
c Roxy
	open(unit=76,form='formatted',status='unknown',access='append',
     1	file=nom_fich2(:long(nom_fich2))//'.ROX')
	agep=-1.d30	!Fichier pour Roxburgh
 
c	write(6,*)'avant resout pms,ipms,un23,dt',pms,ipms,un23,dt
 
	if(pms)then
	 call resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,nuc,lim_ext,tdetau,coeff_diff,perte_m)
	 t_inf0=t_inf		!initialisation de Iben_3, t, ro fictifs
	 call iben_3(t(1),ro(1),xchim,dcomp,bidt,.false.,1,
     1	epsilo,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt)
c	 write(6,2000)t_inf,t_inf0
c	 pause'entre les deux resout'
	 call resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,iben_3,lim_ext,tdetau,coeff_diff,perte_m)
 
	else
	 call resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,nuc,lim_ext,tdetau,coeff_diff,perte_m)
	endif
 
c	quelques constantes
 
	aradias3=aradia/3.
		
	cte3=4.*pi*rsol**3/msol
	cte5=-4.*pi*g/3.*rsol**2
	cte7=rsol/secon6
	cte8=lsol/pi/rsol**2/aradia/clight
 
c	encode(4,3,cm)int(mtot*100.)
        write(cm,3) (int(mtot*100.))
c	encode(3,1,ch)int(x0*1000.)
        write(ch,1) (int(x0*1000.))
c	encode(3,1,cr)int(alpha*100.)
        write(cr,1) (int(alpha*100.))
1	format(i3)
2	format(i2)
3	format(i4)
 
c	identification de la methode
 
c	encode(2,2,c_qs)m_qs
        write(c_qs,2) (m_qs)
c	encode(2,2,d_qs)m_ch
        write(d_qs,2) (m_ch)
 
	methode=version
 
	methode=methode(:long(methode))//' colloc.'//c_qs
	if(diffusion)then
	 methode=methode(:long(methode))//', diffus.'
	else
	 methode=methode(:long(methode))//', no diffus.'
	endif
	
c	fait-on le dessin sur ecran?	
 
	write(6,*)'dessin du modele en cours de l''evolution ? o/n'
	read(5,'(a)')oui
	dessin=oui .eq. 'o'
	
c	on garde le modele d'age 0 sur le fichier xxxx_B.hom ou xxxx_B.pms
 
	if(abs(un23) .ge. 2)then	!cas d'un modele homogene
	 dt_=0.		!initialisation
	 dt=0.
	 if(pms)then
	  open(unit=24,form='unformatted',status='unknown', !1ier modele PMS
     1	file=nom_fich2(:long(nom_fich2))//'_B.pms')
	  modeleb='m'//cm//'X'//ch//'a'//cr//'_B.pms'
	 else
	  open(unit=24,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.hom')
	  modeleb='m'//cm//'X'//ch//'a'//cr//'_B.hom'
	 endif
	 write(24)age,nbelem,nchim,mtot,alpha,modeleb,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	 close(unit=24)
	endif
 
c	garde des points de reprise
 
c	write(6,*)'faut-il garder tous les modeles intermediaires ? o/n'
c	read(5,'(a)')oui
c	statut=oui .eq. 'o'
	statut=.false.		!elimine
 
c	ecriture des solutions a l'ecran
 
c	write(6,*)'faut-il ecrire les modeles sur l''ecran ? o/n'
c	read(5,'(a)')oui
c	llist=oui .eq. 'o'
	llist=.false.		!elimine
 
c	gestion du premier pas temporel
 
	if(ipms .gt. 2)then	!sauf pour les modeles initiaux PMS
	 dt=0.
	 if(un23 .eq. 1)then	!si reprise d'une evolution
	  if(dt_ .gt. 0.d0)then
	   write(6,111)dt_
111	   format(1x,'pas temporel precedent=',1pd10.3,
     1	' 10**6 ans, le garde-t-on ? o/n')
	   read(5,'(a)')oui
	   if(oui .eq. 'o')then
	    dt=dt_
	   else
	    write(6,*)'entrer le nouveau pas temporel unite: 10**6 ans'
	    read(5,*)dt
	    if(dt .gt. dtmax)then
	     write(6,112)dtmax
112	     format(1x,'dt > dtmax, le pas temporel est ramene a dtmax=',
     1	1pd10.3)
	    endif
	   endif
	  else
	   write(6,113)dt0
113	   format(1x,'utilise-t-on dt0=',1pd10.3,' pour le pas temporel? o/n')
	   read(5,'(a)')oui
	   if(oui .eq. 'o')then
	    dt=dt0
	   else
	    write(6,*)'entrer le nouveau pas temporel unite: 10**6 ans'
	    read(5,*)dt	
	    if(dt .gt. dtmax)then
	     write(6,112)dtmax
	    endif
	   endif
	   if(age+dt .gt. agemax)then
	    write(6,53)age+dt,agemax
	    write(2,53)age+dt,agemax	
53	    format(1x,'age+dt=',1pd10.3,' > agemax =',1pd10.3,
     1	' ajustement de dt')
	    dt=max(0.d0,agemax-age)
	   endif
	   write(2,50)dt
	   write(6,50)dt
50	   format(//,t2,'reprise de l''evolution avec le pas temporel',
     1	1pd10.3,' 10**6 ans',//)
	  endif	!dt_>0
	
	 else		!un23 > 1
	  dt_=0.
	  if(agemax .gt. 0.d0)then
	   write(6,113)dt0
	   read(5,'(a)')oui
	   if(oui .eq. 'o')then
	    dt=dt0
	   else	
	    write(6,*)'entrer le pas temporel unite: 10**6 ans'
	    read(5,*)dt
	   endif		!initialisation d'une evolution
	   if(dt .gt. agemax)then
	    write(6,54)dt,agemax
	    write(2,54)dt,agemax	
54	    format(1x,'dt=',1pd10.3,' > agemax =',1pd10.3,
     1	' ajustement de dt')
	    dt=agemax
	   endif
	   write(2,51)dt
	   write(6,51)dt
51	   format(//,t2,'debut d''evolution avec le pas temporel',
     1	1pd10.3,' 10**6 ans',//)
	
	  else
	   print*,'calcul d''un modele d''age 0'
	   dt=0.
	  endif		!agemax > 0
	 endif		!un23 = 1
	endif		!ipms > 2
 
c	open(unit=48,form='formatted',status='unknown',
c	1	file=nom_fich2(:long(nom_fich2))//'.des')
 
c'''''''''''''''''''''' retour ''''''''''''''''''''''''''''''''''''''''''''''
 
c	fichier du point de reprise
 
	reprise=10		!avec 1, 2, ... , 9 ca ne marche pas ????
 
100	do i=1,n			!extraction p,t,r,l,m
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),lx,f,dfdx)
	 p(i)=exp(f(1))	!variable ln p
	 t(i)=exp(f(2))	!variable ln T
	 r(i)=sqrt(abs(f(3)))	!rayon/Rsol
	 l(i)=sqrt(abs(f(4)))**3	!l/Lsol
	 m(i)=sqrt(abs(f(5)))**3	!m/Msol
	 grad_mj(i)=dfdx(2)/dfdx(1)		!le gradient reel pour Marie Jo
c	 write(6,2000)p(i),t(i),r(i),l(i),m(i),grad_mj(i)
 
c	 la composition chimique
 
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lx,xchim,dxchim)
	
	 do j=1,nbelem
	  chimg(j)=xchim(j)
	  dchimg(j)=dxchim(j)
	 enddo
	 call chim_gram_3(chimg,dchimg,nuc_m)	
	 do j=1,nbelem
	  compchim(nbelem*(i-1)+j)=xchim(j)
	  compg(nbelem*(i-1)+j)=chimg(j)	
	  if(f(5) .gt. 0.d0)then
	   dcompchim(nbelem*(i-1)+j)=dxchim(j)*2.d0/3.d0/sqrt(f(5))	
	   dcompg(nbelem*(i-1)+j)=dchimg(j)*2.d0/3.d0/sqrt(f(5))
	  else
	   dcompchim(nbelem*(i-1)+j)=0
	   dcompg(nbelem*(i-1)+j)=0	
	  endif
	  if(i .eq. n)then
	   xchim1(j)=xchim(j)	!pour l'atmosphere
	   chim1g(j)=chimg(j)	!pour l'atmosphere
	  endif	
	 enddo
	
c	 la rotation
	
	 if(iw .gt. 1 .and. f(3) .gt. 0.)then
	  w(i)=xchim(iw)/f(3)
	  dw(i)=dxchim(iw)/f(3)
	 else
	  w(i)=w_rot
	  dw(i)=0.
	 endif
 
	enddo
 
c	on ecrit la solution a l'ecran
 
	if(llist)then
	 write(6,101)age
101	 format(1x,'age=',1pd10.3)
	 write(6,102)
102	 format(t9,'1-m',t24,'p',t36,'t',t49,'r',t60,'l',t71,'X')
	endif
	
	if(iw .gt. 1)then	!rotation au centre
	 w(1)=w(2)
	 dw(1)=dw(2)
	endif
 
	if(llist)then		!liste du modele sur l'ecran
	 do k=n,1,-1		
	  write(6,2001)1.d0-m(k)/mstar,p(k),t(k),r(k),l(k),compg(nbelem*(k-1)+1)
2001	  format((1x,1pd17.10,1p4d12.5,1pd10.3))
	 enddo
	endif
 
c	tabulation du TdS et interpolation en m23,
c	calcul des proportions en energies PP, CNO, 3alpha et gravifiques
 
c	print*,ipms
c	write(6,2000)dt,dt_
c	pause'tabulation TdS'
	e_tot=0.
	e_pp=0.
	e_cno=0.
	e_3a=0.
	e_gr=0.
	do i=1,n
	 if(i .eq. 1)then	!dm sert a integrer les e_pp etc..
	  dm=m(2)-m(1)		!au facteur mstar pres
	 elseif(i .eq. n)then
	  dm=m(n)-m(n-1)
	 else
	  dm=m(i+1)-m(i-1)
	 endif

c-----------------------------------------
c modtds
      if(i .eq. 1) then
        dlum=(l(2)-l(1))*lsol
      elseif(i .eq. n) then
        dlum=(l(n)-l(n-1))*lsol
      else
        dlum=(l(i+1)-l(i-1))*lsol
      endif
  
c      print*, 'm(n)-m(1) ', m(n)-m(1)
c      print*, 'mstar= ', mstar
c      print*, 'msol= ', msol
c      print*, 'l(n)= ', l(n)
c      pause

      dlum_s_dm=dlum/dm/mstar/msol
c----------------------------------------                                                                                                           

	 do j=1,nbelem
	  chimg(j)=compg(nbelem*(i-1)+j)
	  xchim(j)=compchim(nbelem*(i-1)+j)
	 enddo
	 call etat(p(i),t(i),chimg,.false.,	!en t+dt
     1	ro(i),drop,drot,drox,drott,drotp,drotx,
     2	eint(i),dup,dut,dux,dutt,dutp,dutx,nh1(i),
     3	nhe1(i),nhe2(i),bidt)
	 delta(i)=-t(i)/ro(i)*drot
	 cp(i)=dut+p(i)/ro(i)/t(i)*delta(i)
	
	 if(ipms .le. 2)then	!pour l'initialisation PMS
	  call iben_3(t(i),ro(i),xchim,dcomp,bidt,.false.,3,
     1	epsilo,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt)
	  tdst(i)=-epsilo(1)*secon6	!dans iben cT sort dans epsilo(1)
	  epsilo(1)=0	!mis a 0 pour le bilan des sources d'energie
	
	 else		!ipms>2
	  call nuc(t(i),ro(i),xchim,dcomp,bidt,.false.,3,
     1	epsilo,bidt,bidt,depsx,bidt,bidt,bidt,bidt,bidt,bidt)
          do j= 1,5
           epsilon(5*(i-1)+j)=epsilo(j)
          enddo
	  if(dt_ .gt. 0.)then	!1ier passage imps=3, dt_=0 --> ct dans TdS
	   call inter_3('mu',bp_t,q_t,qt_t,n_t,knot_t,m23(i),f,dfdx,r2_t,m23_t)	
	   p_t=exp(f(1))
	   t_t=exp(f(2))
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c------------------------------
c modtds
      tdstest=.true.
      if(tdstest) then
        tdst(i)=-(dlum_s_dm-(epsilon(5*(i-1)+2)+epsilon(5*(i-1)+3)
     +                   +epsilon(5*(i-1)+4)))
      else                                                                                      
	   if(kipp)then	!approximation de Kippenhahan                                                                          
	    tdst(i)=(cp(i)*(t(i)-t_t)-delta(i)/ro(i)*(p(i)-p_t))/dt_	                                                          
	   else                                                                                                                
	    call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,                                                                   
     1	knotc_t,.true.,min(m23(i),mc_t(nc_t)),lx,chimg,dchimg)                                                           
	    call chim_gram_3(chimg,dchimg,nuc_m)		!X, Y, Z                                                                     
c	    write(6,2000)m23(i),p_t,t_t	                                                                                      
	    call etat(p_t,t_t,chimg,.false.,                                                                                   
     1	ro_t,drop,drot,drox,drott,drotp,drotx,                                                                           
     2	u_t,dup,dut,dux,dutt,dutp,dutx,nh1(i),                                                                           
     3	nhe1(i),nhe2(i),bidt)                                                                                            
	    tdst(i)=(eint(i)-u_t-p(i)/ro(i)**2*(ro(i)-ro_t))/dt_                                                               
	   endif	!kipp     

        endif
c------------------------------------
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c	   write(6,2000)tdst(i),eint(i),u_t,p(i),p_t,t_t,ro(i),ro_t
           epsilon(5*(i-1)+5)=tdst(i)
           epsilon(5*(i-1)+1)=epsilon(5*(i-1)+1)+tdst(i)
	  else					!dt_ > 0
	   if(pms)then	!1ier passage apres modeles initiaux PMS (ipms=3)
	    tdst(i)=-c_iben*t(i)*secon6		!PMS
	    epsilo(1)=0.	!initialisation
	   else
	    tdst(i)=0		!modele de ZAMS
	   endif
	  endif		!dt_ > 0
	 endif 		!ipms < 2
	
c	 les energies et proportions en %
	
	 e_pp =e_pp +dm*epsilo(2)
	 e_cno=e_cno+dm*epsilo(3)
	 e_3a =e_3a +dm*epsilo(4)

c-------------
c modtds  
	 if ( tdstest ) then
            e_gr=e_gr+dm*abs(tdst(i))
         else
            e_gr =e_gr +dm*abs(tdst(i))/secon6
         end if 
c-------------

	enddo

	e_tot=e_pp+e_cno+e_3a+e_gr
 
	i_pp=nint(100.*e_pp/e_tot)
	i_cno=nint(100.*e_cno/e_tot)
	i_3a=nint(100.*e_3a/e_tot)
	i_gr=nint(100.*e_gr/e_tot)
c	print*,'energies',i_pp,i_cno,i_gr
c	write(6,2000)e_tot,e_pp,e_cno,e_3a,e_gr
c	pause'apres i_gr'
	 	
	call sbsp1dn(1,tdst,m23,m23t_t,n,m_ch,	!xchim, dxchim: VT
     1	knot23_t,.false.,m23(1),lx,xchim,dxchim)
c	do i=1,n
c	 call sbsp1dn(1,tdst,m23,m23t_t,n,m_ch,			!test
c	1	knot23_t,.true.,m23(i),lx,xchim,dxchim)
c	 print*,lx,m23(i),xchim(1)
c	enddo
c	pause'TdS'
 
c	write(6,2000)compg(nbelem*(n-1)+1),ab_min(1)
	pms_pr=pms
	zams_pr=zams
	post_pr=post
	cohe_pr=cohe
	if(compg(1) .gt. 1.d-3)then
	 if(i_pp+i_cno+i_3a .ge. 99)then
	  chaine='modele de la serie principale '
	  pms=.false.
	  zams=.true.
	  post=.false.
	  cohe=.false.
	 else
	  chaine='modele de pre-serie principale'
	  pms=.true.
	  zams=.false.
	  post=.false.
	  cohe=.false.
	 endif
	elseif(t(1) .lt. 1.d8)then
	 chaine= 'modele post serie principale  '
	 pms=.false.
	 zams=.false.
	 post=.true.
	 cohe=.false.
	else
	 chaine= 'modele avec combustion helium '
	 pms=.false.
	 zams=.false.
	 post=.false.
	 cohe=.true.
	endif
	
c	sauf pour le premier passage (alors passe=.false.)	
 
	zams_e= .not.zams_pr .and. zams .and. passe	!on passe PMS--->ZAMS
	post_e= .not.post_pr .and. post .and. passe	!on passe ZAMS-->POST
	cohe_e= .not.cohe_pr .and. cohe .and. passe	!on passe POST-->COHE
	passe=.true.
	
	if(zams_e)then
	 open(unit=26,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.zams')
	 modeleb='m'//cm//'X'//ch//'a'//cr//'_B.zams'
	 write(26)age,nbelem,nchim,mtot,alpha,modeleb,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	 close(unit=26)	!fermeture du fichier de ZAMS
	endif
	if(post_e)then
	 open(unit=26,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.post')
	 modeleb='m'//cm//'X'//ch//'a'//cr//'_B.post'
	 write(26)age,nbelem,nchim,mtot,alpha,modeleb,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	 close(unit=26)	!fermeture du fichier de POST
	endif
 
	if(cohe_e)then
	 open(unit=26,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.cohe')
	 modeleb='m'//cm//'X'//ch//'a'//cr//'_B.cohe'
	 write(26)age,nbelem,nchim,mtot,alpha,modeleb,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	 close(unit=26)	!fermeture du fichier de COHE
	endif
 
c	write(6,*)'e_tot,e_pp,e_cno,e_3a,e_gr',e_tot,e_pp,e_cno,e_3a,e_gr
c	pause
 
c---------------------------------------------------------------------
 
c		E C R I T U R E S
 
c---------------------------------------------------------------------
c-------------
c modifarret
	logg=log10(msol*mstar/rsol/rsol*g)	!pour ecriture de Log Teff
	
	lteff=log10(cte8*l(n)/r(n)**2)/4    !teff approche pour le test d'arret
	
	sort=		(age .ge. agemax-1.d-5)
     1	.or. (compg(1) .lt. x_stop)
     3	.or. (t(1) .ge. t_stop)
     4	.or. (zams .and. arret .eq. 'zams')
     5	.or. (post .and. arret .eq. 'post')
     6	.or. (cohe .and. arret .eq. 'cohe')
     7  .or. (compg(1) .le. 1.d-3 .AND. compg(2)+compg(3) .le. y_stop ) ! Modif. DC Nov. 98
     8  .or. (.NOT. lcroissant .AND. ( log10(l(n)) .le. log_lslsol) ) ! idem
	
	if(log_teff .gt. 0.)then
	 sort=sort .or. lteff .gt. log_teff
	else
	 sort=sort .or. lteff .lt. -log_teff
	endif	
	
	ecrit=sort .or. zams_e .or. post_e .or. cohe_e .or. (age-agep .ge. dtlist)
		
c	ecriture du *.osc, du *.dat et arret	
	
	if(ecrit)then
c	 write(48,*)'**n/age/m,p,t,ro,r,l,xchim,dxchim/epsi(5),xdot,ydot**'
c	 write(48,*)n
c	 write(48,2000)age
	 do i=1,n
	  do j=1,nbelem
	   xchim(j)=compchim(nbelem*(i-1)+j)
	   dxchim(j)=dcompchim(nbelem*(i-1)+j)
	   chimg(j)=compg(nbelem*(i-1)+j)
	   dchimg(j)=dcompg(nbelem*(i-1)+j)
	  enddo
c	  write(6,2000)m(i),r(i),t(i),p(i),l(i),(chimg(j),j=1,nbelem)
 
c	  elements lourds	
	
	  if(iz .gt. 0)then
	   z(i)=xchim(iz)
	  elseif(ihe4 .gt. 1 .and. .not.z_cte)then
	   z(i)=1.d0-chimg(1)-chimg(ihe4)-chimg(ihe4-1)
	  else
	   z(i)=z0	
	  endif
 
	  if(ipms .le. 2)then			!TdS par Iben
	   call thermo_3(p(i),t(i),xchim,m(i),l(i),r(i),.true.,dxchim,mstar,
     1	ro(i),drop,drot,drox,eint(i),dup,dut,dux,lambda(i),r_zc,r_ov,
     2	lim,gradient(i),dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilo,depsp,depst,depsx,kap(i),dkapp,dkapt,dkapx,
     4	delta(i),deltap,deltat,deltax,cp(i),dcpp,dcpt,dcpx,
     5	hp(i),gradad(i),gradrad(i),gradconv(i),d_grad(i),w(i),dw(i),
     6	nh1(i),nhe1(i),nhe2(i),etat,opa,conv,iben_3)
	   epsilo(1)=0
	  else
	   call thermo_3(p(i),t(i),xchim,m(i),l(i),r(i),.true.,dxchim,mstar,
     1	ro(i),drop,drot,drox,eint(i),dup,dut,dux,lambda(i),r_zc,r_ov,
     2	lim,gradient(i),dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilo,depsp,depst,depsx,kap(i),dkapp,dkapt,dkapx,
     4	delta(i),deltap,deltat,deltax,cp(i),dcpp,dcpt,dcpx,
     5	hp(i),gradad(i),gradrad(i),gradconv(i),d_grad(i),w(i),dw(i),
     6	nh1(i),nhe1(i),nhe2(i),etat,opa,conv,nuc)
	  endif		!ipms
	  alfa(i)=p(i)/ro(i)*drop	!un peu de thermo
	  beta(i)=1.d0-aradias3*t(i)**4/p(i)
	  dcapdr(i)=dkapp/drop
	  dcapdt(i)=dkapt-dcapdr(i)*drot
	  depsdr(i)=depsp/drop
	  depsdt(i)=depst-depsdr(i)*drot
 
c	  calcul d'une des formes de Vaissala:
c	  gm/r^2(1/gamma1 dlnP/dr - dln ro/dr)
c	  uniquement valable, pour les ZR, que dans les parties tot. ionisees	
	
	
	  if(hp(i) .le. 0.d0)then		!au centre
	   vaissala(i)=0.
	
	  else
 
	   vaissala(i)=r(i)*rsol/hp(i)*delta(i)*(gradad(i)-gradient(i))
     1	-cte3*r(i)**3*drox*dchimg(1)
	
c	   write(6,*)'m,r,drot,t,p,gradient,gradad,drox/vai,dchimg,nuc'
c	   write(6,2000)m(i),r(i),drot,t(i),p(i),gradient(i),gradad(i),drox
c	   write(6,2000)vaissala(i),dchimg(1)
 
	  endif
	  convec(i)=d_grad(i) .ge. 0.d0
 
c 	  production de neutrinos
 
	  call nuc(t(i),ro(i),xchim,bidt,bidt,.false.,4,
     1	bidt,bidt,bidt,bidt,anupp(i),anupep(i),anub8(i),anube7(i),
     2	anun13(i),anuo15(i))
 
c	  dX/dt et dY/dt
 
	  if(ipms .le. 2)then
	   call iben_3(t(i),ro(i),xchim,dcomp,bidt,.false.,2,
     1	bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt)
	  else
	   call nuc(t(i),ro(i),xchim,dcomp,bidt,.false.,2,
     1	bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt,bidt)
	  endif
	  xdot(i)=dcomp(1)/secon6*ah
	  if(ihe4 .le. 1)then
	   ydot(i)=0.
	  else
	   ydot(i)=dcomp(ihe4)/secon6*ahe4
	  endif
 
c	  les epsilon et dR/dt
 
	  call sbsp1dn(1,tdst,m23,m23t_t,n,m_ch,
     1	knot23_t,.true.,m23(i),lx,tds,bid)
	  epsilo(5)=tds/secon6
	  if(dt_ .le. 0.)then
	   rdot(i)=0.
	  else
	   call inter_3('mu',bp,q,qt,n,knot,m23(i),f,dfdx,r2,m23)	
	   rdot(i)=(sqrt(abs(bp(ne*(i-1)*m_qs+3)))-sqrt(f(3)))/dt_*cte7
	  endif
	  epsilo(1)=epsilo(1)-epsilo(5)
c	  write(6,*)'epsilo'
c	  write(6,2000)(epsilo(j),j=1,5)
	  do j=1,5
	   epsilon(5*(i-1)+j)=epsilo(j)
	  enddo
	 enddo
	
c	 atmosphere, temperature effective, rayon exterieur
 
c	 write(6,*)'l(n),r(n),xchim en n'
c	 write(6,2000)l(n),r(n),(chim1g(j),j=1,nchim)
	 call lim_ext(.true.,l(n),r(n),chim1g,pext,text,dpdl,dpdr,dtdl,dtdr,
     1	teff,rstar,mext,dml,dmr,p_atm,t_atm,m_atm,tau,r_atm,mstar,
     2	tdetau,etat,opa)
	 do i=1,nchim
	  dxchim(i)=0
	 enddo
	 do i=1,n_atm		!thermo pour l'atmosphere
	  lnp_a(i)=log(p_atm(i))
	  lnt_a(i)=log(t_atm(i))
c	  write(6,2000)p_atm(i),t_atm(i),xchim1(1),m_atm(i),l(n),r(n)
	  if(ipms .le. 2)then		!les reac_nuc par iben
	   call thermo_3(p_atm(i),t_atm(i),xchim1,m_atm(i),l(n),r(n),.false.,
     1	dxchim,mstar,ro_atm(i),drop,drot,drox,u,dup,dut,dux,lamb,r_zc,
     2	r_ov,lim,grad_atm(i),dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsa,depsp,depst,depsx,k_atm(i),dkapp,dkapt,dkapx,
     4	delta_atm(i),deltap,deltat,deltax,cp_atm(i),dcpp,dcpt,dcpx,
     5	hpa,grada_atm(i),gradr_atm(i),gradc_atm(i),bidt,w(n),0.d0,
     6	nh1_atm(i),nhe1_atm(i),nhe2_atm(i),etat,opa,conv,iben_3)
	   epsa(1)=0
	  else
	   call thermo_3(p_atm(i),t_atm(i),xchim1,m_atm(i),l(n),r(n),.false.,
     1	dxchim,mstar,ro_atm(i),drop,drot,drox,u,dup,dut,dux,lamb,r_zc,
     2	r_ov,lim,grad_atm(i),dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsa,depsp,depst,depsx,k_atm(i),dkapp,dkapt,dkapx,
     4	delta_atm(i),deltap,deltat,deltax,cp_atm(i),dcpp,dcpt,dcpx,
     5	hpa,grada_atm(i),gradr_atm(i),gradc_atm(i),bidt,w(n),0.d0,
     6	nh1_atm(i),nhe1_atm(i),nhe2_atm(i),etat,opa,conv,nuc)
	  endif
	  ldcapdr_a(i)=dkapp/drop
	  ldcapdt_a(i)=dkapt-ldcapdr_a(i)*drot
c	  write(6,2000)m_atm(i),r_atm(i)
	  alfa_atm(i)=p_atm(i)/ro_atm(i)*drop
	
c	  dans l'atmosphere pas de gradient de comp. chim	
	 	
	  vais_atm(i)=r_atm(i)*rsol/hpa*delta_atm(i)*(grada_atm(i)-grad_atm(i))
 
c	  write(6,*)'m,r,drot,t,p,gradient,gradad,vai'
c	  write(6,2000)m_atm(i),r_atm(i),drot,t_atm(i),p_atm(i),grad_atm(i),
c	1	grada_atm(i),vais_atm(i)
	 enddo		!le gradient pour Marie Jo
	 if(n_atm .gt. m_ch)then
	  do i=1,n_atm	!on prend l'oppose pour tabuler ln T en fct. de ln P
	   lnp_a(i)=-lnp_a(i)	!les abs. devant etre croissantes
	  enddo
	  call sbsp1d(lnt_a,lnp_a,bidt,n_atm,m_ch,
     1	knota,.false.,lnp_a(1),lx,d2p,grad_mj_a(1))	!d2p: VT
	  do i=1,n_atm
	   call sbsp1d(lnt_a,lnp_a,bidt,n_atm,m_ch,knota,.true.,
     1	lnp_a(i),lx,d2p,grad_mj_a(i))		!d2p: VT
	  enddo
	  do i=1,n_atm
	   lnp_a(i)=-lnp_a(i)
	   grad_mj_a(i)=-grad_mj_a(i)	!on prend l'oppose
c	   write(6,2000)p_atm(i),grad_mj_a(i),grad_atm(i)	
	  enddo
c	  pause'atm cesam_3'
	 endif
	
c	 sorties sur ecran
 
c	 print*,mstar,mtot
	 pertem=mstar/mtot-1.d0
	 if(ihe4 .le. 1)then
	  write(6,205)age,log10(teff),log10(l(n)),log10(r(n)),
     1	logg+log10(m(n)/r(n)**2),
     2	p(1),t(1),ro(1),compg(nbelem*(1-1)+1),
     3	i_pp,i_cno,i_3a,i_gr,pertem,mstar,chaine
205	  format(1x,'*********',/,1x,'age=',1pd10.3,' log Teff=',1pd10.3,
     1	' log L/Lsol=',1pd10.3,' log R/Rsol=',1pd10.3,/,
     2	1x,'Log10 g=',1pd10.3,' Pc=',1pd10.3,' Tc=',1pd10.3,
     3	' Roc=',1pd10.3,' Xc=',1pd10.3,
     4	/,1x,'en. PP=',i3,'%, en. CNO=',i3,'%, en. 3 alpha=',i3,
     5	'%, en. grav=',i3,'%',/,1x,'Var. rel. de masse :',1pd10.3,
     6	', M*=',1pd10.3,'Msol, ',a30,/,1x,'*********')
	 else
	  write(6,206)age,log10(teff),log10(l(n)),log10(r(n)),
     1	logg+log10(m(n)/r(n)**2),
     2	p(1),t(1),ro(1),compg(nbelem*(1-1)+1),
     3	i_pp,i_cno,i_3a,i_gr,compg(nbelem*(1-1)+ihe4),
     4	pertem,mstar,chaine
206	  format(1x,'*********',/,1x,'age=',1pd10.3,' log Teff=',1pd10.3,
     1	' log L/Lsol=',1pd10.3,' log R/Rsol=',1pd10.3,/,
     2	1x,'Log10 g=',1pd10.3,' Pc=',1pd10.3,' Tc=',1pd10.3,
     3	' Roc=',1pd10.3,' Xc=',1pd10.3,
     4	/,1x,'en. PP=',i3,'%, en. CNO=',i3,'%, en. 3 alpha=',i3,
     5	'%, en. grav=',i3,'%',' Yc=',1pd10.3,/,1x,
     6	'Var. rel. de masse :',1pd10.3,
     7	', M*=',1pd10.3,'Msol, ',a30,/,1x,'*********')	
	 endif
	
c	 print*,dessin
c	 pause'avant dessin'
 
	 if(dessin)call des_cesam(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,teff,
     1	r2,m23,age,lim,lconv,m_zc,mstar)
	
c***************************************************
clehr
c Determination des zones de combustion nucleaire
c -----------------------------------------------
 
c lpp : variable logique, lpp=.true. si on se trouve dans une region ou
c       \epsilon_{pp} .ge. 10.**3 erg.s-1.g-1
c lcno, l3a, lgrav : idem avec le cycle CNO, la reaction 3 alpha et energie
c gravitationnelle (+ ou -).
 
c n_pp     : nbre de limites non-combustion/combustion pour PP.
c n_cno    : idem pour le cycle CNO.
c n_3a     : idem pour les chaines 3 alpha.
c n_grav_p : idem pour l'energie gravitationnelle positive.
c n_grav_m : idem pour l'energie gravitationnelle negative.
 
c natpp(10) : tableau de variables logiques, exemple : natpp(3)=.true.
c             la troisieme limite non_combu./combu. est le debut d'une
c             zone de combu. quand on va de la surface vers le centre
c             rq.:valeurs max. : 10 !!
 
c pos_lpp(10)   : position en m/M* des limites des zones de combustion (PP)
c pos_cno(10)   : idem pour CNO
c pos_3a(10)    : idem pour 3 alpha
c pos_grav_p(10): idem energie grav. +
c pos_grav_m(10): idem energie grav. -
 
c t_lpp(10)    : temperature des limites des zones PP
c t_lcno(10)   : idem pour CNO
c t_l3a(10)    : idem pour 3 alpha
c t_lgrav_p(10): idem energie grav. +
c t_lgrav_m(10): idem energie grav. -
 
c Initialisation des differentes quantites :
 
      lpp=.false.
      lcno=.false.
      l3a=.false.
      lgrav_p=.false.
      lgrav_m=.false.
 
      n_pp=0
      n_cno=0
      n_3a=0
      n_grav_p=0
      n_grav_m=0
 
c Determination des zones Epsilon_pp > 10.**3
 
      do indi= n, 1, -1
 
         if( epsilon(5*(indi-1)+2) .lt. 1.d+3) then
           if(lpp)then
              lpp=.false.
              n_pp=n_pp+1
              natpp(n_pp)=lpp
              pos_lpp(n_pp)=1.d+0-m(indi)/mstar
              t_lpp(n_pp)=log10(t(indi))
              num_pp(n_pp)=indi
           end if
         else
           if(.not.lpp)then
              lpp=.true.
              n_pp=n_pp+1
              natpp(n_pp)=lpp
              pos_lpp(n_pp)=1.d+0-m(indi)/mstar
              t_lpp(n_pp)=log10(t(indi))
              num_pp(n_pp)=indi
           end if
         end if
 
c Determination des zones Epsilon_cno > 10.**3
 
         if( epsilon(5*(indi-1)+3) .lt. 1.d+3) then
           if(lcno)then
              lcno=.false.
              n_cno=n_cno+1
              natpp(n_cno)=lcno
              pos_lcno(n_cno)=1.d+0-m(indi)/mstar
              t_lcno(n_cno)=log10(t(indi))

              rho_lcno(n_cno)=log10(ro(indi))
              x_lcno(n_cno)=  compchim(nbelem*(indi-1)+1)*nucleo(1) ! H
              y_lcno(n_cno)=  compchim(nbelem*(indi-1)+2)*nucleo(2) ! He3
     +                      + compchim(nbelem*(indi-1)+3)*nucleo(3) ! He4
              Xcno_lcno(n_cno)=  compchim(nbelem*(indi-1)+4)*nucleo(4) ! C12
     +                      +    compchim(nbelem*(indi-1)+5)*nucleo(5) ! C13
     +                      +    compchim(nbelem*(indi-1)+6)*nucleo(6) ! N13
     +                      +    compchim(nbelem*(indi-1)+7)*nucleo(7) ! N14
     +                      +    compchim(nbelem*(indi-1)+8)*nucleo(8) ! N15
     +                      +    compchim(nbelem*(indi-1)+9)*nucleo(9) ! O16
     +                      +    compchim(nbelem*(indi-1)+10)*nucleo(10) ! O17
     +                      +    compchim(nbelem*(indi-1)+11)*nucleo(11) ! O18


              num_cno(n_cno)=indi
           end if
         else
           if(.not.lcno)then
              lcno=.true.
              n_cno=n_cno+1
              natcno(n_cno)=lcno
              pos_lcno(n_cno)=1.d+0-m(indi)/mstar
              t_lcno(n_cno)=log10(t(indi))

              rho_lcno(n_cno)=log10(ro(indi))
              x_lcno(n_cno)=  compchim(nbelem*(indi-1)+1)*nucleo(1) ! H
              y_lcno(n_cno)=  compchim(nbelem*(indi-1)+2)*nucleo(2) ! He3
     +                      + compchim(nbelem*(indi-1)+3)*nucleo(3) ! He4
     +                      +    compchim(nbelem*(indi-1)+5)*nucleo(5) ! C13
     +                      +    compchim(nbelem*(indi-1)+6)*nucleo(6) ! N13
     +                      +    compchim(nbelem*(indi-1)+7)*nucleo(7) ! N14
     +                      +    compchim(nbelem*(indi-1)+8)*nucleo(8) ! N15
     +                      +    compchim(nbelem*(indi-1)+9)*nucleo(9) ! O16
     +                      +    compchim(nbelem*(indi-1)+10)*nucleo(10) ! O17
     +                      +    compchim(nbelem*(indi-1)+11)*nucleo(11) ! O18

              num_cno(n_cno)=indi
           end if
         end if
 
c Determination des zones Epsilon_3a > 10.**3
 
         if( epsilon(5*(indi-1)+4) .lt. 1.d+3) then
           if(l3a)then
              l3a=.false.
              n_3a=n_3a+1
              nat3a(n_3a)=l3a
              pos_l3a(n_3a)=1.d+0-m(indi)/mstar
              t_l3a(n_3a)=log10(t(indi))
              num_3a(n_3a)=indi
           end if
         else
           if(.not.l3a)then
              l3a=.true.
              n_3a=n_3a+1
              nat3a(n_3a)=l3a
              pos_l3a(n_3a)=1.d+0-m(indi)/mstar
              t_l3a(n_3a)=log10(t(indi))
              num_3a(n_3a)=indi
           end if
         end if
 
c Determination des zones Epsilon_grav > + 10.**2
 
         if( epsilon(5*(indi-1)+5) .lt. 1.d+2) then
           if(lgrav_p)then
              lgrav_p=.false.
              n_grav_p=n_grav_p+1
              natgrav_p(n_grav_p)=lgrav_p
              pos_lgrav_p(n_grav_p)=1.d+0-m(indi)/mstar
              t_lgrav_p(n_grav_p)=log10(t(indi))
              num_grav_p(n_grav_p)=indi
           end if
         else
           if(.not.lgrav_p)then
              lgrav_p=.true.
              n_grav_p=n_grav_p+1
              natgrav_p(n_grav_p)=lgrav_p
              pos_lgrav_p(n_grav_p)=1.d+0-m(indi)/mstar
              t_lgrav_p(n_grav_p)=log10(t(indi))
              num_grav_p(n_grav_p)=indi
           end if
         end if
 
c Determination des zones Epsilon_grav < - 10.**2
 
         if( epsilon(5*(indi-1)+5) .gt. -1.d+2) then
           if(lgrav_m)then
              lgrav_m=.false.
              n_grav_m=n_grav_m+1
              natgrav_m(n_grav_m)=lgrav_m
              pos_lgrav_m(n_grav_m)=1.d+0-m(indi)/mstar
              t_lgrav_m(n_grav_m)=log10(t(indi))
              num_grav_m(n_grav_m)=indi
           end if
         else
           if(.not.lgrav_m)then
              lgrav_m=.true.
              n_grav_m=n_grav_m+1
              natgrav_m(n_grav_m)=lgrav_m
              pos_lgrav_m(n_grav_m)=1.d+0-m(indi)/mstar
              t_lgrav_m(n_grav_m)=log10(t(indi))
              num_grav_m(n_grav_m)=indi
           end if
         end if
 
      end do
 
c   Fin de la determination des zones "energetiques"
c------------------------------------------------------------
c Caracterisation du coeur d'helium, determination de Rc, Mc,
c Lc et Phic

c Rq.: 'Phi_c' n'est pas exactement le potentiel du coeur d'helium
c      mais une quantite qui s'en approche.
      log_lum_c=0.d0
      masse_c=0.d0
      rayon_c=0.d0
      log_phi_c=0.d0	
      l_c=.false.

      xsurf=compchim(nbelem*(n-1)+1)*nucleo(1) ! X a la surface.

      do indi= n, 1, -1
 
      if( compchim(nbelem*(indi-1)+1)*nucleo(1) .le. 0.5d0*xsurf 
     +    .AND. (.NOT. l_c) )then 
        l_c=.true.
        rayon_c=r(indi) ! r est en rsol
        log_lum_c=log10(l(indi))
        masse_c=m(indi) ! m est en msol
        log_phi_c=log10(masse_c/rayon_c)
      end if
 
      end do
 
c------------------------------------------------------------
c Determination des pourcentages de production d'energie
c
c      ener_pp=0.d0
c      ener_cno=0.d0
c      ener_3a=0.d0
c      ener_grav=0.d0
c
c      do indi= 1, n-1
c
c         dmasse=m(indi+1)-m(indi)
c
c         ener_pp=ener_pp+epsilon(5*(indi-1)+2)*dmasse
c         ener_cno=ener_cno+epsilon(5*(indi-1)+3)*dmasse
c         ener_3a=ener_3a+epsilon(5*(indi-1)+4)*dmasse
c         ener_grav=ener_grav+epsilon(5*(indi-1)+5)*dmasse
c
c      end do
c
c      ener_tot=ener_pp+ener_cno+ener_3a+ener_grav
c
c      i_pp=nint(100.*ener_pp/ener_tot)
c      i_cno=nint(100.*ener_cno/ener_tot)
c      i_3a=nint(100.*ener_3a/ener_tot)
c      i_gr=nint(100.*ener_grav/ener_tot)
c
c------------------------------------------------------------------------
c Determination du facteur de degenerescence Lambda au centre
 
      do indi= 1, nbelem  ! Composition chimique au centre
         l_xchi(indi)=compchim(nbelem*(1-1)+indi)*nucleo(indi)
      end do
 
      call etat(p(1),t(1),l_xchi,.false.,
     +          bid1,bid2,bid3,bid4,bid5,bid6,bid7,bid8,bid9,
     +          bid10,bid11,bid12,bid13,bid14,bid15,bid16,bid17,
     +          lambda_c)
 
c-----------------------------------------------------------------------
c Determination de la zone ou Z > 0.1 :
 
      z_lim=0.1d0
 
      do indi= n, 1, -1
         z_local=1.d0-compchim(nbelem*(indi-1)+1)*nucleo(1)
     +           -compchim(nbelem*(indi-1)+2)*nucleo(2)
     +           -compchim(nbelem*(indi-1)+3)*nucleo(3)
         if( z_local .ge. z_lim ) then
           change_compo=n
           goto 1111
         end if
      end do
1111   m_change_compo=1.d0-m(change_compo)/mstar
 
c-----------------------------------------------------------------------
c      Les Ecritures
 
c        10        20        30        40        50        60        70
      write(53,2201)'--------------------------------------------------'
 
      write(53,2101) age, n ! age et nombre de couches
c2101  format(1pd22.15,i6)
 
      write(53,2102) mstar, log10(r(n)), log10(l(n)), log10(teff)
c2102  format(1p4d13.6)
 
      write(53,2112) i_pp, i_cno, i_3a, i_gr ! les % des diff. sources d'energie
c2112  format(1p,4i3)
 
      write(53,2103) log10(t(1)), log10(p(1)), log10(ro(1)),
     +          lambda_c
 
      write(53,2133) (compchim(nbelem*(1-1)+i)*nucleo(i),i=1,18)
c2103  format(1p4d13.6)
c2133  format(1p18d13.6)
 
c ---- Composition en surface ---------------------
 
      write(53,2113) (compchim(nbelem*(n-1)+i)*nucleo(i),i=1,18)
c2113  format(1p18d13.6)
 
c ---- Caracteristiques du coeur d'helium ------------------------------
 
      write(53,2124) rayon_c, masse_c, log_lum_c, log_phi_c, 'Coeur He'
c2124  format(1p4d13.6,A9)
 
c ---- Profondeur a laquelle il y changement de composition chimique ---
c (par rapport a la compo initiale du modele homogene!!)
 
      write(53,2123) m_change_compo
c2123  format(1pd13.6)
 
c ---- Zones de combustion PP ---------------------
 
      write(53,2104) n_pp, (natpp(i),i=1,n_pp)
c2104  format(1p,i3,10l2)
      do i= 1, n_pp
         write(53,2115) pos_lpp(i), t_lpp(i), num_pp(i)
      end do
c2105  format(1p2D13.6)
 
c ---- Zones de combustion CNO -------------------
 
      write(53,2104) n_cno, (natcno(i),i=1,n_cno)
      do i= 1, n_cno
         write(53,2115) pos_lcno(i), t_lcno(i), num_cno(i), rho_lcno(i), 
     +                  x_lcno(i), y_lcno(i), Xcno_lcno(i)
c2115   format(1p2d13.6,i5,4d16.6)
      end do
 
c ---- Zones de combustion 3 alpha -------------------
 
      write(53,2104) n_3a, (nat3a(i),i=1,n_3a)
      do i= 1, n_3a
         write(53,2115) pos_l3a(i), t_l3a(i), num_3a(i)
      end do
 
c ---- Zones d'energie gravita. + -------------------
 
      write(53,2104) n_grav_p, (natgrav_p(i),i=1,n_grav_p)
      do i= 1, n_grav_p
         write(53,2115) pos_lgrav_p(i), t_lgrav_p(i), num_grav_p(i)
      end do
 
c ---- Zones d'energie gravita. - -------------------
 
      write(53,2104) n_grav_m, (natgrav_m(i),i=1,n_grav_m)
      do i= 1, n_grav_m
         write(53,2115) pos_lgrav_m(i), t_lgrav_m(i), num_grav_m(i)
      end do
 
c ---- Zones Convectives --------------------------
 
      write(53,2106) lim, (lconv(i),i=1,lim)
c2106  format(i3,20l2)
      do i= 1, lim
         write(53,2107) 1.d0-m_zc(i), r_zc(i), r_ov(i), logtzc(i), 
     +                  logtov(i), logrozc(i)
c2107  format(1p5d13.6,d16.6)
      end do

c****************************************************
c Ecriture dans le fichier pour Roxburgh
c Roxy
	if ( ROX ) then
	write(76,1076) log10(teff), diff, mroxSmsch, 
     +      compchim(nbelem*(1-1)+1)*nucleo(1), age
 1076	format(5d16.6)
	end if
c****************************************************
c	 derivees secondes au centre
 
	 d2p=cte5*ro(1)**2*rstar**2/p(1)
	 fonc(0)=ro(1)
	 fonc(1)=ro(2)
	 fonc(2)=ro(3)
	 fonc(3)=ro(4)
	 fonc(4)=ro(5)
	 fonc(5)=ro(6)
	 absc(0)=r(1)
	 absc(1)=r(2)
	 absc(2)=r(3)
	 absc(3)=r(4)
	 absc(4)=r(5)
	 absc(5)=r(6)
 
c	 write(6,*)'ro/r/deriv'
c	 write(6,2000)(fonc(i),i=0,5)
c	 write(6,2000)(absc(i),i=0,5)
 
	 call newton(0,m_qs+1,fonc,absc,0.d0,poly,2)
	 d2ro=rstar**2/ro(1)*poly(2)
 
c	 write(6,2000)(poly(i),i=0,2)
c	 write(6,*)'rstar,d2ro,ro,poly(2)'
c	 write(6,2000)rstar,d2ro,ro(n),poly(2)
c	 write(6,*)'d2p,d2ro,rstar,p(n),ro(n)'
c	 write(6,2000)d2p,d2ro,rstar,p(n),ro(n)
 
c	 verification de la solution: estimation des premiers et 2d membres
c	 des equations d'evolution
 
c 	 do i=2,n-1
c	  write(6,2000)(t(i+1)-t(i-1))/(m(i+1)-m(i-1))/(mtot*msol),
c	1	-g*m(i)*mtot*msol/4./pi/(rsol*r(i))**4*gradient(i)*t(i)/p(i),
c	2	(p(i+1)-p(i-1))/(m(i+1)-m(i-1))/(mtot*msol),
c	3	-g*m(i)*mtot*msol/4./pi/(rsol*r(i))**4,
c	4	(r(i+1)-r(i-1))*rsol/(m(i+1)-m(i-1))/(mtot*msol),
c	5	1./4./pi/(rsol*r(i))**2/ro(i),
c	6	(l(i+1)-l(i-1))*lsol/(m(i+1)-m(i-1))/(mtot*msol),
c	7	epsilon(5*(i-1)+1)
c	 enddo
c	 if(.true.)stop
 
c	 sortie
 
	 if(sort)then
	  call list_3(n,m,p,ro,t,r,l,compg,age,teff,modele,rstar,hp,mstar,
     1	methode,precix,.true.,kap,epsilon,eint,xdot,ydot,rdot,
     2	nh1,nhe1,nhe2,alfa,beta,delta,cp,vaissala,convec,dlpdxt,dlpdxr,
     3	gradconv,gradrad,gradad,dcapdt,dcapdr,depsdr,depsdt,d2p,d2ro,
     4	anupp,anupep,anub8,anube7,anun13,anuo15,lambda,w,z,chaine,
     5	p_atm,t_atm,m_atm,tau,r_atm,ro_atm,k_atm,gradr_atm,grada_atm,
     6	gradc_atm,n_atm,i_pp,i_cno,i_3a,i_gr)
	  write(6,*)' '
	  write(6,*)'nom du modele : ',modele
	  write(6,*)' '
	  write(6,*)'nom du fichier des NAMELISTs des donnees : ',
     1	nom_fich2(:long(nom_fich2))//'.don3'
	  write(6,*)' '
	  write(6,*)'nom du fichier du listing des resultats en ASCII : ',
     1	nom_fich2(:long(nom_fich2))//'.lis'
	  write(6,*)' '
	  write(6,*)'nom du fichier du modele d''age 0 PMS en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.pms'
	  write(6,*)' '
	  write(6,*)'nom du fichier du modele d''age 0 homogene en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.hom'
	  write(6,*)' '
	  write(6,*)'nom du fichier du modele de ZAMS en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.zams'
	  write(6,*)' '
	  write(6,*)'nom du fichier du modele de POST en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.post'
	  write(6,*)' '
	  write(6,*)'nom du fichier du modele de COHE en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.cohe'
	  write(6,*)' '
	  write(6,*)'nom du fichier du modele final en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.dat'
	  write(6,*)' '
	  write(6,*)'nom des fichiers des modeles intermediaires en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.rep'
	  write(6,*)' '
	  write(6,*)'nom du fichier pour trace du diag. HR et ZC : ',
     1	nom_fich2(:long(nom_fich2))//'.HR'
	  if(n_atm .gt. 0)then
	   write(6,*)' '
	   write(6,*)'nom du fichier du modele d''atmosphere en binaire : ',
     1	nom_fich2(:long(nom_fich2))//'_B.atm'
	  endif
	  write(6,*)' '
	  close(unit=2)
 
c	  on garde le modele final
 
	  open(unit=25,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.dat')
	  modeleb='m'//cm//'X'//ch//'a'//cr//'_B.dat'
	  write(25)age,nbelem,nchim,mtot,alpha,modeleb,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	
	  close(unit=25)		!fichier *.dat
	  close(unit=53)		!fichier du diagramme HR et des ZC
c Roxy
          if ( ROX ) then
	     close(unit=76)	!fichier pour Roxburgh
	  else
	     close(unit=76,status='delete') 
          end if
 
c	  fichier d'oscillation
 
	  write(6,*)' '
	  write(6,*)'formation du fichier pour oscillations ? o/n'
	  !read(5,'(a)')oui
	  oui='o' ! Modif D.C.
	  if(oui .eq. 'o')then
	   open(unit=30,form='formatted',status='unknown',
     &	file=nom_fich2(:long(nom_fich2))//'.osc')
	   write(6,*)'nom du fichier ASCII pour le calcul des oscillations :',
     &	nom_fich2(:long(nom_fich2))//'.osc'
 
	   open(unit=32,form='formatted',status='unknown',
     &	file=nom_fich2(:long(nom_fich2))//'.struc')
	   write(6,*)' '
	   write(6,*)'nom du fichier ASCII contenant la strucutre de l''étoile :',
     &	nom_fich2(:long(nom_fich2))//'.struc'
 
c------------- date sur SUN, ULTRIX -----------------------
 
c	   call idate(i,j,k)
c	   write(6,*)'idate',i,j,k
c	   pause'idate'
c	   encode(2,2,ci)i
           write(ci,2) (i)
c	   encode(2,2,cj)j
           write(cj,2) (j)
c	   encode(4,3,ck)k
           write(ck,3) (k)
 
c	   titre=modeleb(:long(modeleb))//' du '//ci//'/'//cj//'/'//ck
c	1	//' - P. Morel, Observatoire de Nice'
	
c------------- date sur HP ---------------------------
 
c	   call date(today)
c	   write(6,*)today
c	   pause'date'
 
	   titre=modeleb(:long(modeleb))//' du '//today
     1	//' - P. Morel, Observatoire de Nice'
	
c------------------------------------------------------	
	
	   write(30,*)titre(:long(titre))
	   write(30,*)'fichier pour le calcul des oscillations'
	   write(30,*)' fichier: '//nom_fich2(:long(nom_fich2))//'.osc'
	   write(30,89)methode
89	   format(1x,'methode: CESAM ',a60)
	   write(30,90)nbelem,(nom_elem(i),i=1,nbelem)
90	   format(i3,14(1x,a3))
 
	   glob(1)=mstar*msol
	   glob(2)=rstar*rsol
	   glob(3)=l(n)*lsol
	   glob(4)=z0
	   glob(5)=x0
	   glob(6)=alpha
	   glob(7)=9./4.
	   glob(8)=1./162.
	   glob(9)=compg(nbelem*(n-1)+1)		!X dans ZC
	   glob(10)=compg(nbelem*(n-1)+ihe4)	!Y dans ZC
	   glob(11)=d2p
	   glob(12)=d2ro
	   glob(13)=age
	   glob(14)=vsal
	   glob(15)=w_rot
	
	   nadd=max(0,n_atm-1)
	   if(nadd .gt. 0)then		!on ajoute l'atmosphere sauf la couche 1
	    if(ihe4 .gt. 1)then
	     y=chim1g(ihe4)
	    else	
	     y=1.d0-z0-chim1g(1)
	    endif
	    do i=1,nadd
	     ipn=n_atm-i+1
	     var(1,i)=r_atm(ipn)*rsol
	     var(2,i)=log(m_atm(ipn)/mstar)
	     var(3,i)=t_atm(ipn)
	     var(4,i)=p_atm(ipn)
	     var(5,i)=ro_atm(ipn)
	     var(6,i)=grad_mj_a(ipn)
c	     var(6,i)=grad_atm(ipn)
	     var(7,i)=l(n)*lsol
	     var(8,i)=k_atm(ipn)
	     var(9,i)=0.
	     var(10,i)=1.d0/(alfa_atm(ipn)-delta_atm(ipn)*grada_atm(ipn))
	     var(11,i)=grada_atm(ipn)
	     var(12,i)=delta_atm(ipn)
	     var(13,i)=cp_atm(ipn)
	     var(14,i)=nh1_atm(ipn)*xchim1(1)+(nhe1_atm(ipn)+2.*
     1	nhe2_atm(ipn))*y+z(n)*.5
	     var(15,i)=vais_atm(ipn)
	     var(16,i)=w(n)
	     var(17,i)=ldcapdt_a(ipn)*t_atm(ipn)/k_atm(ipn)
	     var(18,i)=ldcapdr_a(ipn)*ro_atm(ipn)/k_atm(ipn)
	     var(19,i)=0.
	     var(20,i)=0.
	     do j=1,nbelem
	      var(20+j,i)=chim1g(j)
	     enddo
	    enddo
	   endif		!atmosphere
 
c	   la structure interne
 
	   do i=n,1,-1	
	    do j=1,nbelem		!par gramme
	     chimg(j)=compg(nbelem*(i-1)+j)
	    enddo
	    if(ihe4 .gt. 1)then
	     y=chimg(ihe4)
	    else
	     y=1.d0-z0-chimg(1)
	    endif
	    	
	    ipn=n-i+nadd+1
	    var(1,ipn)=r(i)*rsol
	    if(m(i) .ne. 0.)then
	     var(2,ipn)=log(m(i)/mstar)
	    else
	     var(2,ipn)=-1.d38
	    endif
	    var(3,ipn)=t(i)
	    var(4,ipn)=p(i)
	    var(5,ipn)=ro(i)
	    var(6,ipn)=grad_mj(i)
c	    var(6,ipn)=gradient(i)	
	    var(7,ipn)=l(i)*lsol
	    var(8,ipn)=kap(i)
	    var(9,ipn)=epsilon(5*(i-1)+1)
	    var(10,ipn)=1.d0/(alfa(i)-delta(i)*gradad(i))
	    var(11,ipn)=gradad(i)
	    var(12,ipn)=delta(i)
	    var(13,ipn)=cp(i)
	    if(max(nh1(i),nhe1(i),nhe2(i)) .lt. 0.)then	!taux d'ionis. inconnus
	     var(14,ipn)=0
	    else
	     var(14,ipn)=nh1(i)*chimg(1)+(nhe1(i)+2.*
     1	nhe2(i))*y+z(i)*.5	
	    endif
	    if(m(i)*r(i) .ne. 0.)then
	     var(15,ipn)=vaissala(i)
	    else
	     var(15,ipn)=0.
	    endif
	    var(16,ipn)=w(i)		!rotation au lieu de xdot
	    var(17,ipn)=dcapdt(i)*t(i)/kap(i)
	    var(18,ipn)=dcapdr(i)*ro(i)/kap(i)
	    var(19,ipn)=depsdt(i)*t(i)
	    var(20,ipn)=depsdr(i)*ro(i)
	    do j=1,nbelem
	     var(20+j,ipn)=compg(nbelem*(i-1)+j)
	    enddo
	   enddo
	   iconst=15		!nb de global
	   ivar=20		!nb de variables
	   ivers=3		!numero de version
	   write(30,137)n+nadd,iconst,ivar,ivers,nchim,iw,iz
137	   format(7i10)
	   write(30,138)(glob(i),i=1,iconst)
138	   format(1p5d19.12)
	   do j=1,n+nadd
	    write(30,138)(var(i,j),i=1,ivar+nbelem)
	   enddo
	   close(unit=30)
	   ! Ecriture dans le fichier de structure : --------------------------------------
	   write(32,'(A1)') '#'
	   write(32,'(A10,F10.4)') '# Age   = ', age
	   write(32,'(A10,F10.4)') '# Mstar = ', mstar
	   write(32,'(A10,F10.4)') '# Rstar = ', rstar
	   write(32,'(A10,F10.4,A13,F10.4)') '# Lstar = ', l(n), '  Log L    = ', log10(l(n))
	   write(32,'(A10,F10.4,A13,F10.4)') '# Tstar = ', t(n), '  Log Teff = ', log10(t(n))
	   write(32,'(A10,F10.4)') '# Alpha = ', alpha
           write(32,'(A10,F10.4)') '# Ov.sup= ', OVSHTS
           write(32,'(A10,F10.4)') '# Ov.inf= ', OVSHTI
	   write(32,'(A10,F10.4)') '# X0    = ',      x0
           write(32,'(A10,F10.4)') '# Y0    = ',      y0
	   z0=1.d0-x0-y0
	   write(32,'(A10,F10.4)') '# Z0    = ',      z0
	   write(32,'(A1)') '#'
	   write(32,'(A1,A4,10A14)') '#', '    ', 'm/Mstar','logP','logT','logRo',
     &              'r/rsol', 'logL/Lsol', 'Delad', 'Delrad','X', 'Y'
	   write(32,'(A1)') '#'
	   do i= n, 1, -1
	       delrad= 3./16./pi/aradia/clight/g*var(8,i)*var(7,i)*var(4,i)
     &               /10.**var(2,i)/Mstar/msol/var(3,i)**4
	       xxx = var(21,i)
	       yyy = var(23,i)
	       if ( xxx .lt. 1.d-5 ) then
		  xxx = 0.d0
	       end if
	       if ( yyy .lt. 1.d-5 ) then
		  yyy = 0.d0
	       end if
              !                             m/Mstar       , logP           , logT           , logRo
	      write(32,'(I5,10d14.5)') n-i+1,10.**var(2,i),log10(var(4,i)),log10(var(3,i)),log10(var(5,i)),
     &            var(1,i)/rsol, log10(var(7,i)/lsol), var(11,i), delrad
     &            , xxx, yyy
	    end do ! Fin de l'écriture du fichier de structure
	   ! -----------------------------------------------------------------------------
	   close(unit=32)
	  endif
c	  if(dessin)call pgend
 
	  write(6,*)' '
	  write(6,*)'***********************************'
	  write(6,*)' '
	  write(6,*)' Sortie de CESAM version: ',version
	  write(6,*)' '
	  write(6,*)'Ce code d''evolution stellaire a ete elabore',
     1	'dans le cadre du Groupement de Recherche Structure Interne',
     2	'des Etoiles et des Planetes Geantes.'
	  write(6,*)' '
	  write(6,*)'  Si son utilisation vous a donne satisfaction, le but ',
     1	'poursuivi par tous ceux qui ont y ont contribue ',
     2	'aura ete atteint.'
	  write(6,*)' '
	  write(6,*)'P. Morel, ON. Decembre 1989, V 1'
	  write(6,*)'P. Morel, OCA. Octobre 1991, V 2'
	  write(6,*)'P. Morel, OCA. Avril 1993, V 3'
	  write(6,*)'P. Morel, OCA. Janvier 1996, CESAM_3 "pro delphini"'	
	  write(6,*)' '
	  write(6,*)'***********************************'
	  write(6,*)' '
	  stop	!****** F I N  D U  C A L C U L *******
 
	 else			!ecriture
	  call list_3(n,m,p,ro,t,r,l,compg,age,teff,modele,rstar,hp,mstar,
     1	methode,precix,.true.,kap,epsilon,eint,xdot,ydot,rdot,
     2	nh1,nhe1,nhe2,alfa,beta,delta,cp,vaissala,convec,dlpdxt,dlpdxr,
     3	gradconv,gradrad,gradad,dcapdt,dcapdr,depsdr,depsdt,d2p,d2ro,
     4	anupp,anupep,anub8,anube7,anun13,anuo15,lambda,w,z,chaine,
     5	p_atm,t_atm,m_atm,tau,r_atm,ro_atm,k_atm,gradr_atm,grada_atm,
     6	gradc_atm,n_atm,i_pp,i_cno,i_3a,i_gr)
	  agep=age
	 endif
 
c	 ecriture du fichier de reprise car il n'y a pas eu arret
c	 et donc formation d'un *.dat puis ecriture reduite
 
	else
 
c-------------------------------------------------------------------
c
c     Modification pour avoir des fichiers de reprise tous les dt
c
c     Daniel, Novembre 1998
c
c modifrep
 
	 if ( enre_bin .AND. dt_ .ge. dtrep ) then
	  n_fich_rep=n_fich_rep+1

c--------- Cas 0
 
	  if(n_fich_rep .le. 9)then
c	   encode(2,2,l1) n_fich_rep
           write(l2,2) (n_fich_rep)
	   open(unit=26,form='unformatted',status='new',
     1	        file=nom_fich2(:long(nom_fich2))//'_'//
     2          key_word(:long(key_word))//'_'//l2(2:2)//'.repdc')
	   open(unit=66,status='unknown',form='formatted',access='append',
     +          file=key_word(:long(key_word))//'.hr_age')
	   write(66,'(2d22.8,I6,d22.8)') lteff, log10(l(n)), n_fich_rep,
     +           age
	   close(66)
	  elseif(n_fich_rep .le. 99)then
 
c--------- Cas 1

c	   encode(2,2,l2) n_fich_rep
           write(l2,2) (n_fich_rep)
	   open(unit=26,form='unformatted',status='new',
     1	        file=nom_fich2(:long(nom_fich2))//'_'//
     2          key_word(:long(key_word))//'_'//l2//'.repdc')
	   open(unit=66,status='unknown',form='formatted',access='append',
     +          file=key_word(:long(key_word))//'.hr_age')
	   write(66,'(2d22.8,I6,d22.8)') lteff, log10(l(n)), n_fich_rep,
     +           age
	   close(66)
	  elseif(n_fich_rep .le. 999)then
 
c--------- Cas 2
 
c	   encode(3,1,l3)n_fich_rep
           write(l3,1) (n_fich_rep)
	   open(unit=26,form='unformatted',status='new',
     1	        file=nom_fich2(:long(nom_fich2))//'_'//
     2          key_word(:long(key_word))//'_'//l3//'.repdc')
	   open(unit=66,status='unknown',form='formatted',access='append',
     +          file=key_word(:long(key_word))//'.hr_age')
	   write(66,'(2d22.8,I6,d22.8)') lteff, log10(l(n)), n_fich_rep,
     +           age
	   close(66)
	  elseif(n_fich_rep .le. 9999)then
 
c--------- Cas 3
 
c	   encode(4,3,l4)n_fich_rep
           write(l4,3) (n_fich_rep)
	   open(unit=26,form='unformatted',status='new',
     1	        file=nom_fich2(:long(nom_fich2))//'_'//
     2          key_word(:long(key_word))//'_'//l4//'.repdc')
	   open(unit=66,status='unknown',form='formatted',access='append',
     +          file=key_word(:long(key_word))//'.hr_age')
	   write(66,'(2d22.8,I6,d22.8)') lteff, log10(l(n)), n_fich_rep,
     +           age
	   close(66)
	  else
 
           print *, ' '
           print *, '*******************************************'
           print *, ' '
	   write(6,*)' les modeles ne sont plus tous conserves'
           print *, ' '
           print *, '     Nombre de fichiers sup. a 9999 !'
           print *, ' '
           print *, '*******************************************'
           print *, ' '
 
           stop
 
	  endif
	 else
	  open(unit=26,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.rep')
	 endif
	
	 modeleb='m'//cm//'X'//ch//'a'//cr//'_B.rep'
	 write(26)age,nbelem,nchim,mtot,alpha,modeleb,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	 close(unit=26)	!fermeture du fichier du point de reprise
 
c	 ecriture reduite
	
	 call lim_ext(.true.,l(n),r(n),chim1g,pext,text,dpdl,dpdr,dtdl,dtdr,
     1	teff,rstar,mext,dml,dmr,p_atm,t_atm,m_atm,tau,r_atm,mstar,
     2	tdetau,etat,opa)
 
c	 print*,mstar,mtot
	 pertem=mstar/mtot-1.d0
	 if(ihe4 .le. 0)then
	  write(6,205)age,log10(teff),log10(l(n)),log10(r(n)),
     1	logg+log10(m(n)/r(n)**2),
     2	p(1),t(1),ro(1),compg(nbelem*(1-1)+1),
     3	i_pp,i_cno,i_3a,i_gr,pertem,mstar,chaine
	 else
	  write(6,206)age,log10(teff),log10(l(n)),log10(r(n)),
     1	logg+log10(m(n)/r(n)**2),
     2	p(1),t(1),ro(1),compg(nbelem*(1-1)+1),
     3	i_pp,i_cno,i_3a,i_gr,compg(nbelem*(1-1)+ihe4),
     4	pertem,mstar,chaine
	 endif
 
	 if(dessin)call des_cesam(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,teff,
     1	r2,m23,age,lim,lconv,m_zc,mstar)
	
c**********************************************************
clehr
c Determination des zones de combustion nucleaire
c -----------------------------------------------
 
c lpp : variable logique, lpp=.true. si on se trouve dans une region ou
c       \epsilon_{pp} .ge. 10.**3 erg.s-1.g-1
c lcno, l3a, lgrav : idem avec le cycle CNO, la reaction 3 alpha et energie
c gravitationnelle.
 
c n_pp     : nbre de limites non-combustion/combustion pour PP.
c n_cno    : idem pour le cycle CNO.
c n_3a     : idem pour les chaines 3 alpha.
c n_grav_p : idem pour l'energie gravitationnelle positive.
c n_grav_m : idem pour l'energie gravitationnelle negative.
 
c natpp(10) : tableau de variables logiques, exemple : natpp(3)=.true.
c             la troisieme limite non_combu./combu. est le debut d'une
c             zone de combu. quand on va de la surface vers le centre
c             rq.:valeurs max. : 10 !!
 
c pos_lpp(10)   : position en m/M* des limites des zones de combustion (PP)
c pos_cno(10)   : idem pour CNO
c pos_3a(10)    : idem pour 3 alpha
c pos_grav_p(10): idem energie grav. +
c pos_grav_m(10): idem energie grav. -
 
c t_lpp(10)    : temperature des limites des zones PP
c t_lcno(10)   : idem pour CNO
c t_l3a(10)    : idem pour 3 alpha
c t_lgrav_p(10): idem energie grav. +
c t_lgrav_m(10): idem energie grav. -
 
c Initialisation des differentes quantites :
 
      lpp=.false.
      lcno=.false.
      l3a=.false.
      lgrav_p=.false.
      lgrav_m=.false.
 
      n_pp=0
      n_cno=0
      n_3a=0
      n_grav_p=0
      n_grav_m=0
 
c Determination des zones Epsilon_pp > 10.**3
 
      do indi= n, 1, -1
 
         if( epsilon(5*(indi-1)+2) .lt. 1.d+3) then
           if(lpp)then
              lpp=.false.
              n_pp=n_pp+1
              natpp(n_pp)=lpp
              pos_lpp(n_pp)=1.d+0-m(indi)/mstar
              t_lpp(n_pp)=log10(t(indi))
              num_pp(n_pp)=indi
           end if
         else
           if(.not.lpp)then
              lpp=.true.
              n_pp=n_pp+1
              natpp(n_pp)=lpp
              pos_lpp(n_pp)=1.d+0-m(indi)/mstar
              t_lpp(n_pp)=log10(t(indi))
              num_pp(n_pp)=indi
           end if
         end if
 
c Determination des zones Epsilon_cno > 10.**3
 
         if( epsilon(5*(indi-1)+3) .lt. 1.d+3) then
           if(lcno)then
              lcno=.false.
              n_cno=n_cno+1
              natpp(n_cno)=lcno
              pos_lcno(n_cno)=1.d+0-m(indi)/mstar
              t_lcno(n_cno)=log10(t(indi))

              rho_lcno(n_cno)=log10(ro(indi))
              x_lcno(n_cno)=  compchim(nbelem*(indi-1)+1)*nucleo(1) ! H
              y_lcno(n_cno)=  compchim(nbelem*(indi-1)+2)*nucleo(2) ! He3
     +                      + compchim(nbelem*(indi-1)+3)*nucleo(3) ! He4
     +                      +    compchim(nbelem*(indi-1)+5)*nucleo(5) ! C13
     +                      +    compchim(nbelem*(indi-1)+6)*nucleo(6) ! N13
     +                      +    compchim(nbelem*(indi-1)+7)*nucleo(7) ! N14
     +                      +    compchim(nbelem*(indi-1)+8)*nucleo(8) ! N15
     +                      +    compchim(nbelem*(indi-1)+9)*nucleo(9) ! O16
     +                      +    compchim(nbelem*(indi-1)+10)*nucleo(10) ! O17
     +                      +    compchim(nbelem*(indi-1)+11)*nucleo(11) ! O18

              num_cno(n_cno)=indi
           end if
         else
           if(.not.lcno)then
              lcno=.true.
              n_cno=n_cno+1
              natcno(n_cno)=lcno
              pos_lcno(n_cno)=1.d+0-m(indi)/mstar
              t_lcno(n_cno)=log10(t(indi))

              rho_lcno(n_cno)=log10(ro(indi))
              x_lcno(n_cno)=  compchim(nbelem*(indi-1)+1)*nucleo(1) ! H
              y_lcno(n_cno)=  compchim(nbelem*(indi-1)+2)*nucleo(2) ! He3
     +                      + compchim(nbelem*(indi-1)+3)*nucleo(3) ! He4
     +                      +    compchim(nbelem*(indi-1)+5)*nucleo(5) ! C13
     +                      +    compchim(nbelem*(indi-1)+6)*nucleo(6) ! N13
     +                      +    compchim(nbelem*(indi-1)+7)*nucleo(7) ! N14
     +                      +    compchim(nbelem*(indi-1)+8)*nucleo(8) ! N15
     +                      +    compchim(nbelem*(indi-1)+9)*nucleo(9) ! O16
     +                      +    compchim(nbelem*(indi-1)+10)*nucleo(10) ! O17
     +                      +    compchim(nbelem*(indi-1)+11)*nucleo(11) ! O18

              num_cno(n_cno)=indi
           end if
         end if
 
c Determination des zones Epsilon_3a > 10.**3
 
         if( epsilon(5*(indi-1)+4) .lt. 1.d+3) then
           if(l3a)then
              l3a=.false.
              n_3a=n_3a+1
              nat3a(n_3a)=l3a
              pos_l3a(n_3a)=1.d+0-m(indi)/mstar
              t_l3a(n_3a)=log10(t(indi))
              num_3a(n_3a)=indi
           end if
         else
           if(.not.l3a)then
              l3a=.true.
              n_3a=n_3a+1
              nat3a(n_3a)=l3a
              pos_l3a(n_3a)=1.d+0-m(indi)/mstar
              t_l3a(n_3a)=log10(t(indi))
              num_3a(n_3a)=indi
           end if
         end if
 
c Determination des zones Epsilon_grav > + 10.**2
 
         if( epsilon(5*(indi-1)+5) .lt. 1.d+2) then
           if(lgrav_p)then
              lgrav_p=.false.
              n_grav_p=n_grav_p+1
              natgrav_p(n_grav_p)=lgrav_p
              pos_lgrav_p(n_grav_p)=1.d+0-m(indi)/mstar
              t_lgrav_p(n_grav_p)=log10(t(indi))
              num_grav_p(n_grav_p)=indi
           end if
         else
           if(.not.lgrav_p)then
              lgrav_p=.true.
              n_grav_p=n_grav_p+1
              natgrav_p(n_grav_p)=lgrav_p
              pos_lgrav_p(n_grav_p)=1.d+0-m(indi)/mstar
              t_lgrav_p(n_grav_p)=log10(t(indi))
              num_grav_p(n_grav_p)=indi
           end if
         end if
 
c Determination des zones Epsilon_grav < - 10.**2
 
         if( epsilon(5*(indi-1)+5) .gt. 1.d+2) then
           if(lgrav_m)then
              lgrav_m=.false.
              n_grav_m=n_grav_m+1
              natgrav_m(n_grav_m)=lgrav_m
              pos_lgrav_m(n_grav_m)=1.d+0-m(indi)/mstar
              t_lgrav_m(n_grav_m)=log10(t(indi))
              num_grav_m(n_grav_m)=indi
           end if
         else
           if(.not.lgrav_m)then
              lgrav_m=.true.
              n_grav_m=n_grav_m+1
              natgrav_m(n_grav_m)=lgrav_m
              pos_lgrav_m(n_grav_m)=1.d+0-m(indi)/mstar
              t_lgrav_m(n_grav_m)=log10(t(indi))
              num_grav_m(n_grav_m)=indi
           end if
         end if
 
      end do
 
c Fin de la determination des zones "energetiques"
c------------------------------------------------------------
c Caracterisation du coeur d'helium, determination de Rc, Mc,
c Lc et Phic

c Rq.: 'Phi_c' n'est pas exactement le potentiel du coeur d'helium
c      mais une quantite qui s'en approche.
      log_lum_c=0.d0
      masse_c=0.d0
      rayon_c=0.d0
      log_phi_c=0.d0	
      l_c=.false.

      xsurf=compchim(nbelem*(n-1)+1)*nucleo(1) ! X a la surface.

      do indi= n, 1, -1
 
      if( compchim(nbelem*(indi-1)+1)*nucleo(1) .le. 0.5d0*xsurf 
     +    .AND. (.NOT. l_c) )then 
        l_c=.true.
        rayon_c=r(indi) ! r est en rsol
        log_lum_c=log10(l(indi))
        masse_c=m(indi) ! m est en Msol
        log_phi_c=log10(masse_c/rayon_c) 
      end if
 
      end do

c------------------------------------------------------------
c Determination des pourcentages de production d'energie
c
c      ener_pp=0.d0
c      ener_cno=0.d0
c      ener_3a=0.d0
c      ener_grav=0.d0
c
c      do indi= 1, n-1
c
c         dmasse=m(indi+1)-m(indi)
c
c         ener_pp=ener_pp+epsilon(5*(indi-1)+2)*dmasse
c         ener_cno=ener_cno+epsilon(5*(indi-1)+3)*dmasse
c         ener_3a=ener_3a+epsilon(5*(indi-1)+4)*dmasse
c         ener_grav=ener_grav+epsilon(5*(indi-1)+5)*dmasse
c
c      end do
c
c      ener_tot=ener_pp+ener_cno+ener_3a+ener_grav
c
c      i_pp=nint(100.*ener_pp/ener_tot)
c      i_cno=nint(100.*ener_cno/ener_tot)
c      i_3a=nint(100.*ener_3a/ener_tot)
c      i_gr=nint(100.*ener_grav/ener_tot)
c
c------------------------------------------------------------------------
c Determination du facteur de degenerescence Lambda au centre
 
      do indi= 1, nbelem  ! Composition chimique au centre
         l_xchi(indi)=compchim(nbelem*(1-1)+indi)*nucleo(indi)
      end do
 
      call etat(p(1),t(1),l_xchi,.false.,
     +          bid1,bid2,bid3,bid4,bid5,bid6,bid7,bid8,bid9,
     +          bid10,bid11,bid12,bid13,bid14,bid15,bid16,bid17,
     +          lambda_c)
 
c-----------------------------------------------------------------------
c Determination de la zone Z > 0.1 :
 
      z_lim=0.1d0
 
      do indi= n, 1, -1
         z_local=1.d0-compchim(nbelem*(indi-1)+1)*nucleo(1)
     +           -compchim(nbelem*(indi-1)+2)*nucleo(2)
     +           -compchim(nbelem*(indi-1)+3)*nucleo(3)
         if( z_local .ge. z_lim ) then
           change_compo=n
           goto 1211
         end if
      end do
1211   m_change_compo=1.d0-m(change_compo)/mstar
 
c-----------------------------------------------------------------------
c      Les Ecritures
 
c        10        20        30        40        50        60        70
      write(53,2201)'-------------------------------------------------'
2201  format(a)
 
      write(53,2101) age, n ! age et nombre de couches
2101  format(1pd22.15,i6)
 
      write(53,2102) mstar, log10(r(n)), log10(l(n)), log10(teff)
2102  format(1p4d13.6)
 
      write(53,2112) i_pp, i_cno, i_3a, i_gr ! les % des diff. sources d'energie
2112  format(1p,4I3)
 
      write(53,2103) log10(t(1)), log10(p(1)), log10(ro(1)),
     +          lambda_c
 
      write(53,2133) (compchim(nbelem*(1-1)+i)*nucleo(i),i=1,18)
2103  format(1p4d13.6)
2133  format(1p18d13.6)
 
c ---- Composition en surface ---------------------
 
      write(53,2113) (compchim(nbelem*(n-1)+i)*nucleo(i),i=1,18)
2113  format(1p18d13.6)
 
c ---- Caracteristiques du coeur d'helium ------------------------------
 
      write(53,2124) rayon_c, masse_c, log_lum_c, log_phi_c, 'Coeur He'
2124  format(1p4d13.6,A9)
 
c ---- Profondeur a laquelle il y changement de composition chimique ---
c (par rapport a la compo initiale du modele homogene!!)
 
      write(53,2123) m_change_compo
2123  format(1pd13.6)
 
c ---- Zones de combustion PP ---------------------
 
      write(53,2104) n_pp, (natpp(i),i=1,n_pp)
2104  format(1p,i3,40l2)
      do i= 1, n_pp
         write(53,2115) pos_lpp(i), t_lpp(i), num_pp(i)
      end do
2105  format(1p2D13.6)
 
c ---- Zones de combustion CNO -------------------
 
      write(53,2104) n_cno, (natcno(i),i=1,n_cno)
      do i= 1, n_cno
         write(53,2115) pos_lcno(i), t_lcno(i), num_cno(i), rho_lcno(i), 
     +                  x_lcno(i), y_lcno(i), Xcno_lcno(i) 
2115  format(1p2d13.6,i5,4d16.6)
      end do
 
c ---- Zones de combustion 3 alpha -------------------
 
      write(53,2104) n_3a, (nat3a(i),i=1,n_3a)
      do i= 1, n_3a
         write(53,2115) pos_l3a(i), t_l3a(i), num_3a(i)
      end do
 
c ---- Zones d'energie gravita. + -------------------
 
      write(53,2104) n_grav_p, (natgrav_p(i),i=1,n_grav_p)
      do i= 1, n_grav_p
         write(53,2115) pos_lgrav_p(i), t_lgrav_p(i), num_grav_p(i)
      end do
 
c ---- Zones d'energie gravita. - -------------------
 
      write(53,2104) n_grav_m, (natgrav_m(i),i=1,n_grav_m)
      do i= 1, n_grav_m
         write(53,2115) pos_lgrav_m(i), t_lgrav_m(i), num_grav_m(i)
      end do
 
c ---- Zones Convectives --------------------------
 
      write(53,2106) lim, (lconv(i),i=1,lim)
2106  format(i3,20l2)
      do i= 1, lim
         write(53,2107) 1.d0-m_zc(i), r_zc(i), r_ov(i), logtzc(i), 
     +                  logtov(i), logrozc(i)
2107  format(1p5d13.6,d16.6)
      end do

c*********************************************************
c Ecriture dans le fichier pour Roxburgh
c Roxy
        if ( ROX ) then
	write(76,1076) log10(teff), diff, mroxSmsch, 
     +      compchim(nbelem*(1-1)+1)*nucleo(1), age
        end if
c*********************************************************
	 call list_3(n,m,p,ro,t,r,l,compg,age,teff,modele,rstar,hp,mstar,
     1	methode,precix,.false.,kap,epsilon,eint,xdot,ydot,rdot,
     2	nh1,nhe1,nhe2,alfa,beta,delta,cp,vaissala,convec,dlpdxt,dlpdxr,
     3	gradconv,gradrad,gradad,dcapdt,dcapdr,depsdr,depsdt,d2p,d2ro,
     4	anupp,anupep,anub8,anube7,anun13,anuo15,lambda,w,z,chaine,
     5	p_atm,t_atm,m_atm,tau,r_atm,ro_atm,k_atm,gradr_atm,grada_atm,
     6	gradc_atm,n_atm,i_pp,i_cno,i_3a,i_gr)
	endif				!ecritures
 
c---------------------------------------------------------------
 
c		E V O L U T I O N	T E M P O R E L L E
 
c---------------------------------------------------------------
 
c	pause'avant evolution'
	if(ipms .eq. 1)then
	 rp1=r(n)
	 lp1=l(n)
	 write(6,55)c_iben,rp1,lp1,t(1)
55	 format(1x,/,1x,'Pre Sequence Principale premier modele: C=',1pd10.3,/,
     1	1x,'Rext=',1pd10.3,', Lext=',1pd10.3,', Tc=',1pd10.3,/,
     2	1x,' ok? o/n')
	 read(5,'(a)')oui
c	 write(6,*)'oui ',oui
	 if(oui .eq. 'o')then
	  c_iben=c_iben*1.1
	  ipms=2
	  write(6,*)'calcul d''un nouveau modele avec 1.1 C'
	 else
	  write(6,*)' '
	  write(6,*)'entrer une nouvelle valeur pour C'
	  read(5,*)c_iben
	  ipms=1
	 endif
	 dt=0.
	 call resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,iben_3,lim_ext,tdetau,coeff_diff,perte_m)
	
	elseif(ipms .eq. 2)then
	 rp2=r(n)
	 lp2=l(n)
	 dt=mtot**2*g*msol/lsol/rsol*msol*2.*abs(rp1-rp2)/
     1	rp1/rp2/(lp1+lp2)/secon6
	 write(6,57)c_iben,rp2,lp2,t(1),dt
57	 format(1x,/,1x,'Pre Sequence Principale second modele: C=',1pd10.3,/,
     1	1x,'Rext=',1pd10.3,', Lext=',1pd10.3,', Tc=',1pd10.3,
     2	', dt=',1pd10.3,/,1x,' ok? o/n')
	 read(5,'(a)')oui
	 if(oui .eq. 'n')then		!nouvelle cte de contraction
	  write(6,*)' '
	  write(6,*)'entrer une nouvelle valeur pour C'
	  read(5,*)c_iben
	  ipms=1
	  age=0
	  dt=0.
	  call resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,iben_3,lim_ext,tdetau,coeff_diff,perte_m)
 
	 else				!on est content
	  ipms=3
	  write(6,*)' '
	  write(6,*)'debut d''evolution PMS'
	  write(6,*)' '
	  t_inf=t_inf0
	  age=0.
	  dt_=dt	!permet d'utiliser le modele ***_t pour calcul du TdS
	  open(unit=24,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.pms')
	  modeleb='m'//cm//'X'//ch//'a'//cr//'_B.pms'
	  write(24)age,nbelem,nchim,mtot,alpha,modeleb,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	  close(unit=24)
	  call resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,nuc,lim_ext,tdetau,coeff_diff,perte_m)
	 endif
	else
c	 pause'appel a resout_3'
	 call resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,nuc,lim_ext,tdetau,coeff_diff,perte_m)
	endif
 
	age=age+dt
 
	goto 100
 
	end
