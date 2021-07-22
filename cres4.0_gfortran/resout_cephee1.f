c********************************************************************

c Modif. septembre 2003, Varsovie : prise en compte des variations
c de r et L dans la fonction de répartition, introduction de cter
c et ctel.

c Version de 'resout_3' avec des reprises automatiques
c Daniel, Septembre 1997
 
c Version avec possibilite d'enregistrer des modeles tous les dt
c nouveaux parametres d'arret : Y_STOP et LOG_LSLSOL
c Daniel, Novembre 1998

c Modif. de la gestion du pas temporel : on augmente dt de 10%
c si Yc augmente sur une BL, mot cle : dty
c********************************************************************
 
	subroutine resout_cephee1(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,un23,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,m23t_t,knot23_t,tdst,
     6	ctes,etat,opa,conv,nuc,lim_ext,tdetau,coeff_diff,perte_m)
	
c	gestion de la resolution des equations avec fonction de repartition
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	version: 3
c	05 09 96, suppression de evite_pms
c	17 09 96 introduction de la "sa" solar accuracy
c	07 03 97 : suppression de l'ecriture du format 111
 
c Pour le probleme aux limites:
c	resolution du systeme d'equations differentielles non lineaires
c	de la structure interne par developpement sur B-splines
c	par iteration Newton avec derivees analytiques
 
c	les equations de la SI sont ecrites pour tout point dans static_m_3
c	le jacobien de la methode NR est forme et resolu dans resout3 dans
c	l'espace des splines
 
c	ordre des splines pour l'integration par collocation m_qs+1,
 
c	pour la composition chimique m_ch
 
c	Pour le probleme de valeurs initiales appel a evol_3
 
c entree
c	un23 :  1 poursuite d'une evolution, modele repris en binaire
c		2 modele initial en binaire, ZAMS
c		3 modele initial en ASCII, ZAMS
c	       -2 modele initial en binaire, PMS
c	       -3 modele initial en ASCII, PMS
 
c	age : age du modele
 
c entree/sortie
c	n : nombre de points de raccord
c	q,qt,knot : points de raccord et de table et nombre de points
c	bp : solution spline
c	dt : pas temporel
c	m_t,mct_t,nc_t,knotc_t,chim_t:   comp. chim. au temps t
c	mc,mct,nc,knotc,chim: comp. chim. a t+dt
c	n,bp,q,qt,knot,p,t,r,l,m,grad_mj:  var. pples. a t+dt
c	lim : nombre de limites ZR/ZC
c	lim_t : nombre de limites ZR/ZC au temps t
c	jlim : abscisses des limites ZR/ZC
c	mm,mt,mr,nr,knot,tdst: elements pour interpoler le TdS au temps t-dt
c	m_zct : masses des limites des ZR/ZC au temps t
c	lconv_t =.true. : debut de ZC au temps t
c	r_zc,r_ov : rayons de limites ZR/ZC et d'overshoot
	
c sortie
c	m_zc : masses des limites des ZR/ZC au temps t+dt
c	lconv =.true. : debut de ZC au temps t+dt
c	mstar: masse avec perte de masse
 
 
c	dts : estimation du pas temporel suivant
c	reprise=.true. : on poursuit une evolution
 
c routines externes
c	ctes, etat, opa, conv, nuc, lim_ext, tdetau, coeff_diff
 
c	le reseau de points de raccord q, possede n points
 
c	la solution, les coefficients des splines sur le
c	reseau de points de table qt, sont dans bp
 
c	5 09 96, suppression de evite_pms P. Morel
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'ctephy.common'
	include 'atmosphere_3.common'
	
	integer n,ni(pn),spi,ip,i,li,der,spl,indpc(pbp),indcol,mpr1,n_min,ini0,
     1	eq,nl,var,ligne,indice,ipe,vare,compt,cx,rp1,j,k,dim,new_n,
     2	nc_t,knotc_t,nc,knotc,kmax,knot,lim,lim_t,ncouche,iter_max,
     3	knot_t,n_t,m_qs1,m_ch1,nbelem1,un23,lx,bloc,knot23_t,iter_max0,
     4	jlim(1),jlim_t(1),n_,nc_,knot_,knotc_,n__,nc__,knot__,knotc__,
     5	nm,nm_,nm__,knotm,knotm_,knotm__
c	data rp1/2/
 
	integer*4 long

	external long

	real*8	a((pm_qs+1)*pne*pbp),b(pbp),dxchim(pnelem),esp(pn),
     1	xl(pne),xx(pqt),derxx(2*(pm_qs+1)*pqt),alpha1,m(pn),t(pn),
     2	lderxx(2*(pm_qs+1)*pne),estim(pnelem),p(pn),r(pn),l(pn),
     3	ae(2*pne**2),be(pne),y(pne*2),corr,err,er,nrm,err1,
     4	fac(pn),mstar,tdst(1),mstar_t,bid,err4,xchim0(pnelem),
     5	old_m23(pn),new_m23(pn),new_m23t(pqt),nuc_m(pnelem),
     6	old_m23_(pn),new_m23_(pn),new_m23t_(pqt),errp,	
     7	old_m23__(pn),new_m23__(pn),new_m23t__(pqt)	
c	data dxchim,err/pnelem*0.d0,1.d3/
	
	real*8	bp(1),q(1),qt(1),chim(1),mc(1),mct(1),m23t_t(1),
     1	r2(1),m23(1),dts,m_zc(1),r_zc(1),r_ov(1),age,psi0,dpsi,
     2	bp_t(1),q_t(1),qt_t(1),chim_t(1),mc_t(1),mct_t(1),
     3	r2_t(1),m23_t(1),dt,m_zc_t(1),r_zc_t(1),r_ov_t(1),
     4	bp_(1),q_(1),qt_(1),chim_(1),mc_(1),mct_(1),r2_(1),m23_(1),dt_,
     5	bp__(pbp),q__(pn),qt__(pqt),chim__(pchim),mc__(pnch),
     6	mct__(pchimt),r2__(pn),m23__(pn),dt__,mtot1,y01,x01,
     7	xchim(pnelem),f(pne),dfdr(pne),w_rot_p,epsilon(5),ro,dcomp,jac,
     8	depsro,hh,be7,b8,n13,o15,f17,depst,depsx(pnelem)
 
	logical init,g_max,lconv(1),lconv_t(1),logic,
     1	ini,reprise,new,inversible,
     2	diffusion1,ovsht,rot_solid_p,z_cte_p
c	data init/.true./

c dty

	logical aug_dt_y

 
	character*70 sub_phys(2)
	common/subphys/sub_phys
	
	save
 
	character*1 oui
	character*2 ci,cj
	character*3 ch,cr
	character*4 cm,ck
	character*9 today
	character*50 modelep
	character*80 titre
 
	external ctes,etat,opa,conv,nuc,lim_ext,tdetau,coeff_diff,perte_m
 
	data rp1/2/
	data dxchim,err/pnelem*0.d0,1.d3/
	data init/.true./
 
2000	format((1x,1p8d10.3))
 
	if(init)then
	 init=.false.
c	 write(2,141)sub_phys(1)(:long(sub_phys(1))),
c	1	sub_phys(2)(:long(sub_phys(2)))
c	 write(6,141)sub_phys(1)(:long(sub_phys(1))),
c	1	sub_phys(2)(:long(sub_phys(2)))
	 write(2,141)sub_phys(1),sub_phys(2)
	 write(6,141)sub_phys(1),sub_phys(2)
141	 format(t14,'*************************************************',//,
     1	t15,'MODELE DE STRUCTURE INTERNE calcule par cesam_3',//,
     2	t14,'*************************************************',///,
     3	t10,'Arguments de CESAM3:',//,t5,a70,/,t5,a70,//,	
     4	t10,'DONNEES des NAMELISTS',/)
	
c dty

	 aug_dt_y=.false.

	 if(abs(un23) .le. 2)then	!pour identification des modeles repris
	  mtot1=mtot
	  mstar_t=mstar
	  y01=y0
	  x01=x0
	  alpha1=alpha
	  m_qs1=m_qs
	  mpr1=m_qs1+1
	  m_ch1=m_ch
	  nbelem1=nbelem
	  diffusion1=diffusion
	  rot_solid_p=rot_solid
	  w_rot_p=w_rot
	  z_cte_p=z_cte
	 endif
 
c	 lecture des namelists
 
c	 call lit_nl_3	!lecture des NAMELISTs
 
         call lit_nl_cephee1
	
         psi0=psi0_c
 
c	 choix des pas temporels initiaux pour evolution a partir de la ZAMS
 
	 if(mtot .le. 1.1)then
	  dt0=10.
	 else
	  dt0=1.
	 endif
	
c	 on fixe les parametres numeriques	
	 	
	 if(precision .eq. 'sp')then
	  m_qs=2
	  m_ch=3
	  ordre=2
c	  precix=1.d-4
c	  precit=.2
	  ro_test=.1
c	  psi0=.08
c	  d_grav=.5
c	  dtmax=100.
	  ini0=6
	  n_atm=30
	  n23=20
	  kipp=.false.		!TdS=dE+PdV
	  write(6,215)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav
	  write(2,215)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav	
215	  format(/,1x,'----------parametres numeriques utilises-----------',//,
     1	1x,'modele calcule avec super precision',//,
     2	1x,'m_qs=',i2,', m_ch=',i2,', ordre=',i2,', precix=',1pd8.1,/,
     3	1x,'precit=',1pd8.1,', ro_test=',1pd8.1,', psi0=',1pd8.1,
     4	', d_grav=',1pd8.1,/)
 
	 elseif(precision .eq. 'ds')then	!delta scuti
	  m_qs=2
	  m_ch=2
	  ordre=2
c	  precix=1.d-4
c	  precit=.2
	  ro_test=.1
c	  psi0=.09
c	  d_grav=.5
c	  dtmax=100.
	  ini0=5
	  n_atm=30
	  n23=20
	  kipp=.true.	!TdS simplifie de Kippenhahn
	  write(6,301)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav
	  write(2,301)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav	
301	  format(/,1x,'----------parametres numeriques utilises-----------',//,
     1	1x,'modele calcule avec precision delta-scuti',//,
     2	1x,'m_qs=',i2,', m_ch=',i2,', ordre=',i2,', precix=',1pd8.1,/,
     3	1x,'precit=',1pd8.1,', ro_test=',1pd8.1,', psi0=',1pd8.1,
     4	', d_grav=',1pd8.1,/)
	
	 elseif(precision .eq. 'sa')then	!Solar accuracy
	  m_qs=2
	  if(diffusion)then
	   m_ch=3
	  else
	   m_ch=4
	  endif
	  ordre=2
c	  precix=5.d-6
c	  precit=.1
	  ro_test=.01
c	  psi0=.06
c	  d_grav=.1
c	  dtmax=50.
	  ini0=6
	  n_atm=50
	  n23=30
	  dt0=1.
	  kipp=.false.		!TdS=dE+PdV
	  write(6,300)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav
	  write(2,300)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav	
300	  format(/,1x,'----------parametres numeriques utilises-----------',//,
     1	1x,'modele calcule avec precision solaire',//,
     2	1x,'m_qs=',i2,', m_ch=',i2,', ordre=',i2,', precix=',1pd8.1,/,
     3	1x,'precit=',1pd8.1,', ro_test=',1pd8.1,', psi0=',1pd8.1,
     4	', d_grav=',1pd8.1,/)
		
	 elseif(precision .eq. 'hp')then
	  m_qs=2
	  m_ch=4
c	  if(mdot .eq. 0.d0)then	!avec perte de masse ordre .le. 2
c	   ordre=4			!suppression
c	  else
c	   ordre=2
c	  endif
	  if(iw .gt. 1)ordre=2
	  ordre=2
c	  precix=1.d-6
c	  precit=.1
	  ro_test=.01
c	  psi0=.09
c	  d_grav=.1
c	  dtmax=50
	  ini0=6
	  n_atm=50
	  n23=30	
	  kipp=.false.		!TdS=dE+PdV
	  write(6,217)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav
	  write(2,217)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav	
217	  format(/,1x,'----------parametres numeriques utilises-----------',//,
     1	1x,'modele calcule avec hyper precision',//,
     2	1x,'ATTENTION : risques d''instabilites, cout eleve',
     2	1x,'m_qs=',i2,', m_ch=',i2,', ordre=',i2,', precix=',1pd8.1,/,
     3	1x,'precit=',1pd8.1,', ro_test=',1pd8.1,', psi0=',1pd8.1,
     4	', d_grav=',1pd8.1,/)
	
	 elseif(precision .eq. 'cp')then	!precision de cepheide
	  m_qs=1
	  m_ch=2
	  ordre=2
c	  precix=5.d-4
c	  precit=.12
	  ro_test=.001
c	  psi0=.07
c	  d_grav=.075
c	  dtmax=10.
	  dt0=.05
	  ini0=5
	  n_atm=30
	  n23=20
	  kipp=.true.		!TdS simplifie de Kippenhahn
	  write(6,302)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav
	  write(2,302)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav	
302	  format(/,1x,'----------parametres numeriques utilises-----------',//,
     1	1x,'modele avec precision cepheide',/,
     2	1x,'m_qs=',i2,', m_ch=',i2,', ordre=',i2,', precix=',1pd8.1,/,
     3	1x,'precit=',1pd8.1,', ro_test=',1pd8.1,', psi0=',1pd8.1,
     4	', d_grav=',1pd8.1,/)
 
	 elseif(precision .eq. 'c1')then	!precision de cepheide basse
	  m_qs=1

          ! Anciennes valeurs utilisées jusqu'à septembre 2003 :
	  m_ch=2
	  ordre=1
          ! Nouvelles valeurs suite aux mésaventures polonaises avec de
          ! très petits pas de temps (1.d-3 Myrs)
        !m_ch  = 4
        !ordre = 2

c	  precix=5.d-3
c	  precit=.15
	  ro_test=.01
c	  psi0=.1
c	  d_grav=.5
c	  dtmax=10.
	  dt0=.1
	  ini0=5
	  n_atm=30
	  n23=20
	  kipp=.true.		!TdS simplifie de Kippenhahn	
	  write(6,303)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav
	  write(2,303)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav	
303   format(/,1x,'----------parametres numeriques utilises-----------',
     &//,1x,'>>> modele avec precision cepheide basse',/,
     &1x,'m_qs=',i2,', >>m_ch=',i2,', >>ordre=',i2,', precix=',1pd8.1,/,
     &1x,'precit=',1pd8.1,', ro_test=',1pd8.1,', psi0=',1pd8.1,
     &', d_grav=',1pd8.1,/)
	
	 else
	  m_qs=1
	  m_ch=2
	  ordre=1
c	  precix=5.d-4
c	  precit=.3
	  ro_test=1.
c	  psi0=.11
c	  d_grav=1.
c	  dtmax=200.
	  ini0=4
	  n_atm=25
	  n23=15
	  kipp=.true.		!TdS simplifie de Kippenhahn	
	  write(6,216)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav
	  write(2,216)m_qs,m_ch,ordre,precix,precit,ro_test,psi0,d_grav	
216	  format(/,1x,'----------parametres numeriques utilises-----------',//,
     1	1x,'modele calcule avec precision normale',/,
     2	1x,'m_qs=',i2,', m_ch=',i2,', ordre=',i2,', precix=',1pd8.1,/,
     3	1x,'precit=',1pd8.1,', ro_test=',1pd8.1,', psi0=',1pd8.1,
     4	', d_grav=',1pd8.1,/)
	
	 endif
	
 
c         write(6,*) 'repri= ', repri
c         write(6,*) 'd_grav= ', d_grav
c         stop
 
	 call ctes		!initialisation des constantes physiques
	 	
c	 identification du modele : masse 1.0 msol.,X hydrogene 0.70 (z=0.02)
c	 alpha l.mel. 1.5 : m100X700a150
 
	 modelep=modele				!modele initial
c	 encode(4,6,cm)int(mtot*100.)
         write(cm,6) (int(mtot*100.))
c	 encode(3,1,ch)int(x0*1000.)
         write(ch,1) (int(x0*1000.))
c	 encode(3,1,cr)int(alpha*100.)
         write(cr,1) (int(alpha*100.))
1	 format(i3)
6	 format(i4)
 
	 modele='m'//cm//'X'//ch//'a'//cr//'.col'
	
c------------- date sur SUN, ULTRIX -----------------------
 
c	 call idate(i,j,k)
c	 write(6,*)'i,j,k',i,j,k
c	 encode(2,2,ci)i		!fictifs
         write(ci,2) (i)
c	 encode(2,2,cj)j
         write(cj,2) (j)
c	 encode(4,6,ck)k
         write(ck,6) (k)
2	 format(i2)
c	 titre='Modele de SI: '//modele(:long(modele))//
c	1	' du '//ci//'/'//cj//'/'//ck//' - CESAM3 (P. Morel, OCA)'
	
c---------------------------- sur HP ---------------------------
 
c	 call date(today)
c	 write(6,*)today
c	 pause'date'
	
	 titre='Modele de SI: '//modele(:long(modele))//
     1	' du '//today//' - CESAM (P. Morel, OCA)'
	
c-----------------------------------------------------------------
	
	
	 write(6,*)'--------------identification du modele-----------------'
	 write(6,*)
	 write(6,*)titre(:long(titre))
	 write(6,*)
	
	 write(2,*)'--------------identification du modele-----------------'
	 write(2,*)' '	 	
	 write(2,*)titre(:long(titre))
	 write(2,*)' '
 
	 if(m_qs .gt. pm_qs)then
	  write(6,*)'ordre B-splines pour eq. quasi-stat. ramene a ',pm_qs
	  write(6,*)'au besoin changer pm_qs dans le fichier cesam_3.parametres'
	  m_qs=pm_qs
	 endif
	 if(m_ch .gt. pm_ch)then
	  write(6,*)'ordre B-splines pour interp. comp. chim. ramene a ',pm_ch
	  write(6,*)'au besoin changer pm_ch dans le fichier cesam_3.parametres'
	  m_ch=pm_ch
	 endif
	
c	 composition chimique, rotation, Z
c	 initialisation du vecteur de comp. chim. generalise
c	 ses elements purement chimiques sont initialises et identifies
c	 (nom_elem, ihe4) dans la routine de reactions nucleaires (NUC)
c	 au retour dans resout_3, ce vecteur est complete.
 
	 write(2,*)'---------- REACTIONS NUCLEAIRES -------------'
	 write(6,*)'---------- REACTIONS NUCLEAIRES -------------'	
	 call nuc(t(1),ro,xchim0,dcomp,jac,.false.,1,
     1	epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	 nbelem=nchim		!longeur du vecteur de comp. chim. generalisee
 
c	 on complete le vecteur de composition chimique generalisee
c	 pour la rotation
 
	 write(6,*)
	 write(6,*)
	 write(6,*)'-----------ROTATION et ELEMENTS LOURDS (Z) ------------'
	 write(2,*)'-----------ROTATION et ELEMENTS LOURDS (Z) ------------'
 	
c	 rotation initiale w_rot .ge. 0.
 
c	 rot_solid = .true. : vitesse angulaire w_rot=constante, iw<0
c	 en tout point on prend w_rot dans le common /modele_3/
 
c	 rot_solid = .false. : conservation, ou diffusion du moment angulaire
c	 dans ZR, melange dans ZC, iw>1 Mt angulaire dans xchim(iw)
 
c	 si w_rot lue est > 1 on suppose que c'est en jours
 
	 if(w_rot .gt. 1.d0)then
	  w_rot=2.*pi/24./3600./w_rot		!rad/sec si lue en jours
	 elseif(w_rot .lt. 0.)then
	  w_rot=0.
	 endif
	 if(rot_solid)then
	  write(6,22)w_rot
	  write(2,22)w_rot
22	  format(//,t2,'Rotation solide, vitesse angulaire =',1pd10.3,
     1	' rd/sec')
	  iw=-100
	 else
	  write(6,23)w_rot
	  write(2,23)w_rot
23	  format(//,t2,'Melange dans ZC du MW, vitesse angulaire initiale =',
     1	1pd10.3,' rd/sec')
	  if(w_rot .eq. 0.)then
	   write(6,223)
	   write(2,223)
223	   format(1x,'on effectuera le melange du MW bien que w_rot=0',//)
	  endif	
	  nbelem=nbelem+1
	  iw=nbelem
	  ab_ini(iw)=w_rot
	  nucleo(iw)=1.
	  nom_elem(iw)=' MW'
	  ab_min(iw)=1.d-5		!fictif
	 endif
	
c	 on complete le vecteur de composition chimique generalise
c	 pour les elements lourds Z
 
c	 x0, y0 teneur initiale en H et He4
c	 z0: abondance initiale en elements lourds z0=1-x0-y0
c	 avec diffusion de Z on tient alors plus correctement compte de
c	 l'abondance des metaux dans le calcul de l'opacite
c	 si z_cte:
c	   Z=Z0 est constant c'est la valeur prise pour l'opacite
c	 sinon :
c	   avec diffusion Z est diffuse
c	   sans diffusion Z = 1 - X - Y
 
c	 Pour l'equation d'etat on prend Z = 1 - X - Y
 
c	 Pour le mu_mol et le nombre d'electrons libres on prend Y=1-X-Z0
 
	 z0=1.d0-x0-y0	!initialisation de z0
	 if(z0 .lt. 0.)then
	  print*,'ERREUR : Z0 = 1-X0-Y0 < 0, arret'
	  stop
	 endif
	 write(6,21)z0,z0/x0
	 write(2,21)z0,z0/x0
21	 format(//,t3,'Z initial',1pd10.3,', Z / X initial:',1pd10.3)
 
	 if(z_cte)then
	  write(6,*)' pour l''opacite Z = Z0 est garde constant'
	  write(2,*)' pour l''opacite Z = Z0 est garde constant'
	  iz=-100	
	 elseif(diffusion)then
	  write(6,*)' Z est diffuse'
	  write(2,*)' Z est diffuse'
	  if(nom_elem(nchim) .eq. ' Ex')then
	   write(6,*)'comme un ensemble d''elements trace'
	   write(2,*)'comme un ensemble d''elements trace'
	   iz=-100
	  else 	
	   write(6,*)'comme un element trace unique'
	   write(2,*)'comme un element trace unique'
	   nbelem=nbelem+1
	   iz=nbelem
	   nucleo(iz)=1.
	   nom_elem(iz)=' Z '
	   ab_ini(iz)=z0
	   ab_min(iz)=1.d-5		!fictif
	  endif
	 else
	  write(6,*)' pour l''opacite Z = 1 - X - Y'
	  write(2,*)' pour l''opacite Z = 1 - X - Y'
	  iz=-100
	 endif	
 
	 err1=1.d-2		!valeur minimale pour correction NR modifiee
	 err4=1.		!correction max apres la 4ieme iteration
c	 dtmin=1.d-8		!pas temporel minimum (dans *.don3 maintenant (16/02/99)
 
	 write(6,*)' '
	 write(6,*)'------ Evolution de la composition chimique  --------- '
	 write(2,*)' '
	 write(2,*)'------ Evolution de la composition chimique  --------- '
	 if(diffusion)then
	  write(6,*)' '
	  write(6,*)'Diffusion de la comp. chim. dans ZR et ZC'
	  write(2,*)' '
	  write(2,*)'Diffusion de la comp. chim. dans ZR et ZC'
	 else
	  write(6,*)' '
	  write(6,*)'Modele sans diffusion'
	  write(2,*)' '
	  write(2,*)'Modele sans diffusion'	
	 endif
	
	 write(6,*)' '
	 write(6,*)'-----------perte de masse----------------'
	 write(6,*)' '
	 write(2,*)' '
	 write(2,*)'-----------perte de masse----------------'
	 write(2,*)' '
	
	 if(mdot .eq. 0.d0)then	
	  write(6,*)'modele sans perte de masse'
	  write(2,*)'modele sans perte de masse'	
	 elseif(mdot .lt. 0.d0)then
	  write(6,212)mdot
	  write(2,212)mdot	
212	  format(1x,'modele avec perte de masse :',1pd10.3,'Msol/an')	  	
	 else
	  write(6,211)abs(mdot)
	  write(2,211)abs(mdot)	
211	  format(1x,'modele avec gain de masse :',1pd10.3,'Msol/an')
	 endif	  	
	
	 write(6,*)' '
	 write(6,*)'------ Integration du modele quasi statique --------- '
	 write(2,*)' '
	 write(2,*)'------ Integration du modele quasi statique --------- '
	 if(der_num)then
	  write(6,*)' '
	  write(6,*)'Calcul du Jacobien par derivees numeriques'
	  write(2,*)' '
	  write(2,*)'Calcul du Jacobien par derivees numeriques'
	 else
	  write(6,*)' '
	  write(6,*)'Calcul du Jacobien par derivees analytiques'
	  write(2,*)' '
	  write(2,*)'Calcul du Jacobien par derivees analytiques'
	 endif
	
	 nrm=.15	
	 write(6,213)nrm
	 write(2,213)nrm	
213	 format(/,'coefficient pour Newton-Raphson modifie, nrm=',1pd10.3,/)
 
	 dpsi=psi0*0.02
	 n_min=150	!nombre de couches minimal
	 write(6,214)psi0
	 write(2,214)psi0
214	 format(/,'acc. de la fonction de repar. par couche =',1pd10.3,/)
 
	 write(6,*)' '
	 write(6,*)'Utilisation de la methode de collocation'
	 write(2,*)' '
	 write(2,*)'Utilisation de la methode de collocation'
 
c	 constantes de repartition
 
	 ctep= -1.d0
	 cter=  0.d0
	 !ctel= -1.d-25 ! Première valeur qui marche !!!!
	 !ctel= -1.d-24 !: converge plus
	 !ctel= -1.d-26  !: résultats identiques avec la situation SANS dépendance en L
	 !ctel= -5.d-25  !: meilleur valeur que -1.d-25 ; maid arrivé au démarrage de l'He cela se casse la gueule
         !ctel= -1.d-25
	 ctel= 0.d0
	 ctem= 15.d0
	 
	
c	 nombre maximum d'iterations pour la convergence du modele initial
 
	 iter_max0=40
 
	 write(6,*)' '
	 write(2,*)' '
	 write(6,*)'---------------------------------------------------'
	 write(2,*)'---------------------------------------------------'
 
	 ne=6		!nombre d'inconnues
 
	 reprise=un23 .eq. 1	!.true. si poursuite d'une evolution
	
c-----------------------------------------------------------------------
 
	 if(un23 .eq. 1)then
	  write(2,3)modelep
	  write(6,3)modelep
3	  format(//,5x,'On poursuit l''evolution du modele :',/,1x,a31)
 
	  if(rot_solid .neqv. rot_solid_p)then
	   print*,'on doit conserver le meme mode de rotation'
	   print*,'mode du modele repris : rot_solid=',rot_solid_p
	   print*,'mode dans *****.don : rot_solid=',rot_solid	
	   stop
	  elseif(w_rot .ne. w_rot_p)then
	   print*,'le modele repris a une rotation solide differente'
	   print*,'de celle lue dans CESAM_3.don'
	   write(6,25)w_rot_p,w_rot
25	   format(1x,'valeur precedente:',1pd10.3,' nouvelle valeur',1pd10.3)
	   print*,'poursuit-on avec la nouvelle valeur ? entrer o/n'
	   read*,oui
	   if(oui .ne. 'o')then
	    print*,'poursuit-on avec l''ancienne valeur ? entrer o/n'
	    read*,oui
	    if(oui .eq. 'o')then
	     w_rot=w_rot_p
	    else
	     print*,'on arrete'
	     stop
	    endif
	   endif
	  endif	
	  if(z_cte .neqv. z_cte_p)then
	   print*,'on doit conserver le meme type de calcul pour Z'
	   print*,'type de calcul du modele repris : z_cte=',z_cte_p
	   print*,'type de calcul dans *****.don : z_cte=',z_cte	
	   stop
	  endif
 
	  if(diffusion1 .neqv. diffusion)then
	   write(6,*)' '
	   write(2,*)' '
	   if(diffusion)then
	    write(6,*)'le modele repris est SANS diffusion'
	    write(6,*)'on poursuiton l''evolution AVEC diffusion'
	   else
	    write(6,*)'le modele repris est AVEC diffusion'
	    write(6,*)'on poursuit l''evolution SANS diffusion'	
	   endif
	  endif
 
	  bid=max(abs(y0-y01)/y0,abs(mtot-mtot1)/mtot,abs(alpha-alpha1)/alpha,
     1	abs(x0-x01)/x0)
	  if(bid .gt. 1.d-5 .or. nbelem .ne. nbelem1)then
	   write(6,9)y0,nbelem,mtot,x0,alpha,y01,nbelem1,mtot1,x01,alpha1
9	   format(1x,'ATTENTION, ATTENTION, ATTENTION, ATTENTION,',/,1x,
     1	'le modele repris ne correspond pas a celui a poursuivre',
     2	/,6x,'pour les donnees Z0=',1pd10.3,' nbelem=',i3,' mtot=',
     3	1pd10.3,/,1x,' X0=',1pd22.15,' alpha=',1pd22.15,/,
     4	1x,'pour le modele repris Z0=',1pd10.3,' nbelem=',i3,
     5	' mtot=',1pd10.3,/,1x,' X0=',1pd22.15,' alpha=',1pd22.15)
 
	   write(6,*)'poursuit-on le calcul avec les nouvelles donnees ? o/n'
	   read(5,'(a)')oui
	   if(oui .ne. 'o')stop
	  elseif(m_qs .ne. m_qs1 .or. m_ch1 .ne. m_ch)then
	   write(6,*)'conservation de l''ordre des splines du modele repris'
	   write(2,*)'conservation de l''ordre des splines du modele repris'
	   m_qs=m_qs1
	   m_ch=m_ch1		!le nombre de couches sera ajuste dans repartit
	  endif
	  dts=dt			!initialisation
	  mpr=m_qs+1
	  dim=(n-1)*m_qs+1
	  nl=ne*dim		!nombre de lignes
	  bloc=ne*mpr	
 
c--------------------initialisation-----------------------
	
	 elseif(abs(un23) .ge. 2)then		!modele en binaire ou ASCII
	
c	  modele d'initialisation en binaire
	
	  if(abs(un23) .eq. 2)then
	   write(2,4)modelep
	   write(6,4)modelep
4	   format(//,5x,'On utilise comme modele initial le modele :',
     1	/,1x,a50)
	   write(6,12)	y0 ,nbelem ,x0 ,mtot ,alpha ,
     1		y01,nbelem1,x01,mtot1,alpha1,n
	   write(2,12)	y0 ,nbelem ,x0 ,mtot ,alpha ,
     1		y01,nbelem1,x01,mtot1,alpha1,n
12	   format(1x,'parametres du modele a calculer',/,1x,'Y0=',1pd10.3,
     1	' nbelem=',i3,' X0=',1pd10.3,' mtot=',1pd10.3,' alpha=',
     2	1pd10.3,/,
     3	1x,'parametres du modele repris',/,1x,'Y0=',1pd10.3,
     4	' nbelem=',i3,' X0=',1pd10.3,' mtot=',1pd10.3,
     5	' alpha=',1pd10.3,' nb. de couches=',i4)
 
c	   print*,ne,n,mpr1,knot
c	   pause'apres lecture'	
	
c	   modele d'initialisation en ASCII, les masses sont en 1-M/Mtot	
 
	  elseif(abs(un23) .eq. 3)then
	   write(2,16)
16	   format(//,1x,'on utilise le modele d''initialisation en ASCII')
	   n=0
14	   n=n+1	
	   read(4,2001,end=15)bp(ne*(n-1)+5),(bp(ne*(n-1)+j),j=1,4)
2001	   format((1x,1pd17.10,1p4d12.5))
	   bp(ne*(n-1)+5)=(1.d0-bp(ne*(n-1)+5))*mtot
c	   write(6,2000)(bp(ne*(n-1)+j),j=1,5)	
	   goto14
15	   n=n-1
	   write(6,*)' '
	   write(6,*)'on utilise le modele de FOR004, nombre de couches',n
	   write(6,*)' '
c	   pause'lecture'
	
	   do i=1,n/2	!retournement 1 au centre, n a la surface
	    do j=1,5
	     bid=bp(ne*(i-1)+j)
	     bp(ne*(i-1)+j)=bp(ne*(n-i)+j)
	     bp(ne*(n-i)+j)=bid
	    enddo
	   enddo
	
c	   on passe en variables d'integration
c	   m--> m**2/3-1, p--> ln p, t--> ln t, r--> r**2, l-->l**2/3
	
	   do i=1,n
	    bp(ne*(i-1)+1)=log(bp(ne*(i-1)+1))	!pour p
	    bp(ne*(i-1)+2)=log(bp(ne*(i-1)+2))	!pour t	
	    bp(ne*(i-1)+3)=bp(ne*(i-1)+3)**2	!pour r
	    bp(ne*(i-1)+4)=bp(ne*(i-1)+4)**(2.d0/3.d0)	!pour l	
	    bp(ne*(i-1)+5)=bp(ne*(i-1)+5)**(2.d0/3.d0)	!pour m	    	    	
	    bp(ne*(i-1)+6)=1.d0				!pour psi (fictif)
	    q(i)=i
	   enddo
	   mpr1=2		!initialisation sur snoein
	   call sbsp1dn(ne,bp,q,qt,n,mpr1,knot,.false.,q(1),lx,f,dfdr)
	  endif		!un23=2 et 3
 
c	  on fixe le nombre de couches du modele d'initialisation a ncouche=300
c	  300 etant un bon ordre de grandeur pour un modele solaire
c	  ce nombre sera ensuite ajuste de facon a ce que dQ/dq~psi0	
c	  on determine les ncouche points de masse m(ncouche) pour
c	  assurer une repartition approximativement uniforme en
c	  ctep lnP+ ctem mu
 
	  do i=1,n
	   call sbsp1dn(ne,bp,q,qt,n,mpr1,knot,.true.,q(i),lx,f,dfdr)
	   esp(i)=ctep*f(1)+cter*f(3)+ctel*f(4)+ctem*f(5)
	   !print*, 'f(1)=',f(1),'f(4)=', f(4),'f(5)=', f(5)
c	   write(6,2000)q(i),(f(j),j=1,ne),esp(i)
	   if(i .gt. 1)then
	    do while(esp(i) .lt. esp(i-1))
	     esp(i)=esp(i)+0.01
	    enddo
	   endif
	  enddo
c	  write(6,*)'esp'
c	  write(6,2000)(esp(i),i=1,n)
c	
c	  on disposee les ncouche pour assurer une repartition
c	  approximativement uniforme de la fonction de repartition
 
	  ncouche=300
	  call zone_3(n,q,esp,ncouche,q_t)	!choix des nouveaux q
	
c	  write(6,*)'anciens q',n
c	  write(6,2000)(q(i),i=1,n)
c	  write(6,*)'esp'
c	  write(6,2000)(esp(i),i=1,n)
c	  write(6,*)'nouveaux q'
c	  write(6,2000)(q_t(i),i=1,ncouche)
c	  pause'apres zone'
	
c	  on se place aux ncouche points q_t, spline dans bp_t
 
	  do i=1,ncouche
	   call sbsp1dn(ne,bp,q,qt,n,mpr1,knot,.true.,q_t(i),lx,f,dfdr)
	   do j=1,ne-1		!sauf pour espacement
	    bp_t(ne*(i-1)+j)=f(j)
	   enddo
	   esp(i)=ctep*f(1)+cter*f(3)+ctel*f(4)+ctem*f(5)	
	   if(i .gt. 1)then
	    do while(esp(i) .lt. esp(i-1))
	     esp(i)=esp(i)+0.01
	    enddo
	   endif
c	   write(6,2000)(f(j),j=1,ne-1),esp(i)
	  enddo
c	  pause'nouvel esp'
 
c	  bp_t dans la base de snoein d'ordre mpr1 sur ncouche
 
	  n_t=ncouche
	  n=ncouche
	  do i=1,n_t
	   q_t(i)=i
	   q(i)=i
	  enddo
	
c	  pour la fonction d'espacement aux nouveaux points
 
	  call sbsp1dn(1,esp,q_t,qt_t,n_t,mpr1,knot_t,.false.,
     1	q_t(1),lx,f,dfdr)
	
c	  derivee de la fonction d'espacement
 
	  do i=1,n_t
	   call sbsp1dn(1,esp,q_t,qt_t,n_t,mpr1,knot_t,.true.,
     1	q_t(i),lx,f,dfdr)
	   bp_t(ne*(i-1)+ne)=dfdr(1)
	  enddo
 
c	  la nouvelle spline	
 
	  call sbsp1dn(ne,bp_t,q_t,qt_t,n_t,mpr1,knot_t,.false.,
     1	q_t(1),lx,f,dfdr)
	
c	  do i=1,n_t
c	   call sbsp1dn(ne,bp_t,q_t,qt_t,n_t,mpr1,knot_t,.true.,
c	1	q_t(i),lx,f,dfdr)
c	   write(6,2000)(f(j),j=1,ne)
c	  enddo
c	  pause'les bp_t verification'
 
c	  base pour la collocation aux nouveaux points
 
	  call noedif(q,qt,n,m_qs,1,knot)
	
c	  changement de base
	
	  call newspn(ne,q_t,qt_t,knot_t,mpr1,qt,knot,m_qs+1,bp_t,bp,a,indpc)
	
c	  nouvelle distribution
 
	  mpr=m_qs+1		!ordre effectif de la nouvelle spline
	  dim=(n-1)*m_qs+1	!dimension
	  nl=ne*dim		!nombre de lignes
	  bloc=ne*mpr		!longueur d'un bloc
	  do i=3,5		!evite des erreurs d'arrondi
	   bp(i)=0.d0
	  enddo
c	  do i=1,n
c	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),lx,f,dfdr)
c	   write(6,2000)(f(j),j=1,ne)
c	  enddo
c	  print*,ne,knot_t,knot,mpr
c	  pause'nouvelle base'	
 
c---------pour la composition chimique et la rotation---------
 
c	  interpolation en mu=(m/Msol)**2/3
 
c	  initialisation totale pour les vecteurs
c	  de composition chimique generalise elements chimiques + MA + Z
 
	  do i=1,n
	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),k,f,dfdr)
c	   write(6,2000)(f(j),j=1,ne)
	   f(5)=max(f(5),0.d0)
	   mc(i)=f(5)
	   do j=1,nchim
	    chim_t(nbelem*(i-1)+j)=xchim0(j)
	    chim  (nbelem*(i-1)+j)=xchim0(j)
	   enddo	!j
	   if(iw .gt. 1)then	!moment angulaire w*(r/Rsol)**2
	    chim_t(nbelem*(i-1)+iw)=w_rot*f(3)
	    chim  (nbelem*(i-1)+iw)=w_rot*f(3)
	   endif
	   if(iz .gt. 1)then	!Z
	    chim_t(nbelem*(i-1)+iz)=z0
	    chim  (nbelem*(i-1)+iz)=z0
	   endif
	  enddo
	
c	  on initialise le reste du vecteur de compostion chimique generalise
c	  pour eviter un changement de base quand il y a diffusion	
	
	  do i=n+1,pchimt
	   do j=1,nchim
	    chim_t(nbelem*(i-1)+j)=xchim0(j)
	    chim  (nbelem*(i-1)+j)=xchim0(j)
	   enddo	!j
	   if(iw .gt. 1)then	!moment angulaire w*(r/Rsol)**2
	    chim_t(nbelem*(i-1)+iw)=w_rot*f(3)
	    chim  (nbelem*(i-1)+iw)=w_rot*f(3)
	   endif
	   if(iz .gt. 1)then	!Z
	    chim_t(nbelem*(i-1)+iz)=z0
	    chim  (nbelem*(i-1)+iz)=z0
	   endif
	  enddo
 
	  nc=n	!composition chimique
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.false.,mc(1),lx,xchim,dxchim)
	
c	  write(6,*)'interpolation de la comp. chim., nc, nbelem',
c	1	nc,nbelem,nchim,(nom_elem(i),i=1,nbelem)
c	  do i=1,nc
c	   call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
c	1	knotc,.true.,mc(i),lx,xchim,dxchim)
c	   write(6,2000)mc(i),(xchim(j),j=1,min(7,nbelem))
c	  enddo
c	  pause'resout_3 xchim'
c	  write(6,*)'mc'
c	  write(6,2000)(mc(i),i=1,nc)
c	  write(6,*)'chim_t'
c	  do j=1,nbelem
c	   write(6,2000)(chim_t(nbelem*(i-1)+j),i=1,nc)
c	  enddo	!j
c	  write(6,*)'chim'
c	  do j=1,nbelem
c	   write(6,2000)(chim(nbelem*(i-1)+j),i=1,nc)
c	  enddo	!j
	
	  dts=dt		 	!initialisation
	  dt=0.
	 else
	  write(6,*)'dans resout_3, erreur sur un23 =',un23
	  stop
	 endif	!suivant un23
 
c----------------------------------------------------------------------------
 
c	 write(6,*)'appel dans resout_3 avant lim_zc_3,n',n,knot,nc,knotc
c	 write(6,2000)(q(i),i=1,n)
c	 write(6,2000)(qt(i),i=1,knot)
c	 write(6,2000)(mc(i),i=1,nc)
c	 write(6,2000)(mct(i),i=1,knotc)	
c	 pause'avant un23'
 
	 if(un23 .ne. 1)then
	  mstar=mtot
	 else
	  ncouche=n	!on garde le nombre de couches du modele a poursuivre
	 endif
 
	 call lim_zc_3(bp,q,qt,knot,n,jlim,lim,lconv,derxx,lderxx,
     1	xx,xl,ncouche,ni,nl,fac,
     2	mc,mct,nc,knotc,chim,m_zc,mstar,r2,m23,dim,
     3	r_zc,r_ov,bloc,etat,opa,conv,nuc)
	
	 ovsht=max(ovshts,ovshti) .gt. 0.
 
	 if(un23 .eq. 1 .and. age .ne. 0.d0)then !extraction r2, m23
	  do i=1,n		!extraction de r2, m23 pour inter_3
	   r2(i)=bp(ne*(i-1)*m_qs+3)
	   m23(i)=bp(ne*(i-1)*m_qs+5)
	  enddo		!i
	 endif
	
c	 si un23<2 reprise d'une evolution ou modele d'age 0 non en equilibre,
c	 par exemple pre-sequence principale: il y a extraction de
c	 p,t,r,l,m pour ecriture puis retour vers cesam pour ecritures puis
c	 retour pour la reprise de l'evolution
 
	 if(un23 .lt. 2)return
	endif		!initialisation
 
c"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 
c	write(6,*)'resout_3 dt',dt
 
c	on modifie le nb. de couche pour que psi ~ psi0
 
c	write(6,2000)bp(6),psi0
c	pause'psi0'
	if(abs(bp(6)-psi0) .gt. dpsi)then
	 new_n=bp(6)/psi0*n
	 new_n=max(n_min,min(new_n,pn))
	 write(6,218)n,new_n,bp(6),psi0
	 write(2,218)n,new_n,bp(6),psi0	
218	 format(/,1x,'modif. nb. couches pour var. quasi-stat. de',i4,' a',i4,
     1	' psi=',1pd10.3,', psi0=',1pd10.3,/)
	
	else
	 new_n=n
	 write(6,219)n,bp(6),psi0
	 write(2,219)n,bp(6),psi0	
219	 format(/,1x,'no modif. nb. couches pour var. quasi-stat., n=',i4,
     1	' psi=',1pd10.3,', psi0=',1pd10.3,/)	
	endif
 
c	pas suivant
 
c	print*,'reprise',mstar,mstar_t,dt,reprise	
 
	if(dt .gt. 0.d0)then
	 if(reprise)then
	  dts=dt
	  reprise=.false.
	
c	  avant update on calcule la variation de masse entre t+dt et t
	
c	  print*,'dans resout_3 reprise',mstar,mstar_t
	  bid=mstar_t
	
	  call perte_m(bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     1	bid,dt_,age-dt,old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     2	etat,nuc)
	 endif
	
c	 write(6,2000)(q(i),i=1,n)
c	 pause'avant update_3'
	 call update_3(.true.,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,bp__,q__,n__,qt__,knot__,chim__,mc__,
     6	nc__,mct__,knotc__,r2__,m23__,dt__,
     7	old_m23,new_m23,nm,new_m23t,knotm,
     8	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     9	old_m23__,new_m23__,nm__,new_m23t__,knotm__)	
c	 write(6,2000)(q(i),i=1,n)
 
c	 perte et defaut de masse
c	 dans update_3 on vient d'assurer mstar=mstar_t
	
c	 print*,'resout_3,2',mstar,mstar_t
 
c	 l'atmosphere au temps t en fraction de mstar au temps t
 
	 call perte_m(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,mstar,
     1	dt,age,old_m23,new_m23,nm,new_m23t,knotm,etat,nuc)
 
c	 pause'avant lim_zc_3'
	 call lim_zc_3(bp,q,qt,knot,n,jlim,lim,lconv,derxx,lderxx,
     1	xx,xl,new_n,ni,nl,fac,
     2	mc,mct,nc,knotc,chim,m_zc,mstar,r2,m23,dim,
     3	r_zc,r_ov,bloc,etat,opa,conv,nuc)
c	 pause'apres lim_zc_3'
	
	 call evol_3(.true.,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,r2_,m23_,dt_,
     4	bp__,q__,n__,qt__,knot__,chim__,mc__,nc__,mct__,knotc__,
     5	r2__,m23__,dt__,kmax,estim,g_max,
     6	old_m23,new_m23,nm,new_m23t,knotm,
     7	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     8	old_m23__,new_m23__,nm__,new_m23t__,knotm__,
     9	etat,opa,conv,nuc,coeff_diff)	
	 new=.true.		!utilisation des TdS du pas precedent
	
	else			!dt = 0
	 mstar=mtot
	 g_max=.false.
	 new=.false.
	
	 call update_3(.true.,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,bp__,q__,n__,qt__,knot__,chim__,mc__,
     6	nc__,mct__,knotc__,r2__,m23__,dt__,
     7	old_m23,new_m23,nm,new_m23t,knotm,
     8	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     9	old_m23__,new_m23__,nm__,new_m23t__,knotm__)
	
	 call lim_zc_3(bp,q,qt,knot,n,jlim,lim,lconv,derxx,lderxx,
     1	xx,xl,new_n,ni,nl,fac,
     2	mc,mct,nc,knotc,chim,m_zc,mstar,r2,m23,dim,
     3	r_zc,r_ov,bloc,etat,opa,conv,nuc)
	 dt=0.
	endif			!dt > 0
 
c	reprise du modele au temps precedent en cas de difficulte
c	ie. TdS trop grand ou non CV soit dans temporel soit dans Newton-Raphson
 
19	do while(g_max)
c dty
	 if ( .NOT. dty ) then
	 dt=dt/2.
	 else
	    if ( aug_dt_y ) then
	       dt=dt*1.10d0
	    else
	       dt=dt/2.
	    end if
	 end if
 
         write(6,*) ' '
         write(6,*) 'Valeur du pas temporel : ', dt
         write(6,*) ' '
 
	 if(dt .lt. dtmin)then
	    if(repri)then
	       write(6,*) ' '
	       write(6,*) ' '
	       write(6,*) 'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
	       write(6,*) ' '
	       write(6,*) 'dt= ', dt , ' < dt min.= ',
     +                    dtmin, 'reprise automatique, avec :'
	       write(6,*) ' '
	       write(6,*) 'dt = dtmax = ', dtmax
	       write(6,*) ' '
	       write(6,*) 'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
	       write(6,*) ' '

	       dt=dtmax
	    else
	       write(6,*)' '
	       write(6,*)'dt= ', dt , ' < dt min.= ',
     +                    dtmin, ' on abandonne!'
	       !pause 'ABANDON'
	       write(6,*)'ABANDON'
	       stop
	    end if
	 endif
 
	 call update_3(.false.,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,bp__,q__,n__,qt__,knot__,chim__,mc__,
     6	nc__,mct__,knotc__,r2__,m23__,dt__,
     7	old_m23,new_m23,nm,new_m23t,knotm,
     8	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     9	old_m23__,new_m23__,nm__,new_m23t__,knotm__)
	
	 call lim_zc_3(bp,q,qt,knot,n,jlim,lim,lconv,derxx,lderxx,
     1	xx,xl,new_n,ni,nl,fac,
     2	mc,mct,nc,knotc,chim,m_zc,mstar,r2,m23,dim,
     3	r_zc,r_ov,bloc,etat,opa,conv,nuc)
c	 write(6,*)'RESOUT_3 (1) dt=',dt
 
	 call evol_3(.true.,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,r2_,m23_,dt_,
     4	bp__,q__,n__,qt__,knot__,chim__,mc__,nc__,mct__,knotc__,
     5	r2__,m23__,dt__,kmax,estim,g_max,
     6	old_m23,new_m23,nm,new_m23t,knotm,
     7	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     8	old_m23__,new_m23__,nm__,new_m23t__,knotm__,
     9	etat,opa,conv,nuc,coeff_diff)	
 
	 new=.true.		!utilisation des TdS du pas precedent
	enddo		!g_max
 
	compt=0
	err=1.d3
 
c	write(6,*)'n,mpr',n,mpr
 
c--------------------iterations Newton Raphson----------------------
 
	compt=0
20	li=0				!pour le decompte des limites
	ligne=0
	do ip=1,n-1
	 spi=(ip-1)*m_qs			!indice-1 de la 1-iere B-spline
	 indcol=ne*spi+1
	 do i=1,m_qs
	  cx=(n-1)*(i-1)+ip
 
c	  calcul des coefficients des equations et des seconds membres
c	  au point xx, a cause de la non linearite il faut
c	  interpoler la solution provisoire pour obtenir la fonction et ses
c	  derivees y(variable,derivee) au point xx ( ou xl(li) )
 
c	  ae(ne*(ne*(der-1)+var-1)+eq)=coefficient de la (der-1)-ieme
c	  derivee de la var-ieme variable de la eq-ieme equation
c	  y(variable,derivee)=y(ne*(derivee-1)+variable)
 
	  do var=1,ne
	   do der=1,rp1
	    y(ne*(der-1)+var)=0.
	    do spl=1,mpr
	     y(ne*(der-1)+var)=y(ne*(der-1)+var)+
     1    derxx(rp1*(mpr*(spi+i-1)+spl-1)+der)*bp(ne*(spi+spl-1)+var)
	    enddo	!spl
	   enddo	!der
	  enddo	!var
	  do der=1,ne*ne*rp1	!mise a 0 des derivees
	   ae(der)=0.
	  enddo	!der
c	  write(6,*)'mpr,m_qs,ne',mpr,m_qs,ne
c	  write(6,2000)(derxx(j),j=1,8)
c	  write(6,2000)(bp(j),j=1,ne)
c	  print*,'static_m_3 1',cx,li
c	  write(6,2000)(y(j),j=1,ne),xx(cx)
c	  pause'avant static_m_3'
	  call static_m_3(1,xx,cx,li,y,be,ae,g_max,xl,compt,new,fac(ip),
     1	mc,mct,nc,knotc,chim,r_zc,r_ov,lim,mstar,age,dt,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	m23_t,r2_t,m23t_t,knot23_t,tdst,
     4	old_m23,new_m23,nm,new_m23t,knotm,	
     5	etat,opa,conv,nuc,lim_ext,tdetau)
		
	  if(g_max)then		!test sur TdS
	   write(6,*)'TdS varie trop, couche:',ip,' reinitialisation'
	   goto19
	  endif
 
c	  write(6,*)'indices colloc, couche',cx,ip
c	  do eq=1,ne
c	   write(6,*)'equation',eq
c	   do der=1,rp1
c	    write(6,2000)(ae(ne*(ne*(der-1)+var-1)+eq),var=1,ne)
c	   enddo
c	  enddo
 
	  do eq=1,ne		!eq-ieme equation
	   ligne=ligne+1
	   indpc(ligne)=indcol
	   b(ligne)=be(eq)
 
	   do spl=1,mpr	!spl - ieme spline
	    do var=1,ne
	     indice=nl*(ne*(spl-1)+var-1)+ligne
	     a(indice)=0.
	     do der=1,rp1 !der-1-eme derivee
	      a(indice)=a(indice)+ae(ne*(ne*(der-1)+var-1)+eq)*
     1		derxx(rp1*(mpr*(spi+i-1)+spl-1)+der)
	     enddo	!der
	    enddo	!var
	   enddo	!spl
 
c	   verification que la ligne n'est pas = 0	
	
	   bid=0.
	   do spl=1,mpr
	    do var=1,ne
	     indice=nl*(ne*(spl-1)+var-1)+ligne
	     bid=max(bid,abs(a(indice)))
	    enddo
	   enddo
	   if(bid .le. 0.)then
	    print*,'ligne',ligne,' nulle, equation',eq,' cx',cx
	    write(6,2000)xx(cx)
	    write(6,2000)(y(var),var=1,ne)
	    write(6,2000)(y(var),var=1+ne,2*ne)
	    g_max=.true.
	    goto 19
c	    pause'ligne nulle'
	   endif
	  enddo		!eq
	 enddo		!i
 
	 if(ni(ip) .eq. 0)goto 30		!limites
	 do i=1,ni(ip)
	  li=li+1
	  do var=1,ne
	   do der=1,rp1
	    y(ne*(der-1)+var)=0.
	    do spl=1,mpr
	     y(ne*(der-1)+var)=y(ne*(der-1)+var)+
     1    lderxx(rp1*(mpr*(li-1)+spl-1)+der)*bp(ne*(spi+spl-1)+var)
	    enddo	!spl
	   enddo	!der
	  enddo		!var
	  do der=1,ne*ne*rp1	!mise a 0 des derivees
	   ae(der)=0.
	  enddo	!der
c	  print*,'static_m_3 2',cx,li
c	  write(6,2000)(lderxx(j),j=1,8)
c	  write(6,2000)(y(j),j=1,ne),xx(cx)
	
	  call static_m_3(2,xx,cx,li,y,be,ae,g_max,xl,compt,new,fac(ip),
     1	mc,mct,nc,knotc,chim,r_zc,r_ov,lim,mstar,age,dt,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	m23_t,r2_t,m23t_t,knot23_t,tdst,	
     4	old_m23,new_m23,nm,new_m23t,knotm,	
     5	etat,opa,conv,nuc,lim_ext,tdetau)
	
	  ligne=ligne+1
	  indpc(ligne)=indcol
	  b(ligne)=be(1)
 
	  do spl=1,mpr
	   do var=1,ne
	    indice=nl*(ne*(spl-1)+var-1)+ligne
	    a(indice)=0.
	    do der=1,rp1    !derivee <r pour cond.lim.
	     a(indice)=a(indice)+ae(ne*(ne*(der-1)+var-1)+1)*
     1    lderxx(rp1*(mpr*(li-1)+spl-1)+der)
	    enddo	!der
	   enddo	!var
	  enddo	!spl
	 enddo		!i
30	enddo		!ip
 
c	do ligne=1,nl,8
c	 write(6,1000)(b(ligne+i),i=0,7)
c	enddo	!ip
c	write(6,*)' '
 
c	do spl=1,dim
c	 write(6,2000)(b(ne*(spl-1)+var),var=1,ne)
c	enddo
c	print*,'nl,dim',nl,dim
c	pause'avant solution'
 
c	do ligne=1,nl
c	 write(6,*)'ligne, couche',ligne,ligne/mpr/ne
c	 do i=1,mpr
c	  write(6,1000)(a(nl*(ne*(i-1)+var-1)+ligne),var=1,ne)
c	 enddo	!i
c	enddo	!ligne
1000	format((1x,1p10e10.3))
c	write(6,*)'resout_3; nl,n,m_qs,ne,mpr,bloc',nl,n,m_qs,ne,mpr,bloc
c	write(6,2000)(derxx(i),i=1,8)
	
	call gausdp_g(a,b,indpc,nl,ne*dim,bloc,inversible)
	if(.not. inversible)stop
	
c	do spl=1,dim
c	 write(6,2000)(b(ne*(spl-1)+var),var=1,ne)
c	enddo
c	write(6,2000)(derxx(i),i=1,8)
c	pause'apres solution'
 
c	limitation des corrections
 
	corr=1.
	do spl=1,dim
	 do var=1,ne-1		!sauf pour la repartition
	  indice=ne*(spl-1)+var
c	  write(6,*)indice,corr
	  if(abs(bp(indice)) .gt. err1)then
	   bid=nrm*abs(bp(indice))
	   do while(corr*abs(b(indice)) .gt. bid)
c	    print*,indice,spl,var
c	    write(6,2000)corr,bp(indice),bid,nrm
	    corr=corr/2.
c	    pause
	   enddo
	  endif
	 enddo	!var
	enddo	!spl
 
c	corrections et b --> bp
 
	errp=err
	err=0.
	do spl=1,dim
	 do var=1,ne
	  er=0.
	  indice=ne*(spl-1)+var
	  if(abs(bp(indice)) .gt. 1.d0)then
	   er=abs(b(indice)/bp(indice))
	  else
	   er=abs(b(indice))
	  endif
	  bp(indice)=bp(indice)-b(indice)*corr
	  err=max(er,err)	!erreur max
	  if(er .eq. err)then
	   ipe=spl
	   vare=var
	  endif
	 enddo	!var
	enddo	!spl
 
	ipe=ipe/mpr
	compt=compt+1
	write(6,*)' '
	write(6,100)compt,err,vare,ipe,1./corr
100	format(1x,'iter. globale:',i3,' err. max:',1pd8.1,
     1	', variable:',i2,', couche: ',i3,', corr:',1pd8.1)
	ipe=min(n,max(ipe,1))
	call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(ipe),lx,f,dfdr)
	write(6,101)(exp(f(i)),i=1,2),sqrt(abs(f(3))),
     1	(sqrt(abs(f(i)))**3,i=4,5)
101	format(1x,'variables P, T, R, L, M :',1p5d10.3)	
	call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lx,xchim,dxchim)
	call chim_gram_3(xchim,dxchim,nuc_m)
	write(6,102)(xchim(i),i=1,min(nbelem,6))
102	format(1x,'ab. 6 1-iers elem. :',1p6d10.3)	
 
c	pour l'interpolation inverse m23 ou r2 ----> lnp, lnt, r2, m23, l23
c	avec inter_3, extraction de r2 et m23 de la solution
 
	do i=3,5		!on evite ainsi les erreurs d'arrondi
	 bp(i)=0.d0
	enddo
	
	do i=1,n		!initialisation
	 r2(i)=bp(ne*(i-1)*m_qs+3)
	 m23(i)=bp(ne*(i-1)*m_qs+5)
	enddo
	
	if(corr .eq. 1.d0)then
	 iter_max=40
	elseif(corr .eq. 2.d0)then
	 iter_max=30
	else
	 iter_max=20
	endif
	
c	write(6,2000)bp(6)
c	pause'psi'
 
c	if(.true.)then
	if(.false.)then
	 do i=1,n			!convergence totale si dt>0
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),lx,f,dfdr)
	  p(i)=exp(f(1))	!variable ln p
	  t(i)=exp(f(2))	!variable ln T
	  r(i)=sqrt(abs(f(3)))	!rayon/Rsol
	  l(i)=sqrt(abs(f(4)))**3	!l/Lsol
	  m(i)=sqrt(abs(f(5)))**3	!m/Mtot
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lx,xchim,dxchim)
c	  write(6,2000)p(i),t(i),r(i),l(i),m(i)
c	  write(6,2000)(xchim(j),j=1,min(8,nbelem))
	  write(6,2000)p(i),t(i),r(i),l(i),m(i),f(6),xchim(1),dble(i)
	 enddo
	 pause'ecriture'
	endif		!ecriture de la solution
 
c	gestion des iterations et de l'evolution temporelle
 
	logic=compt .ge. 20 .and. err .lt. 10.*precix
     1	.and. abs(err-errp) .lt. precix !stagnation
	if((compt .ge. 2 .and. err .lt. precix) .or. logic)then

c dty, modif. de la gestion du pas de temsp pour avoir Yc decroissant sur
c les Blue Loop :

	   if ( dty ) then
	      if ( chim_t(ihe4) .ge. chim(ihe4) ) then
		 g_max=.true.
		 aug_dt_y=.true. ! Indique qu'on doit augmenter dt !
		 write(6,*) 'Ajustement pas temporel (+10%) pour Yc decroissant'
		 write(6,*) ' Yc(t+dt)= ', chim_t(ihe4), ' Yc(t)= ', chim(ihe4)
		 write(6,*) 'Nouveau pas temporel a essayer : ', dt*1.1
		 go to 19
	      else
		 aug_dt_y=.false.
	      end if
	  end if


	 if(err .gt. precix)then
	  write(6,*)' '
	  write(6,*)'on force la convergence'
	 else
	  write(6,*)' '
	  write(6,*)'modele converge'
	 endif
c	 if(dt .gt. 0. .and. kmax .ne. 0)write(2,111)kmax,dt,
c	1	(estim(i),i=1,nchim)
111	 format(/,1x,'integration temporelle de la composition chimique',
     1	' erreur maximale couche',i4,/,1x,'pas temporel',1pd10.3,
     2	' erreur par element:',/,1x,(1p12d8.1))
 
c	 modele totalement convectif: lim=1,jlim(1)=n,lconv(1)=.false.
c	 modele totalement radiatif: lim=0,jlim(i)=-100,lconv(i)=.false.
 
	 if(lim .gt. 0)then
	  if(lim .eq. 1 .and. jlim(1) .eq. 1 .and. .not. lconv(1))then
	   write(2,*)'modele completement convectif'
	  else
	   write(2,*)'limites Zones radiatives/convectives'
	   do j=1,lim
	    if(ovsht)then
	     if(lconv(j))then	
	      write(2,122)jlim(j),m_zc(j),r_zc(j),r_ov(j)
	     else
	      write(2,222)jlim(j),m_zc(j),r_zc(j),r_ov(j)
	     endif
	    else
	     if(lconv(j))then	
	      write(2,1221)jlim(j),m_zc(j),r_zc(j)
	     else
	      write(2,2221)jlim(j),m_zc(j),r_zc(j)
	     endif
	    endif
	   enddo
1221	   format(1x,'couche:',i4,' masse:',1pd10.3,' rayon=',1pd10.3,
     1	' debut de ZC')
2221	   format(1x,'couche:',i4,' masse:',1pd10.3,' rayon=',1pd10.3,
     1	' fin de ZC')
122	   format(1x,'couche:',i4,' masse:',1pd10.3,' rayon=',1pd10.3,
     1	' rayon_ov=',1pd10.3,' debut de ZC')
222	   format(1x,'couche:',i4,' masse:',1pd10.3,' rayon=',1pd10.3,
     1	' rayon_ov=',1pd10.3,' fin de ZC')
	  endif
	 else
	  write(2,*)'modele completement radiatif'
	 endif
	 write(6,*)' '
	
c	 modele initiaux : on replace la comp. chim., Z et Mw
c			   initialisation de ***_
	
	 if(abs(un23) .gt. 1 .and. dt .eq. 0.d0)then
	  nc=n
	  do i=1,n
	   mc(i)=m23(i)
	   do j=1,nchim
	    chim(nbelem*(i-1)+j)=xchim0(j)	!initialisation pour contour
	   enddo	!j
	   if(iz .gt. 1)chim(nbelem*(i-1)+iz)=z0
	   if(iw .gt. 1)chim(nbelem*(i-1)+iw)=w_rot*bp(ne*(i-1)*m_qs+3)
	  enddo		!i
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.false.,mc(1),lx,xchim,dxchim)
	 endif
 
c	 print*,'avant return',mstar,mstar_t,mtot
 
	 return
 
	elseif(dt .gt. 0.d0)then
	 if(err .gt. 1.d3 .and. compt .gt. 1)then
	  write(6,*)'pas de convergence dans resout_3 modele trop eloigne'
	  write(6,*)'diminution du pas temporel, reinitialisation'
	  g_max=.true.
	  goto19
	 elseif(compt .ge. 8 .and. err .gt. err4 .and. age .gt. 0.d0)then
	  write(6,*)'mauvaise convergence: diminution du pas temporel'
	  g_max=.true.
	  goto19
	 elseif(compt .ge. iter_max .and. age .gt. 0.d0)then
	  write(6,*)'pas de convergence dans resout_3'
	  write(6,*)'diminution du pas temporel, reinitialisation'
	  g_max=.true.
	  goto19
	 elseif(compt .ge. iter_max0)then
	  write(6,*)'pas de convergence dans resout_3'
	  write(6,*)'diminution du pas temporel, reinitialisation'
	  g_max=.true.
	  goto19
	 endif
	elseif(compt .ge. iter_max0)then
	 write(6,*)'pas de convergence dans resout_3'
	 stop
	endif
	
	new=.false.		!nouveau TdS
 
c	au plus 4 iterations de la composition chimique
 
	ini=compt .lt. ini0
	if(ini)then
c	 write(6,2000)(derxx(i),i=1,8)
 
	 call lim_zc_3(bp,q,qt,knot,n,jlim,lim,lconv,derxx,lderxx,
     1	xx,xl,n,ni,nl,fac,
     2	mc,mct,nc,knotc,chim,m_zc,mstar,r2,m23,dim,
     3	r_zc,r_ov,bloc,etat,opa,conv,nuc)
c	 write(6,2000)(derxx(i),i=1,8)
c	 write(6,2000)(lderxx(i),i=1,8)	
c	 pause'apres lim_zc_3 en bas'
 
	 if(dt .gt. 0.)then
 
c	  reevaluation de la perte et du defaut de masse	
	
	  mstar=mstar_t
	
	  call perte_m(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,mstar,
     1	dt,age,old_m23,new_m23,nm,new_m23t,knotm,etat,nuc)
	
	  call evol_3(.false.,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,r2_,m23_,dt_,
     4	bp__,q__,n__,qt__,knot__,chim__,mc__,nc__,mct__,knotc__,
     5	r2__,m23__,dt__,kmax,estim,g_max,
     6	old_m23,new_m23,nm,new_m23t,knotm,
     7	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     8	old_m23__,new_m23__,nm__,new_m23t__,knotm__,
     9	etat,opa,conv,nuc,coeff_diff)	
 
	  if(g_max)goto19		!dt-->dt/2
	 endif
	endif
	goto20			!nouveau modele statique
 
	end
 
