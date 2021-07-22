c********************************************************************
 
	subroutine static_m_3(fait,xx,cx,li,y,be,ae,g_max,xl,compt,new,fac,
     1	mc,mct,nc,knotc,chim,r_zc,r_ov,lim,mstar,age,dt,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	m23_t,r2_t,m23t_t,knot23_t,tdst,	
     4	old_m23,new_m23,nm,new_m23t,knotm,
     5	etat,opa,conv,nuc,lim_ext,tdetau)

c -------------------------------------------------------------------------
c Septembre 2003, Varsovie, D.C.
c
c Je réintroduis la dépendance en L (luminosité of course) de la fonction
c de répartition.
c
c Voir mot clé "depl"
c -------------------------------------------------------------------------

c	Les quantites caracteristiques du probleme sont calculees dans
c	cette routine
 
c	spline collocation avec fonction de repartition
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	08 10 96 introduction de la rotation
c	08 01 96 rectification de ctem*fac
c	05 05 97 2 possibilites de calcul du TdS
 
c	CESAM3
 
c entrees
c	fait=1
c	-calcule les residus be(i) des equations au point xx i=1,ne
c	fait=2
c	-calcule le "residu" be(1) de la condition au point limite li
 
c	xx : table des points de collocation
c	cx : indice du point de collocation
c	li : numero de la limite
c	y : variables au point de collocation xx(cx)
c	age: age du modele au temps t
c	dt: pas temporel
c	xl : table des points limites
c	fac : facteur de repartition pour la couche
c	mc,mct,nc,knotc,chim :	comp. chim. a t+dt
c	n: nombre de couches a t+dt
c	compt: compteur du nb. iter. Newton Raphson
c	lnp,lnt,m23,q_t,qt_t,knot_t,nt: pour interpolation a t
c	r_zc,r_ov,lim : rayons de limites ZR/ZC et d'overshoot, nb de limites
c	mc_t,mct_t,nc_t,knotc_t,chim_t: comp. chim. a t
c	mm,mt,mr,nr,knotm,tdst: elements pour interpoler le TdS au temps t
c	mstar: masse totale au temps t+dt avec perte de masse
c	old_m23,new_m23,nm,new_m23t,knotm : interpolation m(t+dt)-->m(t)
 
c sorties
c	be : residus
c	ae : termes du jacobien pour le point de collocation
c	g_max=.true. : variation relative de U ou Ro trop forte
 
c routines externes
c	etat, opa, conv, nuc, lim_ext,tdetau
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer fait,li,cx,i,j,l,lim,knot23_t,nm,knotm,id,
     &	nc,knotc,nc_t,knotc_t,knot_t,n_t,compt
 
	real*8	be(1),y(1),xx(1),duv,ae(1),dt,chim(1),mc(1),mct(1),
     &	chim_t(1),mc_t(1),mct_t(1),r2_t(1),ay3,ay4,ay5,
     &	u_t,dup_t,dut_t,dux_t,xchim_t(pnelem),xl(1),mk_t,
     &	dxchim_t(pnelem),ro_t,drop_t,drot_t,drox_t,du_tm,dro_tm,
     &	tdst(1),r_zc(1),r_ov(1),gradconv,f(pne),dfdx(pne),
     &	bp_t(1),q_t(1),qt_t(1),m23_t(1),m23t_t(1),	
     &	old_m23(1),new_m23(1),new_m23t(1),nuc_m(pnelem),dpn,dtn
 
	real*8 	prn,tn,pext,dp_tm,dt_tm,age,lamb,tds0_t,depsm,mstar,
     &	ro,gradient,epsilon(5),text,ray,teff,rtot,teta0,dtdsm,
     &	drop,drot,ln,l13,m13,mn,nh1,nhe1,nhe2,drox,eps0,dy0,
     &	xchim(pnelem),dxchim(pnelem),u,dup,dut,dux,dpdl,dpdr,
     &	dtdl,dtdr,mext,dmdl,dmdr,d_grad,w,dwdm,dwdr,v,dvdr,dvdm,dyv,
     &	dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,p_t0,t_t0,
     &	depsp,depst,kap,dkapp,dkapt,dkapx,hp,dxchiml(pnelem),
     &	gradad,gradrad,dy1,dy2,dy3,dy4,dy5,teta,dtetap,dtetar,
     &	dtetam,psist,dpsitp,dpsitr,dpsitm,dpsitt,dum,
     &	drom,psist0,grad0,ro0,u0,ro_t0,u_t0,x00,dyvdr,dyvdm,
     &	delta,deltap,deltat,deltax,deltam,cp,dcpp,dcpt,dcpx,dcpm
 
	real*8	depsx(pnelem),fac,cte11,cte13,cte14,cte14t,cte21,
     &	p_t,t_t,tds,tds0,dtdsp,dtdst,dtdsl,
     &	tds_t,dtds_tm,drott_t,drotp_t,drotx_t,dutt_t,dutp_t,dutx_t,
     &	pext0,text0,mext0,cte15,cte25,cte16,iw2,iw2_t,diw2

      ! depl ---
      real*8 cte_1, cte_2, cte_3, cte_21, cte_22, cte_23, epsilon_n,
     &       dtetal, dpsitl
      ! --------

	logical g_max,deriv,new,init
c	data init/.true./,deriv/.false./
 
	real*8 stor,stor0,dstor,dd,unpdd,bes(pne),dern(pne)
c	data dd/0.000001d0/
 
	external etat,opa,conv,nuc,lim_ext,tdetau
 
	data init/.true./,deriv/.false./
	data dd/0.000001d0/
 
	save
 
2000	format((1x,1p8d10.3))
 
c	print*,'entree static_m_3',fait,new,cx
c	write(6,2000)dt,xx(cx),(y(i),i=1,ne)
 
	if(init)then
	 init=.false.
	 cte11=-3.d0*g/8.d0/pi*(msol/rsol**2)**2		!pour d ln P /d q
	 cte13=msol/rsol*3.d0/4.d0/pi/rsol/rsol		!pour d zeta /d q
	 cte14=msol/lsol				!pour d lambda/d q
	 cte14t=cte14/secon6				! "    " dt .ne. 0	
	 cte15=msol/rsol/4.d0/pi				!pour  la rotation
	 cte16=2.d0/3.d0/secon6*rsol**2		!pour energie cinetique

	 ! depl
	 cte_1 = -3.d0*g/8.d0/pi*(msol/rsol**2)**2
	 cte_2 = 3.d0/4.d0/pi * msol / rsol**3
	 cte_3 = msol/rsol

	 if(kipp)then
	  write(2,10)
	  write(6,10)	
10	  format(1x,/,1x,
     1	'Approximation de Kippenhahn pour le calcul du TdS',/)
	 else
	  write(6,11)	
11	  format(1x,/,1x,
     1	'Calcul exact du TdS',/)
	 endif
 
	endif
 
c	write(6,*)'fait,cx',fait,cx
 
	goto (100,200),fait
 
 
c------------------------------------------------------------------
 
c	ensemble des variables
 
c		y(1)=ln p		y(7)=(ln p)'
c		y(2)=ln t		y(8)=(ln t)'	
c		y(3)=r**2		y(9)=(r**2)'
c		y(4)=l**2/3		y(10)=(l**2/3)'
c		y(5)=m**2/3		y(11)=(m**2/3)'
c		y(6)=psi=dQ/dq(=cte)	y(12)=psi'(=0)
 
c	xx(cx) valeur de q au point de collocation q=1, 2, 3, .. , n
 
c	ae(ne(ne(der-1)+var-1)+eq)=derivee de la eq-ieme equation par
c	rapport a la (der-1)-ieme derivee de la var-ieme variable
 
100	cte21=cte11*ctep*fac	!constantes servant au calcul
	cte25=cte15*ctep*fac	!pour la rotation
	! depl
	cte_21 = cte_1 * ctep * fac
	cte_22 = cte_2 * cter * fac
	cte_23 = cte_3 * ctel * fac

	prn=exp(y(1))		!pression cgs
	tn=exp(y(2))		!temperature K
	ay3=abs(y(3))
	ray=sqrt(ay3)		!rayon
	ay4=abs(y(4))
	l13=sqrt(ay4)		!l**1/3
	ln=l13**3		!l/lsol
	ay5=abs(y(5))
	m13=sqrt(ay5)		!m**1/3
	mn=m13**3		!masse/mtot
 
c	write(6,*)'mn=',mn
 
c	si, a cause d'une erreur d'arrondi, ce qui peut arriver au voisinage
c	du centre r**2, l**2/3 ou m**23 est negatif, on tente de forcer la
c	convergence en utilisant les abs(**), les derivees sont inexactes et
c	la physique est violee... mais si ca passe les erreurs d'arrondi
c	disparaissent, parfois, au cours des iterations
 
	if(y(3)*y(5)*y(4) .le. 0.)then
	 print*
	 print*,'static_m_3, r2, l23, ou m23 < 0, on tente de CV',cx,xx(cx)
	 write(6,2000)y(3),y(5),y(4)
c	 print*,(y(i),i=1,ne)
c	 pause'erreur, static_m3, m, r ou l < 0'
	 g_max=.true.
	 return
	endif
 
c	composition chimique au temps t+dt
 
	call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(ay5,mc(nc)),l,xchim,dxchim)
	do i=1,nbelem
	 xchim(i)=abs(xchim(i))
	 dxchiml(i)=dxchim(i)*2.d0/3.d0/m13	!derivee / m pour Ledoux
	enddo
 
c	la rotation
 
	if(iw .gt. 1 .and. y(3) .ne. 0.)then		!rotation
	 w=xchim(iw)/y(3)
	 iw2=xchim(iw)*w	
	 dwdm=dxchim(iw)/y(3)	!manque la dependance en r**2
	else
	 w=w_rot
	 iw2=(y(3)*w_rot)**2
	 dwdm=0
	endif
	
c	write(6,*)'cx,nbelem/prn,tn,mn,ray / xchim / y  /dy',cx,nbelem
c	write(6,2000)prn,tn,mn,ray,xchim(1),dxchiml(1)
c	write(6,2000)(xchim(i),i=1,nbelem)
c	write(6,2000)(dxchim(i),i=1,nbelem)
c	write(6,2000)(y(i),i=1,ne)
c	write(6,2000)(y(i),i=1+ne,2*ne)
c	write(6,*)' '
 
	call thermo_3(prn,tn,xchim,mn,ln,ray,.not.der_num,dxchiml,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,lim,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp,gradad,gradrad,gradconv,d_grad,w,dwdm,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc)

      ! depl
      epsilon_n = epsilon(1) ! Energie nuc. tot. en CGS

c	write(6,*)'mn,prn,tn,(xchim(i),i=1,nbelem),epsilon,depsp,depst,depsx'
c	write(6,2000)mn,prn,tn,xchim(1),dxchiml(1),ro,
c	1	(xchim(i),i=1,nbelem),epsilon(1),depsp,depst,depsx
c	write(6,*)'couche',cx
c	write(6,2000)mn,xchim(1),prn,tn,epsilon(1),depsp,depst,depsx(1)
 
	if(.not.der_num)then
 
c	 derivees par rapport a : ln p, ln t, zeta, lambda, mu
 
	 drop=drop*prn
	 drot=drot*tn
	 drom=drox*dxchim(1)
 
	 dcpp=dcpp*prn
	 dcpt=dcpt*tn
	 dcpm=dcpx*dxchim(1)
	
	 deltap=deltap*prn
	 deltat=deltat*tn
	 deltam=deltax*dxchim(1)
	
	 dup=dup*prn
	 dut=dut*tn
	 dum=dux*dxchim(1)
 
	 dgradp=dgradp*prn
	 dgradt=dgradt*tn
	 dgradr=dgradr/2.d0/ray				!d grad / d r**2+1
	 dgradl=dgradl*3.d0/2.d0*l13			!d grad / d l**2/3+1
	 dgradm=dgradm*3.d0/2.d0*m13+dgradx*dxchim(1)	!d grad / d m**2/3+1
 
	 if(epsilon(1) .gt. 1.d-25)then
	  depsp=depsp*prn/epsilon(1)			!a 1/epsilon pres
	  depst=depst*tn/epsilon(1)
	  depsm=0.
	  do i=1,nchim
	   depsm=depsm+depsx(i)*dxchim(i)
	  enddo
	  depsm=depsm/epsilon(1)
	 else
	  depsp=0.
	  depst=0.
	  depsm=0.
	 endif
	endif		!not der_num
 
	if(dt .ne. 0.)then
	
	 if(new)then			!interpolation du TdS au temps t
	  call sbsp1dn(1,tdst,m23_t,m23t_t,n_t,m_ch,
     1	knot23_t,.true.,min(ay5,m23_t(n_t)),l,tds_t,dtds_tm)
c	  write(6,*)l,n_t,m_ch,knot23_t
 
c	  write(6,2000)tds_t,dtds_tm,mn,y(5),m23_t(n_t)
c	  pause'TdS interpole'
	 else				!TdS au temps t+dt
	  call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	min(ay5,new_m23(nm)),l,f,dfdx)
	  mk_t=f(1)		
	  call inter_3('mu',bp_t,q_t,qt_t,n_t,knot_t,mk_t,f,dfdx,r2_t,m23_t)
	  p_t=exp(f(1))
	  t_t=exp(f(2))
	  dp_tm=dfdx(1)*p_t	!d p_t / d mu au temps t
	  dt_tm=dfdx(2)*t_t	!d t_t / d mu au temps t
 
	  if(kipp)then	!approximation de Kippenhahan
	   dpn=prn-p_t
	   dtn=tn-t_t
	   if(tn .gt. t_inf)then	!controle de la variation
	    duv=max(abs(dpn/prn),abs(dtn/tn))
	   else
	    duv=0.
	   endif	
	   g_max=(duv .gt. d_grav) .and. (compt .gt. 3)
	   if(g_max)then		!TdS varie trop
	    write(6,*)	
	    write(6,*)'P ou T varie trop indice',cx
	    write(6,*)'P, T, dP, dT, duv / P_t, T_t'
	    write(6,2000)prn,tn,dpn,dtn,duv
	    write(6,2000)p_t,t_t	    	
	    write(6,*)'xchim'	
	    write(6,2000)(xchim(i)*nucleo(i),i=1,min(8,nchim))
	    write(6,*)
	    return
	   else
	    duv=0.
	   endif	!g_max
	
	   tds_t=(cp*dtn-delta/ro*dpn)/dt
	
	  else		!TdS=dU+PdV	
	   call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,
     1	knotc_t,.true.,min(mk_t,mc_t(nc_t)),l,xchim_t,dxchim_t)
	   call chim_gram_3(xchim_t,dxchim_t,nuc_m)	
c	   write(6,2000)p_t,t_t,xchim_t(1),xx(cx)
	   call etat(p_t,t_t,xchim_t,.false.,		!u_t et ro_t
     1	ro_t,drop_t,drot_t,drox_t,drott_t,drotp_t,drotx_t,
     2	u_t,dup_t,dut_t,dux_t,dutt_t,dutp_t,dutx_t,nh1,nhe1,nhe2,lamb)
	
	   if(ro .gt. ro_test .and. tn .gt. t_inf)then
	    duv=max(abs(u-u_t)/u,abs(ro-ro_t)/ro)
	   else
	    duv=0.
	   endif
	   g_max=(duv .gt. d_grav) .and. (compt .gt. 3)
	   if(g_max)then		!TdS varie trop
	    write(6,*)'trop de var. de U ou RO: cx/u,ro,tds,duv/u_t,ro_t',cx
	    write(6,2000)u,ro,tds_t,duv
	    write(6,2000)u_t,ro_t
	    write(6,*)'P, T / P_t, T_t'
	    write(6,2000)prn,tn
	    write(6,2000)p_t,t_t
	    write(6,*)'xchim/xchim_t'
	    write(6,2000)(xchim(i)*nucleo(i),i=1,min(8,nchim))
	    write(6,2000)(xchim_t(i),i=1,min(8,nchim))
	    write(6,*)' '
	    return
	   endif	!g_max
	
	   du_tm =dup_t *dp_tm+dut_t *dt_tm+dux_t *dxchim_t(1)
	   dro_tm=drop_t*dp_tm+drot_t*dt_tm+drox_t*dxchim_t(1)
	   tds_t=(u-u_t-prn/ro**2*(ro-ro_t))/dt
	
c	   write(6,*)'cx/p,t,p_t,t_t,u,u_t,ro,ro_t,tds',cx
c	   write(6,2000)prn,tn,p_t,t_t,u,u_t,ro,ro_t,tds_t
 
	  endif	
	
	  if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	   iw2_t=xchim_t(iw)**2/f(3)		!energie cinetique
	  else
	   iw2_t=f(3)*w_rot**2
	  endif
	  diw2=cte16*(iw2-iw2_t)/dt
	  diw2=0.d0
	 endif	!new
	endif	!dt .ne. 0.
	
c	write(6,2000)y(5),w,y(3)
	if(w .eq. 0.)then
	 v=0.
	 dvdr=0.d0
	 dvdm=0.d0
	else
	 v=sqrt(ay5/ay3)/prn*w**2	
	 dwdr=-w/ay3
	 dvdr=-v/2./ay3+2.*v/w*dwdr
	 dvdm= v/2./ay5+2.*v/w*dwdm	
	endif
 
	dyv=cte25*v		!derivee/lnP comme celle de dy1: dyv/dlnp=-dyv
	dyvdr=cte25*dvdr
	dyvdm=cte25*dvdm

c ----------------------------------
c	la fonction de repartition :
        ! depl

	dy1=  cte_21 * (ay5/ay3)**2 / prn
     &      - cte_22 / ro * sqrt(ay5/ay3)
     &      - cte_23 * sqrt(ay5/ay4) * epsilon_n

	dy0=dy1+dyv
	teta=dy0		!teta=ctep dksi/dmu
	if(.not.der_num)then
	   ! depl
	   dtetap = - cte_21 * (ay5/ay3)**2 / prn                    !derivee/ln p (khi)
	   dtetar = - 2.d0 * (cte_21 * (ay5/ay3)**2 / prn) / ay3
     &              + 1.d0/2.d0 * (cte_22/ro * sqrt(ay5/ay3)) / ay3  !derivee/r**2 (zeta)
           dtetal = 1.d0/2.d0 * (cte_23 * sqrt(ay5/ay4) * epsilon_n) !derivee/l23 (lambda)
     &              / ay4
           dtetam =   2.d0 * (cte_21 * (ay5/ay3)**2 / prn) /ay5
     &              - 1.d0/2.d0 * (cte_22 / ro * sqrt(ay5/ay3)) / ay5
     &              - 1.d0/2.d0 * (cte_23 * sqrt(ay5/ay4) 
     &                * epsilon_n) / ay5                             !derivee/m23 (mu)
	endif
 
	teta=teta+fac*ctem			!teta+d mu/ d mu
	psist=y(6)/teta			!psi/teta
 
	if(.not.der_num)then
	   ! depl
	   dpsitp = - psist / teta * dtetap !d (psi/teta) / d ln p
	   dpsitr = - psist / teta * dtetar !d (psi/teta) / d r**2
	   dpsitl = - psist / teta * dtetal !d (psi/teta) / d l**2/3
	   dpsitm = - psist / teta * dtetam !d (psi/teta) / d m**2/3
	endif

c ----------------------------------
c	equations :
 
	dyv=cte15*v		!derivee/lnP comme celle de dy1: dyv/dlnp=-dyv
	dyvdr=cte15*dvdr
	dyvdm=cte15*dvdm
	dy1=cte11*(ay5/ay3)**2/prn
	dy0=dy1+dyv
	be(1)=y(ne+1)-dy0*psist		!dln p/dmu=3g/8 pi (mu+1/zeta+1)**2 / p
	if(.not.der_num)then
	 ae(ne*(ne*(1-1)+1-1)+1)=dy0*(psist-dpsitp)	!derivee/ln p
	 ae(ne*(ne*(1-1)+3-1)+1)=-(dyvdr-dy1*2.d0/ay3)*psist
     &				-dy0*dpsitr	!derivee/r**2
	 ! depl
	 ae(ne*(ne*(1-1)+4-1)+1)=-dy0 * dpsitl !derivee/l**2/3
	 ae(ne*(ne*(1-1)+5-1)+1)=-(dyvdm+dy1*2.d0/ay5)*psist
     &				-dy0*dpsitm	!der. /m**2/3
	 ae(ne*(ne*(1-1)+6-1)+1)=-dy0/teta		!derivee/psi
	endif
	ae(ne*(ne*(2-1)+1-1)+1)=1.d0			!derivee/(ln p)'
 
	dy2=dy0*gradient
	be(2)=y(ne+2)-dy2*psist			!dln t/dmu=dln p/dmu grad
	if(.not.der_num)then
	 ae(ne*(ne*(1-1)+1-1)+2)=ae(ne*(ne*(1-1)+1-1)+1)*gradient
     &			-dy0*dgradp*psist	!d/dlnp
	 ae(ne*(ne*(1-1)+2-1)+2)=ae(ne*(ne*(1-1)+2-1)+1)*gradient
     &			-dy0*dgradt*psist	!derivee/ln t	
	 ae(ne*(ne*(1-1)+3-1)+2)=ae(ne*(ne*(1-1)+3-1)+1)*gradient	
     &			-dy0*dgradr*psist	!derivee/r**2
	 ! depl
	 ae(ne*(ne*(1-1)+4-1)+2)=ae(ne*(ne*(1-1)+4-1)+1)*gradient	
     &			-dy0 * dgradl * psist	
     &                  -dy0 * gradient * dpsitl !derivee/l**2/3
	 ae(ne*(ne*(1-1)+5-1)+2)=ae(ne*(ne*(1-1)+5-1)+1)*gradient	
     &			-dy0*dgradm*psist	!derivee/m**2/3		
	 ae(ne*(ne*(1-1)+6-1)+2)=ae(ne*(ne*(1-1)+6-1)+1)*gradient !derivee/psi
	endif
	ae(ne*(ne*(2-1)+2-1)+2)=1.d0			!derivee/(ln t)'
 
	dy3=cte13/ro*sqrt(ay5/ay3)
	be(3)=y(ne+3)-dy3*psist		!dr**2/dmu=-3/4pi /ro sqrt(mu1/r**2)
	if(.not.der_num)then
	 ae(ne*(ne*(1-1)+1-1)+3)=dy3*(psist*drop/ro-dpsitp)  !derivee/ln p
	 ae(ne*(ne*(1-1)+2-1)+3)=dy3*psist*drot/ro  !derivee/ln t
	 ae(ne*(ne*(1-1)+3-1)+3)=dy3*(psist/ay3/2.d0-dpsitr) !derivee/r**2
	 ! depl
	 ae(ne*(ne*(1-1)+4-1)+3)= -dy3 * dpsitl !derivee/l**2/3
	 ae(ne*(ne*(1-1)+5-1)+3)=dy3*(psist*(-.5d0/ay5+drom/ro)
     &			-dpsitm)  		!derivee/m**2/3
	 ae(ne*(ne*(1-1)+6-1)+3)=-dy3/teta		!derivee/psi
	endif
	ae(ne*(ne*(2-1)+3-1)+3)=1.d0			!derivee/(r**2)'
	
	dy4=cte14*epsilon(1)*m13/l13		!dl/dmu=epsilon m13/l13 +Tds/dt
	be(4)=y(ne+4)-dy4*psist		!contribution de epsilon
	if(.not.der_num)then
	 ae(ne*(ne*(1-1)+1-1)+4)=-dy4*(psist*depsp+dpsitp)	!d/dlnp
	 ae(ne*(ne*(1-1)+2-1)+4)=-dy4*psist*depst		!derivee/ln t
	 ae(ne*(ne*(1-1)+3-1)+4)=-dy4*dpsitr			!derivee/r**2
	 ! depl
	 ae(ne*(ne*(1-1)+4-1)+4)= dy4 * (psist/2.d0/ay4 - dpsitl) !der./l**2/3	
	 ae(ne*(ne*(1-1)+5-1)+4)=-dy4*(psist*(.5d0/ay5+depsm)+dpsitm)	!d/dmu	
	 ae(ne*(ne*(1-1)+6-1)+4)=-dy4/teta		!derivee/psi
	endif
	ae(ne*(ne*(2-1)+4-1)+4)=1.d0			!derivee/(l**2/3)'
 
	if(dt .ne. 0.)then		!contribution du Tds/dt
	 dy5=cte14t*m13/l13
	 tds=dy5*tds_t
	 if(.not.der_num)then
	  if(new)then
	   dtdsl=-tds/2.d0/ay4				!derivee/l**2/3
	   dtdsm= tds/2.d0/ay5+dy5*dtds_tm	 	!derivee/m**2/3
	   ae(ne*(ne*(1-1)+1-1)+4)=ae(ne*(ne*(1-1)+1-1)+4)
     &	+tds*dpsitp				!derivee/ ln p
	   ae(ne*(ne*(1-1)+3-1)+4)=ae(ne*(ne*(1-1)+3-1)+4)
     &	+tds*dpsitr				!derivee/r**2
	   ! depl
	   ae(ne*(ne*(1-1)+4-1)+4)=ae(ne*(ne*(1-1)+4-1)+4)
     &	+psist*dtdsl + tds * dpsitl		!derivee/ l**2/3
	   ae(ne*(ne*(1-1)+5-1)+4)=ae(ne*(ne*(1-1)+5-1)+4)
     &	+tds*dpsitm+psist*dtdsm			!derivee/ m**2/3	
	   ae(ne*(ne*(1-1)+6-1)+4)=ae(ne*(ne*(1-1)+6-1)+4)
     &	+tds/teta					!derivee/ psi
	
	  else
	   dtdsl=-tds/2./ay4				!derivee/l**2/3-1
	   dtdsm= tds/2.d0/ay5	   	
	   if(kipp)then
	    dtdsp=dy5*(dcpp*dtn-delta/ro*((deltap/delta
     &	-drop/ro)*dpn+prn))/dt				!d/lnp	
	    dtdst=dy5*(dcpt*dtn+cp*tn
     &	-delta/ro*dpn*(deltat/delta-drot/ro))/dt	!d/lnt	
	    dtdsm= dtdsm+dy5*(dcpm*dtn-cp*dt_tm 		!derivee/m**2/3
     &	-delta/ro*((deltam/delta-drom/ro)*dpn-dp_tm))/dt
	   else
	    dtdsp=dy5*(dup-prn/ro**2*((ro-ro_t)*(1.d0-2.d0*drop/ro)+drop))/dt !d/lnp
	    dtdst=dy5*(dut-prn/ro**2*drot*(1.d0-2.d0*(ro-ro_t)/ro))/dt	!d/ln t
	    dtdsm= dtdsm+dy5*(dum-du_tm-prn/ro**2*(-2.d0*drom/ro*(ro-ro_t)
     &	+drom-dro_tm))/dt
	   endif		 	!derivee/m**2/3
	   ae(ne*(ne*(1-1)+1-1)+4)=ae(ne*(ne*(1-1)+1-1)+4)
     &	+tds*dpsitp+psist*dtdsp				!derivee/ ln p
	   ae(ne*(ne*(1-1)+2-1)+4)=ae(ne*(ne*(1-1)+2-1)+4)
     &	+psist*dtdst				!derivee/ ln t
	   ae(ne*(ne*(1-1)+3-1)+4)=ae(ne*(ne*(1-1)+3-1)+4)
     &	+tds*dpsitr				!derivee/r**2
	   ! depl
	   ae(ne*(ne*(1-1)+4-1)+4)=ae(ne*(ne*(1-1)+4-1)+4)
     &	+psist*dtdsl + tds * dpsitl		!derivee/ l**2/3
	   ae(ne*(ne*(1-1)+5-1)+4)=ae(ne*(ne*(1-1)+5-1)+4)
     &	+tds*dpsitm+psist*dtdsm			!derivee/ m**2/3
	   ae(ne*(ne*(1-1)+6-1)+4)=ae(ne*(ne*(1-1)+6-1)+4)
     &	+tds/teta					!derivee/ psi
	  endif		!new
	 endif		!der_num
	 be(4)=be(4)+(tds-diw2)*psist		!epsilon-tdS+d Iw2 / dt
	
	
	endif		!dt
 
	be(5)=y(ne+5)-psist				!dmu/dq=psi/teta
	if(.not.der_num)then
	 ae(ne*(ne*(1-1)+1-1)+5)=-dpsitp		!derivee/ln p
	 ae(ne*(ne*(1-1)+3-1)+5)=-dpsitr		!derivee/r**2
         ! depl
         ae(ne*(ne*(1-1)+4-1)+5)=-dpsitl                !derivee/ l**2/3
	 ae(ne*(ne*(1-1)+5-1)+5)=-dpsitm		!derivee/m**2/3
	 ae(ne*(ne*(1-1)+6-1)+5)=-1.d0/teta		!derivee/psi
	endif
	ae(ne*(ne*(2-1)+5-1)+5)= 1.d0			!derivee/(m**2/3)'
 
	be(6)=y(ne+6)					!dpsi/dq=0
	ae(ne*(ne*(2-1)+6-1)+6)=1.d0			!derivee/psi'
 
	if(der_num)then		!calcul des derivees numeriques
 
	 if(mn .gt. .5)then
	  unpdd=1.d0-dd
	 else
	  unpdd=1.d0+dd
	 endif
 
	 do id=1,ne
	  stor0=y(id)
	  stor=stor0*unpdd
	  if(stor .eq. 0.)stor=dd
	  dstor=stor-stor0
	  y(id)=stor
	  prn=exp(y(1))		!pression cgs
	  tn=exp(y(2))		!temperature K
	  ay3=abs(y(3))
	  ray=sqrt(ay3)		!rayon
	  ay4=abs(y(4))
	  l13=sqrt(ay4)		!l**1/3
	  ln=l13**3		!l/lsol
	  ay5=abs(y(5))
	  m13=sqrt(ay5)		!m**1/3
	  mn=m13**3		!masse/mtot
 
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     &	knotc,.true.,min(ay5,mc(nc)),l,xchim,dxchim)
 
	  do j=1,nchim		!dx/dm pour Ledoux
	   xchim(j)=abs(xchim(j))
	   dxchiml(j)=dxchim(j)*2.d0/3.d0/m13
	  enddo
 
	  if(iw .gt. 1 .and. y(3) .ne. 0.)then		!rotation
	   w=xchim(iw)/y(3)
	  else
	   w=w_rot
	  endif		!dwdm inutile
	
	  call thermo_3(prn,tn,xchim,mn,ln,ray,.false.,dxchiml,mstar,
     &	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,lim,
     &	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     &	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     &	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     &	hp,gradad,gradrad,gradconv,d_grad,w,0.d0,nh1,nhe1,nhe2,
     &	etat,opa,conv,nuc)
	  dy1=cte21*(ay5/ay3)**2/prn	
 
	  v=sqrt(ay5/ay3)/prn*w**2
	  dyv=cte25*v				!rotation	
 
	  dy0=dy1+dyv
	  teta=dy0+fac*ctem
	  psist=y(6)/teta			!psi/teta
	  	
	  if(dt .ne. 0.)then
	   if(new)then			!interpolation du TdS au temps t
	    call sbsp1dn(1,tdst,m23_t,m23t_t,n_t,m_ch,
     1	knot23_t,.true.,min(ay5,m23_t(n_t)),l,tds_t,dtds_tm)
c	    write(6,*)'tds_t,dtds_tm,mn,l',tds_t,dtds_tm,mn,l
	   else				!TdS au temps t+dt
	    call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	min(ay5,new_m23(nm)),l,f,dfdx)
	    mk_t=f(1)		
	    call inter_3('mu',bp_t,q_t,qt_t,n_t,knot_t,mk_t,f,dfdx,r2_t,m23_t)
 	    p_t=exp(f(1))
	    t_t=exp(f(2))
	
	    if(kipp)then	!approximation de Kippenhahan
	     dpn=prn-p_t
	     dtn=tn-t_t
	     if(tn .gt. t_inf)then
	      duv=max(abs(dpn/prn),abs(dtn/tn))	!controle de la variation
	     else
	      duv=0.
	     endif	
	     g_max=(duv .gt. d_grav) .and. (compt .gt. 3)
	     if(g_max)then		!TdS varie trop
	      write(6,*)	
	      write(6,*)'P ou T varie trop indice',cx
	      write(6,*)'P, T, dP, dT, duv / P_t, T_t'
	      write(6,2000)prn,tn,dpn,dtn,duv
	      write(6,2000)p_t,t_t	    	
	      write(6,*)'xchim, xchim_t'	
	      write(6,2000)(xchim(i)*nucleo(i),i=1,min(8,nchim))
	      write(6,2000)(xchim_t(i),i=1,min(8,nchim))
	      write(6,*)
	      return
	     else
	      duv=0.
	     endif	!g_max
	
	     tds_t=(cp*dtn-delta/ro*dpn)/dt
	
	    else		!TdS=dU+PdV
	     call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,
     1	knotc_t,.true.,min(mk_t,mc_t(nc_t)),l,xchim_t,dxchim_t)
	     call chim_gram_3(xchim_t,dxchim_t,nuc_m)
	     call etat(p_t,t_t,xchim_t,.false.,	!.true. calcul derivees 2d
     1	ro_t,drop_t,drot_t,drox_t,drott_t,drotp_t,drotx_t,
     2	u_t,dup_t,dut_t,dux_t,dutt_t,dutp_t,dutx_t,nh1,nhe1,nhe2,lamb)
	
	     if(ro .gt. ro_test .and. tn .gt. t_inf)then
	      duv=max(abs(u-u_t)/u,abs(ro-ro_t)/ro)
	     else
	      duv=0.
	     endif
	     g_max=(duv .gt. d_grav) .and. (compt .gt. 3)
	     if(g_max)then		!TdS varie trop
	      write(6,*)'trop de var. de U ou RO: cx/u,ro,tds,duv/u_t,ro_t',cx
	      write(6,2000)u,ro,tds_t,duv
	      write(6,2000)u_t,ro_t
	      write(6,*)'P, T / P_t, T_t'
	      write(6,2000)prn,tn
	      write(6,2000)p_t,t_t
	      write(6,*)'xchim/xchim_t'
	      write(6,2000)(xchim(i)*nucleo(i),i=1,min(8,nchim))
	      write(6,2000)(xchim_t(i),i=1,min(8,nchim))
	      write(6,*)' '
	      return
	     endif	!g_max
	
	     tds_t=(u-u_t-prn/ro**2*(ro-ro_t))/dt   !+Tds
	
	    endif		!kipp
	
	    if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	     iw2_t=xchim_t(iw)**2/f(3)		!energie cinetique
	    else
	     iw2_t=f(3)*w_rot**2
	    endif		!kipp
	
	    diw2=cte16*(iw2-iw2_t)/dt
	    diw2=0
	   endif	!new
	  endif		!dt
 
	  dy1=cte11*(ay5/ay3)**2/prn
	  dyv=cte15*v				!rotation	
	  dy0=dy1+dyv
	  bes(1)=y(ne+1)-dy0*psist	
	  bes(2)=y(ne+2)-dy0*gradient*psist	!dln t/dmu=dln p/dmu grad
	  bes(3)=y(ne+3)-cte13/ro*sqrt(ay5/ay3)*psist
	  bes(4)=y(ne+4)-cte14*epsilon(1)*m13/l13*psist	!dl/dmu=epsilon m13
	  if(dt .ne. 0.)then
	   dy5=cte14t*m13/l13
	   tds=dy5*tds_t
	   bes(4)=bes(4)+(tds-diw2)*psist				!epsilon-tdS
	  endif				!dt
 
	  bes(5)=y(ne+5)-psist
	  bes(6)=y(ne+6)
 
	  do j=1,ne
	   ae(ne*(ne*(1-1)+id-1)+j)=(bes(j)-be(j))/dstor
	  enddo	!j
	  y(id)=stor0
	 enddo	!i
	endif
 
c	deriv=.true.
c	deriv=cx .le. 2
c	deriv=cx .le. 3
c	deriv=cx .eq. 1
c	1	.or. cx .eq. 50
c	2	.or. cx .eq. 90
c	3	.or. (cx .gt. 74 .and. cx .lt. 76)
c	3	.or. (cx .gt. 43 .and. cx .lt. 48)
c	4	.or. (cx .gt. 88 .and. cx .lt. 92)
c	5	.or. cx .eq. 149
c	6	.or. cx .eq. 103
c	7	.or. cx. ge. 148
c	8	.or. cx. ge. 149
c	9	.or. (cx. ge. 148 .and. cx .lt.150)
 
	if(.not. deriv)return	!tests de mise au point pour derivees
 
	if(mn .gt. .5)then
	 unpdd=1.d0-dd
	else
	 unpdd=1.d0+dd
	endif
	write(6,*)'ctep,ctem,fac,Q'
	write(6,2000)ctep,ctem,fac,fac*(y(1)*ctep+ay5*ctem)
	write(6,*)'cx,new,be,xx(cx),mstar',cx,new
	write(6,2000)(be(i),i=1,ne),xx(cx),mstar
	write(6,*)'y / dy'
	write(6,2000)(y(i),i=1,ne)
	write(6,2000)(y(i),i=1+ne,2*ne)
	write(6,*)'mn,prn,tn,ray,ln,dt,epsilon,ro'
	write(6,2000)mn,prn,tn,ray,ln,dt,epsilon(1),ro
	write(6,*)'xchim(i)'
	write(6,2000)(xchim(i),i=1,nbelem)
	write(6,*)'dxchim(i)'
	write(6,2000)(dxchim(i),i=1,nbelem)
	write(6,*)'despx(i)'
	write(6,2000)(depsx(i),i=1,nchim)
	write(6,*)'gradient,dgradp,dgradt,dgradr,dgradl,dgradm,t_t,p_t'
	write(6,2000)gradient,dgradp,dgradt,dgradr,dgradl,dgradm,t_t,p_t
	write(6,*)'psist,teta,dpsitp,dpsitr,dpsitm,dpsistteta'
	write(6,2000)psist,teta,dpsitp,dpsitt,dpsitr,dpsitm,1./teta
	write(6,*)'dtetap,dtetar,dtetam,tds_t,dp_tm,dt_tm'
	write(6,2000)dtetap,dtetar,dtetam,tds_t,dp_tm,dt_tm
	write(6,*)'depsp,depst,depsm,tds,dtdsp,dtdst,dtdsl,dtdsm'
	write(6,2000)depsp*epsilon(1),depst*epsilon(1),depsm*epsilon(1),tds,
     1	dtdsp,dtdst,dtdsl,dtdsm
	write(6,*)'u,dup,dut,dux,dum,u_t,du_tm,dtds_tm'
	write(6,2000)u,dup,dut,dux,dum,u_t,du_tm,dtds_tm
	write(6,*)'ro,drop,drot,drox,drom,ro_t,dro_tm'
	write(6,2000)ro,drop,drot,drox,drom,ro_t,dro_tm
	psist0=psist
	teta0=teta
	grad0=gradient
	eps0=epsilon(1)
	ro0=ro
	u0=u
	p_t0=p_t
	t_t0=t_t
	ro_t0=ro_t
	u_t0=u_t
	tds0=tds
	tds0_t=tds_t
	x00=xchim(1)
	do i=1,ne
	 stor0=y(i)
	 stor=stor0*unpdd
	 if(stor .eq. 0.)stor=dd
	 dstor=stor-stor0
	 y(i)=stor
	 ay3=abs(y(3))
	 ray=sqrt(ay3)		!rayon
	 ay4=abs(y(4))
	 l13=sqrt(ay4)		!l**1/3
	 ln=l13**3		!l/lsol
	 ay5=abs(y(5))
	 m13=sqrt(ay5)		!m**1/3
	 mn=m13**3		!masse/mtot
	 prn=exp(y(1))		!pression cgs
	 tn=exp(y(2))		!temperature K
c	 pause'avant spsp1dn1'
	
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
 
     1	knotc,.true.,min(ay5,mc(nc)),l,xchim,dxchim)
	 do j=1,nchim
	  xchim(j)=abs(xchim(j))
	  dxchiml(j)=dxchim(j)*2.d0/3.d0/m13		!dX/dm pour Ledoux
	 enddo
	
c	 la rotation
 
	 if(iw .gt. 1 .and. y(3) .ne. 0.)then		!rotation
	  w=xchim(iw)/y(3)
	  iw2=xchim(iw)**2/y(3)
	 else
	  w=w_rot
	  iw2=y(3)*w**2
	 endif		!dwdm inutile
	
	 call thermo_3(prn,tn,xchim,mn,ln,ray,.false.,dxchiml,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,lim,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp,gradad,gradrad,gradconv,d_grad,w,0.d0,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc)
	
	 dy1=cte21*(ay5/ay3)**2/prn
	
	 v=sqrt(ay5/ay3)/prn*w**2
	 dyv=cte25*v					!rotation	
	 dy0=dy1+dyv
	 teta=dy0+fac*ctem
	 psist=y(6)/teta			!psi/teta
	 if(dt .ne. 0.)then
c	  print*,dt
	  if(new)then			!interpolation du TdS au temps t
c	   pause'sbsp'
	   call sbsp1dn(1,tdst,m23_t,m23t_t,n_t,m_ch,
     1	knot23_t,.true.,min(ay5,m23_t(n_t)),l,tds_t,dtds_tm)
c	   write(6,*)'tds_t,dtds_tm,mn,l',tds_t,dtds_tm,mn,l
 
	  else				!TdS au temps t+dt
	   call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	min(ay5,new_m23(nm)),l,f,dfdx)
	   mk_t=f(1)	 	
	   call inter_3('mu',bp_t,q_t,qt_t,n_t,knot_t,mk_t,f,dfdx,r2_t,m23_t)
	   p_t=exp(f(1))
	   t_t=exp(f(2))
	
	   if(kipp)then	!approximation de Kippenhahan
	    dpn=prn-p_t
	    dtn=tn-t_t
	    tds_t=(cp*dtn-delta/ro*dpn)/dt
	   else		!TdS=dU+PdV
	    call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,
     1	knotc_t,.true.,min(mk_t,mc_t(nc_t)),l,xchim_t,dxchim_t)
	    call chim_gram_3(xchim_t,dxchim_t,nuc_m)
 	    p_t=exp(f(1))
	    t_t=exp(f(2))
	    call etat(p_t,t_t,xchim_t,.false.,	!.true. calcul derivees 2d
     1	ro_t,drop_t,drot_t,drox_t,drott_t,drotp_t,drotx_t,
     2	u_t,dup_t,dut_t,dux_t,dutt_t,dutp_t,dutx_t,nh1,nhe1,nhe2,lamb)
	    tds_t=(u-u_t-prn/ro**2*(ro-ro_t))/dt   !+Tds
	   endif		!kipp
	
	   if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	    iw2_t=xchim_t(iw)**2/f(3)		!energie cinetique
	   else
	    iw2_t=f(3)*w_rot**2
	   endif
	   diw2=cte16*(iw2-iw2_t)/dt
	   diw2=0
	  endif			!new
	 endif		!dt
 
	 dy1=cte11*(ay5/ay3)**2/prn
	 dyv=cte15*v				!rotation	
	 dy0=dy1+dyv
	 bes(1)=y(ne+1)-dy0*psist
	 bes(2)=y(ne+2)-dy0*gradient*psist	!dln t/dmu=dln p/dmu grad
	 dy3=cte13/ro*sqrt(ay5/ay3)
	 bes(3)=y(ne+3)-dy3*psist
	 bes(4)=y(ne+4)-cte14*epsilon(1)*m13/l13*psist	!dl/dmu=epsilon m13
	 if(dt .ne. 0.)then
	  dy5=cte14t*m13/l13
	  tds=dy5*tds_t
	  bes(4)=bes(4)+(tds-diw2)*psist		!epsilon-tdS
	 endif				!dt
 
	 bes(5)=y(ne+5)-psist
	 bes(6)=y(ne+6)
 
	 do j=1,ne
	  dern(j)=(bes(j)-be(j))/dstor
	 enddo	!j
	 y(i)=stor0
	 write(6,*)'der,psist,grad,teta,eps,ro,ro_t,u,u_t,tds,tds_t,p_t,t_t,X'
	 write(6,2000)(dern(j),ae(ne*(ne*(1-1)+i-1)+j),j=1,ne),
     1	(psist-psist0)/dstor,(gradient-grad0)/dstor,
     2	(teta-teta0)/dstor,(epsilon(1)-eps0)/dstor,
     3	(ro-ro0)/dstor,(ro_t-ro_t0)/dstor,(u-u0)/dstor,
     4	(u_t-u_t0)/dstor,(tds-tds0)/dstor,(tds_t-tds0_t)/dstor,
     5	(p_t-p_t0)/dstor,(t_t-t_t0)/dstor,(xchim(1)-x00)/dstor
	enddo	!i
 
	pause'pause entrer c/return'
 
	return
 
 
c----------------------------------------------------------------------------
 
 
c	conditions aux limites : residu be(1) au point limite li
 
 
c	li=1 : limite sur r au centre
c	li=2 : limite sur l au centre
c	li=3 : limite sur m au centre
c	li=4 : limite sur P a l'exterieur
c	li=5 : limite sur T a l'exterieur
c	li=6 : limite sur m a l'exterieur
 
 
200	if(li .eq. 1)then		!au centre
	 be(1)=y(3)			!en q=1 r**2=0
	 ae(ne*(ne*(1-1)+3-1)+1)=1.d0		!derivee/r**2
c	 write(6,*)'limite au centre',li
c	 write(6,2000)(y(i),i=1,ne),be(1)
 
	elseif(li .eq. 2)then		!au centre en l
	 be(1)=y(4)			!en q=1 l**2/3=0
	 ae(ne*(ne*(1-1)+4-1)+1)=1.d0		!derivee/l**2/3
c	 write(6,*)'limite au centre',li
c	 write(6,2000)(y(i),i=1,ne),be(1)
 
	elseif(li .eq. 3)then		!au centre en m
	 be(1)=y(5)			!en q=1 m**2/3=0
	 ae(ne*(ne*(1-1)+5-1)+1)=1.d0	!derivee/m**2/3
c	 write(6,*)'limite au centre',li
c	 write(6,2000)(y(i),i=1,ne),be(1)
c	 pause'au centre'
 
	elseif(li .eq. 4)then		!P a l'exterieur
	 mn=ay5**(3.d0/2.d0)
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(ay5,mc(nc)),l,xchim,dxchim)
	 call chim_gram_3(xchim,dxchim,nuc_m)
	 ray=sqrt(y(3))
	 l13=sqrt(y(4))
	 ln=l13**3
	 call lim_ext(.false.,ln,ray,xchim,pext,text,dpdl,dpdr,dtdl,dtdr,teff,
     1	rtot,mext,dmdl,dmdr,dy1,dy1,dy1,dy1,dy1,mstar,
     2	tdetau,etat,opa)
	 be(1)=y(1)-log(pext)	!condition sur p a l'exterieur
	 ae(ne*(ne*(1-1)+1-1)+1)=1.d0			!derivee/ln p
	 ae(ne*(ne*(1-1)+3-1)+1)=-dpdr/pext/2.d0/ray	!derivee/r**2
	 ae(ne*(ne*(1-1)+4-1)+1)=-dpdl/pext*3.d0/2.d0*l13	!derivee/l**2/3
 
c	 deriv=.true.
	 if(deriv)then		!test derivee
	  unpdd=1.d0-1.d-6
	  write(6,*)'limite',li
	  write(6,2000)(y(i),i=1,ne)
	  write(6,*)'exterieur lnp,ln pext,be(1),ray,lamb,xchim'
	  write(6,2000)y(1),log(pext),be(1),ray,ln,(xchim(i),i=1,nbelem)
	  pext0=pext
	  text0=text
	  mext0=mext
	  do i=1,ne
	   stor0=y(i)
	   stor=stor0*unpdd
	   if(stor .eq. 0.)stor=dd
	   dstor=stor-stor0
	   y(i)=stor
	   ray=sqrt(y(3))
	   ln=sqrt(y(4))**3
	   call lim_ext(.false.,ln,ray,xchim,pext,text,dpdl,dpdr,dtdl,dtdr,
     1	teff,rtot,mext,dmdl,dmdr,dy1,dy1,dy1,dy1,dy1,mstar,
     2	tdetau,etat,opa)
	   bes(1)=y(1)-log(pext)	!condition sur p a l'exterieur
	   write(6,*)log(pext),y(1),bes(1),be(1)
	   dern(1)=(bes(1)-be(1))/dstor
	   y(i)=stor0
	   write(6,*)'derivee en pext',i
	   write(6,2000)dern(1),ae(ne*(ne*(1-1)+i-1)+1),bes(1)-be(1)
	  enddo	!i
	  pext=pext0
	  text=text0
	  mext=mext0
	 endif		!deriv
 
	elseif(li .eq. 5)then		!T a l'exterieur
	 be(1)=y(2)-log(text)	!condition sur T a l'exterieur
	 ae(ne*(ne*(1-1)+2-1)+1)=1.d0			!derivee/ln t
	 ae(ne*(ne*(1-1)+3-1)+1)=-dtdr/text/2./ray	!derivee/r**2
	 ae(ne*(ne*(1-1)+4-1)+1)=-dtdl/text*3.d0/2.d0*l13	!derivee/l**2/3
	 if(deriv)then 	!test de derivation
	  write(6,*)'limite y(2),ln text',li
	  write(6,2000)(y(i),i=1,ne)
	  write(6,*)'y(2),log(text),be(1)'
	  write(6,2000)y(2),log(text),be(1)
	  do i=1,ne
	   stor0=y(i)
	   stor=stor0*unpdd
	   if(stor .eq. 0.)stor=dd
	   dstor=stor-stor0
	   y(i)=stor
	   ray=sqrt(y(3))
	   ln=sqrt(y(4))**3
	   call lim_ext(.false.,ln,ray,xchim,pext,text,dpdl,dpdr,dtdl,dtdr,
     1	teff,rtot,mext,dmdl,dmdr,dy1,dy1,dy1,dy1,dy1,mstar,
     2	tdetau,etat,opa)
	   bes(1)=y(2)-log(text)	!condition sur T a l'exterieur
	   dern(1)=(bes(1)-be(1))/dstor
	   y(i)=stor0
	   write(6,*)'derivee pour li=2',i
	   write(6,2000)dern(1),ae(ne*(ne*(1-1)+i-1)+1)
	  enddo	!i
	  pext=pext0
	  text=text0
	  mext=mext0
	 endif		!deriv
	elseif(li .eq. 6)then	!condition sur M a l'exterieur
c	 write(6,*)'limite',li
c	 write(6,2000)(y(i),i=1,ne)
	 be(1)=y(5)-mext**(2.d0/3.d0)			!en q=n m=mext(r,l)
	 ae(ne*(ne*(1-1)+3-1)+1)=-dmdr/3.d0/ray/m13	!derivee/r**2
	 ae(ne*(ne*(1-1)+4-1)+1)=-dmdl*l13/m13		!derivee/l**2/3
	 ae(ne*(ne*(1-1)+5-1)+1)=1.d0			!derivee/m**2/3
c	 write(6,2000)y(5)
	 if(deriv)then 	!test de derivation
	  write(6,*)'limite y(5),ln text',li
	  write(6,2000)(y(i),i=1,ne)
	  write(6,*)'y(5),mext,be(1)'
	  write(6,2000)y(5),mext,be(1)
	  do i=1,ne
	   stor0=y(i)
	   stor=stor0*unpdd
	   if(stor .eq. 0.)stor=dd
	   dstor=stor-stor0
	   y(i)=stor
	   ray=sqrt(y(3))
	   ln=sqrt(y(4))**3
	   call lim_ext(.false.,ln,ray,xchim,pext,text,dpdl,dpdr,dtdl,dtdr,
     1	teff,rtot,mext,dmdl,dmdr,dy1,dy1,dy1,dy1,dy1,mstar,
     2	tdetau,etat,opa)
	   bes(1)=y(5)-mext**(2.d0/3.d0)
	   dern(1)=(bes(1)-be(1))/dstor
	   y(i)=stor0
	   write(6,*)'derivee pour li=3',i
	   write(6,2000)dern(1),ae(ne*(ne*(1-1)+i-1)+1)
	  enddo	!i
	 endif		!deriv
c	 pause'limite li= '
	endif		!li
 
c	write(6,*)'li,be(1),y(1),y(2),y(5)',li
c	write(6,2000)be(1),y(1),y(2),y(5)
 
	return
 
	end
