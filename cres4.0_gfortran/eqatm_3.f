 
c************************************************************************
 
	subroutine eqatm_3(fait,knot,li,q,qt,xx,xl,mstar,
     1	derxx,lderxx,bp,xchim,cx,y,be,ae,ray,lum,tdetau,etat,opa,lu)
 
c	equations de l'atmosphere pour la methode de collocation
c	le rayon, la masse sont mesures en tau23 fixe au point n23
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	CESAM version 3
 
c entree:
c	fait=1 : initialisation
c	     2 : point courant	
c	     3 : limite
c	mstar: masse totale, temps t+dt, avec perte de masse
 
c entree/sortie
c	knot : nombre de points de table
c	li : indice de la limite
c	xx : points de collocation
c	xl : limites
c	q, qt: abscisses et points de table
c	derxx,lderxx : derivees aux points de collocation
c	bp : solution
c	xchim : composition chimique / gramme
c	cx : indice du point de collocation
c	y : variables
c	be,ae : equations et derivees
c	ray :rayon au fond de l'atmosphere
c	lum : luminosite=cte
c	lu=.true. : derivees inutiles
 
c routines externes
c	tdetau,etat,opa,lu,
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
	include 'evol_chim_3.common'
	
	integer pnpt
	parameter(pnpt=30)
 
	integer i,fait,knot,li,cx,j,ll,knot0,knot1,npt
 
	real*8 xx(1),xl(1),derxx(1),lderxx(1),bp(1),y(1),xchim(1),
     1	ro,drop,drot,drox,drott,drotp,drotx,ro0,tau,mstar,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,ae(1),be(1),
     3	cte1,cte2,cte3,cte4,grav,teff,dtsdtau,dtsdteff,dtsdg,
     4	dtbr23,to(pnpt),pd(pnpt),trn,kap0,tb0,w,
     6	kap,dkapt,dkapro,dkapx,tb,ray,lum,bid,dkapp,
     7	dgravre,dteffre,tau23,delfim,delfip,
     8	dy1,dy3,dy5,prn,bspint,stor,stor0,lamb,ltauf,ltau23,ltaue,
     9	dd,unpdd,dstor,bes(pne_atm),dern(pne_atm)
c	data dd,unpdd/0.00001,1.00001/
 
c	real*8 fx(7),dfxdx(7)
 
	real*8	vt1(pn_atm*pm_qs+2),vt2((pn_atm-1)*pm_qs+1),
     1	vt3((pn_atm*(pm_qs+1))**2),delfi,ddelfi,ltau(pn_atm),lnp(pn_atm),
     2	lnt(pn_atm),q(1),q1(pn_atm*pm_qs+2),qt(1),
     3	psi0(0:2),psix(0:2),psi(0:2),ro_ext,dro_grav,dro_teff
 
	logical deriv,lu,para
 
	external tdetau,etat,opa
 
	data dd,unpdd/0.00001,1.00001/
 
c	pressions en valeurs naturelles
 
	data pd
     &       /6.300D+02,6.360D+02,6.474D+02,6.687D+02,7.088D+02,
     &        7.827D+02,9.164D+02,1.152D+03,1.556D+03,2.229D+03,
     &        3.328D+03,5.092D+03,7.893D+03,1.230D+04,1.914D+04,
     &        2.960D+04,4.525D+04,6.781D+04,9.843D+04,1.367D+05,
     &        1.790D+05,2.177D+05,2.495D+05,2.719D+05,2.856D+05,
     &        2.933D+05,2.972D+05,2.992D+05,3.002D+05,3.007D+05/
 
c	profondeurs optiques en ln
 
	data to
     &      /-1.451D+01,-1.317D+01,-1.253D+01,-1.188D+01,-1.124D+01,
     &       -1.059D+01,-9.946D+00,-9.301D+00,-8.656D+00,-8.011D+00,
     &       -7.366D+00,-6.721D+00,-6.076D+00,-5.431D+00,-4.786D+00,
     &       -4.141D+00,-3.496D+00,-2.851D+00,-2.206D+00,-1.561D+00,
     &       -9.163D-01,-3.028D-01, 3.107D-01, 9.242D-01, 1.538D+00,
     &        2.151D+00, 2.765D+00, 3.378D+00, 3.992D+00, 5.011D+00/
 
c	external tdetau,etat,opa
 
	save cte1,cte2,cte3,grav,teff,pd,to,psix,psi0,dd,unpdd
 
2000	format(1x,1p8d10.3)
 
c	write(6,*)'entree eqatm_3, tau_min,tau_max,n_atm,n23,nea',
c	1	tau_min,tau_max,n_atm,n23,nea
 
	goto(100,200,300),fait
 
100	cte1=g*msol/rsol**2
	cte2=lsol/pi/rsol**2/aradia/clight	!sigma=a c /4
	cte3=4*pi*rsol**2/msol
	cte4=2./3.*rsol
 
	nea=7		!nombre d'inconnues
 
c	a l'exterieur
 
	teff=(cte2*lum/ray**2)**.25
	grav=cte1*mstar/ray**2
	
	call taueff_3(teff,grav,tdetau,tau23)	!initialisation de tau23
c	write(6,2000)teff,grav,tau23,lum,ray
c	pause'initialisation'
 
	npt=pnpt		!nb. de points tabules
	do i=1,npt		!pressions en ln
         pd(i)=log(pd(i))
	enddo
 
c	initialisation des tau=exp(fi)
 
c	write(6,*)'eqatm les tau'
c	write(6,2000)tau_max,tau_min,tau23
 
c	on fait le changement de variable tau --> psi lineaire par morceaux
c	si (para=.false.) ou quadratique si (para=.true.)
c	avec le cht. quadratique on risque que psi ne soit pas monotone,
c	avec le lineaire par morceaux la CV de NR est moins bonne
 
	para=.false.
	if(para)then
	 psi0(0)=log(tau_max)
	 psi0(1)=log(tau23)
	 psi0(2)=log(tau_min)
	 psix(0)=1
	 psix(1)=n23
	 psix(2)=n_atm
	 do i=1,n_atm
	  call newton(0,2,psi0,psix,dfloat(i),psi,0)
	  ltau(i)=psi(0)
	 enddo
	else
	 ltauf=log(tau_max)
	 ltaue=log(tau_min)
	 ltau23=log(tau23)
	 delfim=1.d0/float(n23-1)
	 delfip=1.d0/float(n23-n_atm)
 
	 do i=1,n_atm
	  if(i .le. n23)then
           delfi=(ltau23-ltauf)*delfim
	   ltau(i)=ltauf+delfi*(i-1)
	  else
           delfi=(ltau23-ltaue)*delfip
	   ltau(i)=ltaue+delfi*(i-n_atm)
	  endif
	 enddo
	endif
 
c	write(6,*)'n_atm,n23,m_qs',n_atm,n23,m_qs
c	write(6,*)'tau_max,tau23,tau_min'
c	write(6,2000)tau_max,tau23,tau_min
c	write(6,*)'ln tau'
c	write(6,2000)(ltau(i),i=1,n_atm)
c	write(6,*)'tau'
c	write(6,2000)(exp(ltau(i)),i=1,n_atm)
c	pause 'les tau'
 
c	initialisation des coefficients d'interpolation des pressions
 
        bid=bspint(pd,to,vt1,npt,2,knot0,.false.,ltau(1),ll)
 
c	calcul des pressions aux points ltau
 
	do i=1,n_atm
	 lnp(i)=bspint(pd,to,vt1,npt,2,knot0,.false.,ltau(i),ll)
	enddo
c	write(6,*)'pressions initiales'
c	write(6,2000)(exp(lnp(i)),i=1,n_atm)
 
c	projection des pressions et temperatures sur une base de
c	B-splines en ltau
c	ltau etant decroissant on ecrase -ltau sur derxx encore inutile
 
	do i=1,n_atm
	 derxx(i)=-ltau(i)
	enddo
 
c	projection des pressions
 
	bid=bspint(lnp,derxx,vt1,n_atm,2,knot0,.false.,-ltau(1),ll)
c	write(6,*)'second',bid
 
c	temperatures aux points ltau
 
	do i=1,n_atm
	 call tdetau(exp(ltau(i)),teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
	 lnt(i)=log(tb)
c	 write(6,2000)ltau(i),exp(ltau(i)),teff,grav,tb
	enddo
c	write(6,*)'temperatures initiales',teff,grav
c	write(6,2000)(exp(ltau(i)),i=1,n_atm)
c	write(6,2000)(exp(lnt(i)),i=1,n_atm)
 
c	projection des temperatures
 
	bid=bspint(lnt,derxx,vt1,n_atm,2,knot0,.false.,-ltau(1),ll)
c	print*,knot0,ll
c	write(6,2000)-ltau(1),bid
c	write(6,2000)(lnt(i),i=1,n_atm)
c	write(6,2000)(derxx(i),i=1,n_atm)
c	write(6,2000)(vt1(i),i=1,knot0)
c	pause'projection des temperatures'
		
c	abscisses pour l'integration
 
	do i=1,n_atm
	 q(i)=i					!les abscisses
	enddo
 
c	projection des ltau sur une base de B-splines en q
 
	bid=bspint(ltau,q,q1,n_atm,2,knot1,.false.,q(1),ll)
c	write(6,*)'4',bid
 
c	vecteur nodal en ltau pour l'integration
 
        call  noedif(derxx,xx,n_atm,m_qs,1,j)		!derxx,xx,j:VT
 
c	on passe sur la base de B-splines en q
 
	call newsp1(derxx,vt1,knot0,2,xx,j,m_qs+1,lnp,vt2,y,vt3)
	do i=1,(n_atm-1)*m_qs+1
	 bp(nea*(i-1)+1)=vt2(i)		!pour la pression
	enddo
c	write(6,*)'splines en lnp'
c	write(6,2000)(bp(nea*(i-1)+1),i=1,(n_atm-1)*m_qs+1,mm_qs)
c	pause'apres newsp1'
 
	call newsp1(derxx,vt1,knot0,2,xx,j,m_qs+1,lnt,vt2,y,vt3)
	do i=1,(n_atm-1)*m_qs+1
	 bp(nea*(i-1)+2)=vt2(i)		!pour la temperature
	enddo
c	write(6,*)'splines en lnt'
c	write(6,2000)(bp(nea*(i-1)+2),i=1,(n_atm-1)*m_qs+1,m_qs)
c	pause'apres newsp1 2'
 
c	vecteur nodal pour l'integration
 
        call  noedif(q,qt,n_atm,m_qs,1,knot)
c	write(6,2000)(qt(i),i=1,knot)
c	pause'vecteur nodal'
 
c	les ltau dans la base des qt
 
	call newsp1(q,q1,knot1,2,qt,knot,m_qs+1,ltau,vt2,y,vt3)
	do i=1,(n_atm-1)*m_qs+1
	 bp(nea*(i-1)+7)=vt2(i)		!pour les profondeurs optiques
	enddo
c	write(6,*)'splines en ltau'
c	write(6,2000)(bp(nea*(i-1)+7),i=1,(n_atm-1)*m_qs+1,m_qs)
c	pause'apres newsp1 3'
 
c	pour les autres variables
 
	do i=1,(n_atm-1)*m_qs+1
	 bp(nea*(i-1)+3)=ray		!r/Rsol
	 bp(nea*(i-1)+4)=ray		!rn23/Rsol
	 bp(nea*(i-1)+5)=mstar		!m/Msol
	 bp(nea*(i-1)+6)=log(tau23)	!ln Tau23
	enddo
 
c	limites
 
	xl(1)=1		!en tau=tau_max, q=1 R/Rsol=y(3)=Rb/Rsol
	xl(2)=n23	!en tau=tau*, q=n23 ln T =y(2)=ln T(tau,Teff)
	xl(3)=n23	!en tau=tau*, q=n23 ln T =y(2)=ln Teff
	xl(4)=n23	!en tau=tau*, q=n23 ln tau =y(7)=ln tau*=y(6)
	xl(5)=n23	!en tau=tau*, q=n23 R/Rsol=y(3)=y(4)=R23/Rsol
	xl(6)=n23	!en tau=tau*, q=n23 M=y(5)=M*
	xl(7)=n_atm	!en tau=tau_min, q=N ro=ro_ext
c	write(6,*)'limites',nea,(xl(i),i=1,nea)
 
	call colloc(m_qs,1,n_atm,q,qt,knot,xx,vt1,derxx,nea,lderxx,xl)
 
c	do i=1,n_atm
c	 call sbsp1dn(nea,bp,q,qt,n_atm,m_qs+1,knot,.true.,q(i),ll,fx,dfxdx)
c	 write(6,2000)(exp(fx(j)),j=1,2),(fx(j),j=3,5),(exp(fx(j)),j=6,7)
c	enddo
 
	return
 
c------------------	equations	--------------------------
 
c              y(1)=ln p               y(8)=(ln P)'
c              y(2)=ln T               y(9)=(ln T)'
c              y(3)=R/Rsol             y(10)=(R/Rsol)'
c              y(4)=R23/Rsol           y(11)=(R23/Rsol)'
c              y(5)=M/Msol             y(12)=(M/Msol)'
c              y(6)=ln tau23           y(13)=(ln tau23)'
c              y(7)=ln tau             y(14)=(ln tau)'
 
c	xx(cx) valeur de tau au point de collocation
 
c	ae(ne(ne(der-1)+var-1)+eq)=derivee de la eq-ieme equation par
c	rapport a la (der-1)-ieme derivee de la var-ieme variable
 
200	prn=exp(y(1))		!pression cgs
	trn=exp(y(2))		!temperature cgs
	tau=exp(y(7))		!profondeur optique
	
	call etat(prn,trn,xchim,.not.lu,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	call opa(xchim,trn,ro,kap,dkapt,dkapro,dkapx)
 
c	print*,'prn,trn,tau,xchim(1),ro,kap'
c	write(6,2000)prn,trn,tau,xchim(1),ro,kap
	
	if(.not.lu)then
	 drop=drop*prn		!d ro /d ln p
	 dkapp=drop*dkapro/kap	!d kap/d ln p a kap pres
	 drop=drop/ro		!d ro /d ln p a ro pres
	 drot=drot*trn				!d ro /d ln t
	 dkapt=(drot*dkapro+dkapt*trn)/kap	!d kap/d ln t a kap pres
	 drot=drot/ro				!d ro /d ln t a ro pres
	endif
 
c	Definition de la fonction fi
 
	if(para)then
	 psi0(1)=y(6)
	 call newton(0,2,psi0,psix,xx(cx),psi,1)
	 delfi=psi(1)
	 ddelfi=(2*xx(cx)-n_atm-1)/(n23-1)/(n23-n_atm)/delfi	!a delfi pres
	else
         if (xx(cx) .le. n23) then
          delfi=(y(6)-ltauf)*delfim
          ddelfi=delfim   !derivee fi /y(6)
         else
          delfi=(y(6)-ltaue)*delfip
          ddelfi=delfip   !derivee fi /y(6)
         endif
	 ddelfi=ddelfi/delfi		!a delfi pres
	endif
 
c	equations
 
	if(iw .gt. 0)then		!la rotation
	 w=xchim(iw)/y(3)**2
	else
	 w=w_rot
	endif
	w=w**2
 
	dy1=(cte1*y(5)/y(3)**2-cte4*w*y(3))*delfi/kap*exp(y(7)-y(1))
	be(1)=y(8)-dy1				!d log p/d tau=gm/kr**2/p
	if(.not.lu)then
	 ae(nea*(nea*(1-1)+1-1)+1)=dy1*(1.d0+dkapp)		!derivee /log p
	 ae(nea*(nea*(2-1)+1-1)+1)=1.d0				!derivee /log p'
	 ae(nea*(nea*(1-1)+2-1)+1)=dy1*dkapt			!derivee /log t
	 ae(nea*(nea*(1-1)+3-1)+1)=(2.*cte1*y(5)/y(3)**3+cte4*w)*
     1	delfi/kap*exp(y(7)-y(1))			!derivee /r
	 ae(nea*(nea*(1-1)+5-1)+1)=-cte1/y(3)**2*
     1			delfi/kap*exp(y(7)-y(1))	!derivee /m
         ae(nea*(nea*(1-1)+6-1)+1)=-ddelfi*dy1               !derive/logtau23
         ae(nea*(nea*(1-1)+7-1)+1)=-dy1                      ! derivee /logtau
	endif
 
	teff=(cte2*lum/y(4)**2)**.25
	grav=cte1*mstar/y(4)**2
	
c	print*,'tau,teff,grav'
c	write(6,2000)tau,teff,grav
	call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
c	write(6,2000)tb,ro_ext
	be(2)=y(2)-log(tb)
        if(.not.lu)then
	 dteffre=-teff/y(4)/2.				!derivee teff / d r23
	 dgravre=-2.*grav/y(4)				!derivee grav / d r23
	 dtbr23=dtsdteff*dteffre+dtsdg*dgravre		!d tb /d r23
         ae(nea*(nea*(1-1)+2-1)+2)=1.d0 			!derivee /log t
         ae(nea*(nea*(1-1)+4-1)+2)=-dtbr23/tb 		!derivee /r23/Rsol
         ae(nea*(nea*(1-1)+7-1)+2)=-dtsdtau*tau/tb	 !derivee /log tau
        endif
 
        dy3=-1.d0/(kap*ro*Rsol)*tau*delfi
        be(3)=y(10)-dy3
	if(.not.lu)then
	 ae(nea*(nea*(1-1)+1-1)+3)=dy3*(drop+dkapp)		!derivee /log p
         ae(nea*(nea*(1-1)+2-1)+3)=dy3*(drot+dkapt)             !derivee /log t
	 ae(nea*(nea*(2-1)+3-1)+3)=1.d0			!derivee /r/rsol'
         ae(nea*(nea*(1-1)+6-1)+3)=-dy3*ddelfi          !derivee /logtau23
         ae(nea*(nea*(1-1)+7-1)+3)=-dy3                 !derivee /logtau
	endif
 
        be(4)=y(11)
        ae(nea*(nea*(2-1)+4-1)+4)=1.d0                    !derivee /r23/Rsol'
 
        dy5=-cte3*y(3)**2*tau*delfi/kap
        be(5)=y(12)-dy5
	if(.not.lu)then
         ae(nea*(nea*(1-1)+1-1)+5)=dy5*dkapp                    !derivee /log(p)
         ae(nea*(nea*(1-1)+2-1)+5)=dy5*dkapt                    !derivee /log(t)
         ae(nea*(nea*(1-1)+3-1)+5)=-2.d0*dy5/y(3)               !derivee /R/Rsol
         ae(nea*(nea*(2-1)+5-1)+5)=1.d0                     !derivee /M/Mstar'
         ae(nea*(nea*(1-1)+6-1)+5)=-dy5*ddelfi              !derivee /log(tau23)
         ae(nea*(nea*(1-1)+7-1)+5)=-dy5		            !derivee /log(tau)
	endif
 
        be(6)=y(13)
	ae(nea*(nea*(2-1)+6-1)+6)=1.d0                      !derivee /log(tau23)'
 
        be(7)=y(14)-delfi
        ae(nea*(nea*(1-1)+6-1)+7)=-ddelfi*delfi          !derivee /log(tau23)
        ae(nea*(nea*(2-1)+7-1)+7)=1.d0                     !derivee /log(tau)'
 
c	deriv=cx .le. 3
c	1	.or. (cx .ge. n23-1 .and. cx .le. n23+1)
c	2	.or. (cx .ge. n_atm-1 .and. cx .le. n_atm)
 
c	deriv=.true.
c	deriv=.false.
	if(.not. deriv)return	!tests de mise au point pour derivees
 
	write(6,*)'y(i)/dy'
	write(6,2000)(y(i),i=1,nea)
	write(6,2000)(y(i),i=nea+1,2*nea)
	write(6,*)'prn,teff,grav,dy1,dy3,dy5,kap,tau'
	write(6,2000)prn,teff,grav,dy1,dy3,dy5,kap,tau
	write(6,*)'iw,iz/delfi,w',iw,iz
	write(6,2000)delfi,w,cte1*y(5)/y(3)**2,cte4*w*y(3)
c	write(6,2000)(be(i),i=1,nea)
c       write(6,*)'cte2',cte1,cte2,lum,lsol,pi,rsol,aradia,clight,nea
	write(6,*)'dkapp,dkapt,dtbr23,drop,drot,dtsdtau',xx(cx),cx
	write(6,2000)dkapp,dkapt,dtbr23,drop,drot,dtsdtau*tau/tb
	write(6,*)' '
	kap0=kap
	tb0=tb
	ro0=ro
        do i=1,nea
         stor0=y(i)
	 stor=stor0*unpdd	
         if(stor .eq. 0.)stor=dd
	 dstor=stor-stor0
	 y(i)=stor
 
	 if(para)then
	  psi0(1)=y(6)
	  call newton(0,2,psi0,psix,xx(cx),psi,1)
	  delfi=psi(1)
	 else
          if (xx(cx) .le. n23) then
           delfi=(y(6)-ltauf)*delfim
          else
           delfi=(y(6)-ltaue)*delfip
          endif
	 endif
	 prn=exp(y(1))		!pression cgs
	 trn=exp(y(2))		!temperature cgs
	 tau=exp(y(7))		!profondeur optique
	 call etat(prn,trn,xchim,.false.,
     & 	 ro,drop,drot,drox,drott,drotp,drotx,
     &    u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	 call opa(xchim,trn,ro,kap,dkapt,dkapro,dkapx)
 
	 dy1=(cte1*y(5)/y(3)**2-cte4*w*y(3))*delfi/kap*exp(y(7))/exp(y(1))
	 bes(1)=y(8)-dy1
 
	 teff=(cte2*lum/y(4)**2)**.25
	 grav=cte1*mstar/y(4)**2
	 call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
	 bes(2)=y(2)-log(tb)
 
	 dy3=-1.d0/(kap*ro*Rsol)*tau*delfi
	 bes(3)=y(10)-dy3
 
	 bes(4)=y(11)
 
         dy5=-cte3*y(3)**2*tau*delfi/kap
         bes(5)=y(12)-dy5
 
         bes(6)=y(13)
 
         bes(7)=y(14)-delfi
 
	 do j=1,nea
	  dern(j)=(bes(j)-be(j))/dstor
	 enddo	!j
 
	 y(i)=stor0
	 write(6,*)'derivee',i
	 write(6,2000)(dern(j),ae(nea*(nea*(1-1)+i-1)+j),j=1,nea)
     1	,(kap-kap0)/dstor/kap,(tb-tb0)/dstor,(ro-ro0)/dstor/ro
 
	enddo	!i
 
	call tdetau(tau,teff,grav,tb0,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
        stor0=tau
	stor=stor0*unpdd	
        if(stor .eq. 0.)stor=dd
	dstor=stor-stor0
	tau=stor
	call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
	write(6,2000)(tb-tb0)/dstor,dtsdtau
	tau=stor0
 
        stor0=teff
	stor=stor0*unpdd	
        if(stor .eq. 0.)stor=dd
	dstor=stor-stor0
	teff=stor
	call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
	write(6,2000)(tb-tb0)/dstor,dtsdteff
c	write(6,*)stor,stor0,tb,tb0
	teff=stor0
 
        stor0=grav
	stor=stor0*unpdd	
        if(stor .eq. 0.)stor=dd
	dstor=stor-stor0
	grav=stor
	call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
	write(6,2000)(tb-tb0)/dstor,dtsdg
	grav=stor0
 
	pause
 
	return
 
c-------------	limites	---------------------------
 
300	if(li .eq. 1)then	!tau=tau_max, q=1 R/Rsol=y(3)=Rb/Rsol=ray
         be(1)=y(3)-ray
         ae(nea*(nea*(1-1)+3-1)+1)=1.d0             !derivee /R/Rsol
 
	elseif(li .eq. 2)then 	!tau=tau_max, q=1, ln T=y(2)=ln T(tau,Teff,..)
	 tau=exp(y(7))		!profondeur optique
	 teff=(cte2*lum/y(4)**2)**.25
	 grav=cte1*mstar/y(4)**2
	 call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
  	 be(1)=y(2)-log(tb)
	 dteffre=-teff/y(4)/2.				!derivee teff / d r23
	 dgravre=-2.*grav/y(4)				!derivee grav / d r23
	 dtbr23=dtsdteff*dteffre+dtsdg*dgravre		!d tb /d r23
	 ae(nea*(nea*(1-1)+2-1)+1)=1.d0 			!derivee /log t
	 ae(nea*(nea*(1-1)+4-1)+1)=-dtbr23/tb 		!derivee /r23/Rsol
	 ae(nea*(nea*(1-1)+7-1)+1)=-dtsdtau*tau/tb !derivee /log tau
 
	 deriv=.false.
c	 deriv=.true.
	 if(deriv)then
	  deriv=.false.
	  write(6,*)'limite',li
          do i=1,nea
           stor0=y(i)
	   stor=stor0*unpdd	
           if(stor .eq. 0.)stor=dd
	   dstor=stor-stor0
	   y(i)=stor
	   tau=exp(y(7))		!profondeur optique
	   teff=(cte2*lum/y(4)**2)**.25
	   grav=cte1*mstar/y(4)**2
	   call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
	   bes(1)=y(2)-log(tb)
	   y(i)=stor0
	   write(6,*)'derivee',i
	   write(6,2000)(bes(1)-be(1))/dstor,ae(nea*(nea*(1-1)+i-1)+1)
	  enddo
	  pause
	 endif
 
	elseif(li .eq. 3)then 	!tau=tau*, q=n23, T=Teff
	 teff=(cte2*lum/y(4)**2)**.25
	 be(1)=y(2)-log(teff)
	 dteffre=-teff/y(4)/2.				!derivee teff / d r23
	 ae(nea*(nea*(1-1)+2-1)+1)=1.d0 			!derivee /log t
	 ae(nea*(nea*(1-1)+4-1)+1)=-dteffre/teff 	!derivee /r23/Rsol
 
	 deriv=.false.
c	 deriv=.true.
	 if(deriv)then
	  deriv=.false.
	  write(6,*)'limite',li,xl(li)
          do i=1,nea
           stor0=y(i)
	   stor=stor0*unpdd	
           if(stor .eq. 0.)stor=dd
	   dstor=stor-stor0
	   y(i)=stor
	   teff=(cte2*lum/y(4)**2)**.25
	   bes(1)=y(2)-log(teff)
	   y(i)=stor0
	   write(6,*)'derivee',i
	   write(6,2000)(bes(1)-be(1))/dstor,ae(nea*(nea*(1-1)+i-1)+1)
	  enddo
	  pause
	 endif
 
	elseif(li .eq. 4)then 	!tau=tau*
	 be(1)=y(6)-y(7)
	 ae(nea*(nea*(1-1)+6-1)+1)=1.d0			!derivee /log tau23
	 ae(nea*(nea*(1-1)+7-1)+1)=-1.d0		!derivee /log tau
 
	elseif(li .eq. 5)then	!tau=tau*, q=n23 R/Rsol=y(3)=y(4)=R23/Rsol
	 be(1)=y(3)-y(4)
         ae(nea*(nea*(1-1)+3-1)+1)=1.d0           !derivee /R/Rsol
         ae(nea*(nea*(1-1)+4-1)+1)=-1.d0          !derivee /R23/Rsol
 
	elseif(li .eq. 6)then	!tau=tau*, q=n23 M=y(5)=M*/Msol
	 be(1)=y(5)-mstar
         ae(nea*(nea*(1-1)+5-1)+1)=1.d0           !derivee /M/Mstar
 
	elseif(li .eq. 7)then	!en tau=tau_min, q=N ro=ro_ext
	 teff=(cte2*lum/y(4)**2)**.25
	 dteffre=-teff/y(4)/2.				!derivee teff / d r23
	 grav=cte1*mstar/y(4)**2
	 dgravre=-2.*grav/y(4)				!derivee grav / d r23
	 call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
         prn=exp(y(1))
         trn=exp(y(2))
         call etat(prn,trn,xchim,.not.lu,
     &                ro,drop,drot,drox,drott,drotp,drotx,
     &                u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	 be(1)=ro-ro_ext                         !ro(n)=ro_ext
	 ae(nea*(nea*(1-1)+1-1)+1)=drop*prn     !derivee /log p
         ae(nea*(nea*(1-1)+2-1)+1)=drot*trn     !derivee /log t
         ae(nea*(nea*(1-1)+4-1)+1)=-dro_teff*dteffre
     1	-dro_grav*dgravre    		 !derivee /r23
 
	 deriv=.false.
c	 deriv=.true.
	 if(deriv)then
	  deriv=.false.
	  write(6,*)'limite',li,xl(li)
          do i=1,nea
           stor0=y(i)
	   stor=stor0*unpdd	
           if(stor .eq. 0.)stor=dd
	   dstor=stor-stor0
	   y(i)=stor
	   teff=(cte2*lum/y(4)**2)**.25
	   grav=cte1*mstar/y(4)**2
	   call tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
           prn=exp(y(1))
           trn=exp(y(2))
           call etat(prn,trn,xchim,.not.lu,
     &                ro,drop,drot,drox,drott,drotp,drotx,
     &                u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	   bes(1)=ro-ro_ext                         !ro(n)=ro_ext
	   y(i)=stor0
	   write(6,*)'derivee',i
	   write(6,2000)(bes(1)-be(1))/dstor,ae(nea*(nea*(1-1)+i-1)+1)
	  enddo
	  pause
	 endif
	else
	 write(6,*)'erreur dans li, li=',li
	 stop
	endif
 
	return
 
	end
 
