c****************************************************************

	subroutine conv_cm_new_3(krad,grav,delta,cp,ro,hp,taur,gradrad,gradad,
     &  der,grad,dgradk,dgradgr,dgradel,dgradcp,dgradro,
     &  dgradtaur,dgradhp,dgradrad,dgradad)
	
c	calcul du gradient convectif selon Canuto Mazitelli ApJ 370, 295, 1991

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM4
c
c       Updated version (30.04.02) by 
c	Reza Samadi, LESIA , Observatoire de Paris 
c       New Version in which Eq.(6) of CM is modified so as to take into
c       account the quantity delta=- ( d ln ro / d ln T )p

c ------------------------------------
c Adapté à CRES, D.C. 7 novembre 2002.
c ------------------------------------

c       
c entree :
c	krad : conductivite radiative 4 ac T**3 / 3 kap ro
c	gravite : G m / r**2
c	delta : - ( d ln ro / d ln T )p
c	cp : chaleur specifique
c	ro : densite
c	hp : echelle de hauteur de pression
c	gradrad : gradient radiatif
c	gradad : gradient adiabatique
c	taur : ep. opt. de la bulle
c	der=.true. : calcul des derivees

c	alpha : longueur de melange (dans le common /modele/)

c sortie :
c	gradconv : gradient convectif
c	dgrad* : derivees
c	gam : le gamma de la convection
c	dgam* : derivees

	implicit none

	include 'modele_3.common'
	include 'ctephy.common'
	!include 'convection.common'

	integer iter

	real*8	krad,grav,delta,cp,ro,hp,taur,gradrad,gradad,
	1	grad,dgradk,dgradgr,dgradel,dgradcp,dgradro,
	2	dgradhp,dgradtaur,dgradrad,dgradad,
	3	gam,dgamk,dgamgr,dgamdel,dgamcp,dgamro,
	4	dgamhp,dgamtaur,dgamrad,dgamad
	
	real*8	a1,a2,m,n,p,l,ki,a,corr,b,a2s,da2s,
	1	sig,dsig,a2sn,da2sn,a2s1,phi,phip,dphi,f,df
	
	real*8	dkik,dkicp,dkiro,dbhp,dbk,dbcp,dbro,dbgr,dahp,dak,dacp,daro,
	1	dagr,dsighp,dsigk,dsigcp,dsigro,dsiggr,dsigad,da2shp,da2sk,
	2	da2scp,da2sro,da2sgr,da2sad,da2snhp,da2snk,da2sncp,
	3	da2snro,da2sngr,da2snad,dphihp,dphik,dphicp,dphiro,dphigr,
	4	dphiad,dgrad,sg12,dgams
		
	logical init,der,conv
	data init/.false./
	
	save
	
2000	format(1x,1p8d10.3)

	if(.not. init)then	!vsal: V/Al
	 init=.true.
	 	 
	 write(2,1)
	 write(6,1)
1	 format(/,1x,'gradient dans zone convective calcule selon',/,
	1	1x,' Canuto-Mazitelli ApJ 370, 295, 1991')
	
	 a1=24.868d0
	 a2=9.766d-2
	 m=0.14972d0
	 n=0.18931d0
	 p=1.8503d0
	endif
	
	l=alpha*hp		!longueur de melange
	ki=krad/cp/ro		!conductivite thermometrique
	b=2.d0*l**2/9.d0/ki*sqrt(grav/2.d0/hp*delta)		!2 x (6)
	! eq.(6) corrected so as to include delta:
	! delta= - (dln rho / ln T)_P                     [RS 30.04.02]
	a=b**2			!c'est le 4a**2 de CM 
	
c	initialisations
		
	conv=.false.
	grad=gradad*1.1
	iter=0
	phi=1.d30
	
c	Newton-Raphson "d" signifie derivee/ grad	
	
	do while(.not.conv)
	 iter=iter+1
	 if(iter .gt. 30)then
	  write(6,*)'pas de convergence dans conv_cm'
	  print*,'donnees : krad,grav,cp,ro,hp,gradrad,gradad'	  
	  write(6,2000)krad,grav,cp,ro,hp,gradrad,gradad
	  print*,'non convergence : phi,phip,abs(phi-phip)/phi'
	  write(6,2000)phi,phip,abs(phi-phip)/phi
	  stop
	 endif
	 phip=phi
	 sig=a*(grad-gradad)	!(5)
	 dsig=a
	 a2s=a2*sig		!(32)
	 da2s=a2*dsig
	 a2s1=1.d0+a2s
	 a2sn=a2s1**n
	 da2sn=n*a2sn/a2s1*da2s
	 phi=a1*sig**m*(a2sn-1.d0)**p
	 dphi=phi*(m*dsig/sig+p*da2sn/(a2sn-1.d0))
	 f=grad+phi*(grad-gradad)-gradrad		!(63)
	 df=1.d0+dphi*(grad-gradad)+phi
	 corr=f/df
c	 print*,iter
c	 write(6,2000)corr,abs(phi-phip)/phi
	 grad=grad-corr
	 conv=abs(phi-phip)/phi .le. 1.d-10
	enddo
c	pause
	
	sg12=sqrt(sig+1.d0)
	gam=(sg12-1.d0)/2.d0
	
	if(der)then	
	 dkik=ki/krad
	 dkicp=-ki/cp
	 dkiro=-ki/ro
	
	 dbhp=b*(2.d0*alpha/l-.5d0/hp)
	 dbk= -b/ki*dkik
	 dbcp=-b/ki*dkicp
	 dbro=-b/ki*dkiro
	 dbgr=b*.5/grav
	
	 dahp=2.d0*b*dbhp
	 dak= 2.d0*b*dbk
	 dacp=2.d0*b*dbcp
	 daro=2.d0*b*dbro
	 dagr=2.d0*b*dbgr
	 
	 sig=a*(grad-gradad)	!(5)
	 dsig=a
	 dsighp=sig/a*dahp
	 dsigk= sig/a*dak
	 dsigcp=sig/a*dacp
	 dsigro=sig/a*daro	
	 dsiggr=sig/a*dagr
	 dsigad=-a
	
	 a2s=sig*a2		!(32)
	 a2s1=1.d0+a2s
	 da2s=    dsig*a2	 
	 da2shp=dsighp*a2
	 da2sk= dsigk *a2	
	 da2scp=dsigcp*a2
	 da2sro=dsigro*a2
	 da2sgr=dsiggr*a2
	 da2sad=dsigad*a2
	 
	 a2sn=a2s1**n
	 da2sn=  n*a2sn/a2s1*da2s
	 da2snhp=n*a2sn/a2s1*da2shp	
	 da2snk= n*a2sn/a2s1*da2sk	
	 da2sncp=n*a2sn/a2s1*da2scp	
	 da2snro=n*a2sn/a2s1*da2sro	
	 da2sngr=n*a2sn/a2s1*da2sgr	
	 da2snad=n*a2sn/a2s1*da2sad
	 
	 phi=a1*sig**m*(a2sn-1.d0)**p
	 dphi=  phi*(m*dsig  /sig+p*da2sn  /(a2sn-1.d0))
	 dphihp=phi*(m*dsighp/sig+p*da2snhp/(a2sn-1.d0))
	 dphik= phi*(m*dsigk /sig+p*da2snk /(a2sn-1.d0))		
	 dphicp=phi*(m*dsigcp/sig+p*da2sncp/(a2sn-1.d0))	
	 dphiro=phi*(m*dsigro/sig+p*da2snro/(a2sn-1.d0))	
	 dphigr=phi*(m*dsiggr/sig+p*da2sngr/(a2sn-1.d0))
	 dphiad=phi*(m*dsigad/sig+p*da2snad/(a2sn-1.d0))
	
	 dgrad=grad-gradad
	
	 dgradhp=-dphihp*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradk=-dphik*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradcp=-dphicp*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradro=-dphiro*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradgr=-dphigr*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradad=-(dphiad*dgrad-phi)/(1.d0+phi+dphi*dgrad)
	 dgradrad=1.d0/(1.d0+phi+dphi*dgrad)	 
	 dgradel=0.d0
	 dgradtaur=0.d0
	 
	 dgams=.25d0/sg12
	 dgamhp=dgams*(dsighp+dsig*dgradhp)
	 dgamk= dgams*(dsigk+dsig*dgradk)	 
	 dgamcp=dgams*(dsigcp+dsig*dgradcp)
	 dgamro=dgams*(dsigro+dsig*dgradro)	 
	 dgamgr=dgams*(dsiggr+dsig*dgradgr)	 
	 dgamad=dgams*(dsigad+dsig*dgradad)
	 dgamrad=dgams*dsig*dgradrad	 	 
	 dgamdel=0.d0
	 dgradtaur=0.d0	 	 	 
	endif
		
	return

	end

