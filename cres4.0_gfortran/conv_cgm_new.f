
c****************************************************************

	subroutine conv_cgm_new(krad,grav,delta,cp,ro,hp,taur,gradrad,gradad,der,
	1	grad,dgradk,dgradgr,dgradel,dgradcp,dgradro,
	2	dgradhp,dgradtaur,dgradrad,dgradad,
	3	gam,dgamk,dgamgr,dgamdel,dgamcp,dgamro,
	4	dgamhp,dgamtaur,dgamrad,dgamad)
	
c	calcul du gradient convectif selon CGM : Canuto Goldman Mazitelli ApJ 473, 550-559, 1996
c       et en adoptant la prescription de Bernkopf
c       on se refere ici au 2 articles suivant :
c         H02 : Heiter et al, A&A, 2002
c         CM  : Canuto Mazitelli ApJ 370, 295, 1991

c       Date:  1/03/2002  version pour CESAM4
c       Reza Samadi, LESIA , Observatoire de Paris 
c	adapte programme conv_cm de P. Morel  (Cassini, O.C.A., Observatoire de Nice)
c
c       Updated version (30.04.02) by Reza Samadi
c       New Version in which Eq.(6) of CM is modified so as to take into
c       account the quantity delta=- ( d ln ro / d ln T )p

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

	include 'modele.common'
	include 'ctephy.common'
	include 'convection.common'

	integer iter

	real*8	krad,grav,delta,cp,ro,hp,taur,gradrad,gradad,
	1	grad,dgradk,dgradgr,dgradel,dgradcp,dgradro,
	2	dgradhp,dgradtaur,dgradrad,dgradad,
	3	gam,dgamk,dgamgr,dgamdel,dgamcp,dgamro,
	4	dgamhp,dgamtaur,dgamrad,dgamad

	real*8	K0,ca,cb,ck,cm,cn,cc,cd,ce,cf,cpp,cq,cr,ct,
	1 l,ki,a,corr,b,phi,phip,dphi,f,df, sig,S,dS,S_sig,F1,F2,
	2 dF1,dF2,cSp,dSq1,eSr,fSt1,bS1,bS1m

	real*8	dkik,dkicp,dkiro,dbhp,dbk,dbcp,dbro,dbgr,dahp,dak,dacp,daro,
	1	dagr,dphihp,dphik,dphicp,dphiro,dphigr,dphiad,dgrad,
	2	dsig,dsighp,dsigk,dsigcp,dsigro,dsiggr,dsigad,sg12,dgams
		
	logical init,der,conv
	data init/.false./
	
	save
	
2000	format(1x,1p8d10.3)

	if(.not. init)then	!vsal: V/Al
	 init=.true.
	 	 
	 write(2,1)
	 write(6,1)
1	 format(/,1x,'gradient dans zone convective calcule selon',/,
	1	1x,' Canuto Goldman Mazitelli (CGM) ApJ 473, 550, 1996') 
c       CGM coefficients :
	 K0=1.7d0
	 S_sig=81d0/2d0
cc	Eq.(19) H02 
         ca=10.8654d0
	 cb=0.00489073d0
         ck=0.149888d0       
	 cm=0.189238d0
	 cn=1.85011d0	
cc       Eq.(21) H02
	 cc=0.0108071d0		
	 cd=0.00301208d0
	 ce=0.000334441d0
	 cf=0.000125d0
	 cpp=0.72d0
	 cq=0.92d0
	 cr=1.2d0
	 ct=1.5d0
	endif
	
	l=alpha*hp		!longueur de melange
	ki=krad/cp/ro		!conductivite thermometrique
	b=2.d0*l**2/9.d0/ki*sqrt(grav/2.d0/hp*delta)		!2 x Eq(6) CM
	! eq.(6) corrected so as to include delta:
	! delta= - (dln rho / ln T)_P                     [RS 30.04.02]
	a=b**2			!c'est le 4a**2 de H02 (Eq. 6)
	
c	initialisations
		
	conv=.false.
	grad=gradad*1.1
	iter=0
	phi=1.d30
	
c	Newton-Raphson "d" signifie derivee/ grad	
	
	do while(.not.conv)
	 iter=iter+1
	 if(iter .gt. 30)then
	  write(6,*)'pas de convergence dans conv_cgm'
	  print*,'donnees : krad,grav,cp,ro,hp,gradrad,gradad'	  
	  write(6,2000)krad,grav,cp,ro,hp,gradrad,gradad
	  print*,'non convergence : phi,phip,abs(phi-phip)/phi'
	  write(6,2000)phi,phip,abs(phi-phip)/phi
	  stop
	 endif
	 phip=phi
	 S=S_sig*a*(grad-gradad)   ! Eq.(6) H02
	 dS=S_sig*a
c	 F1  (Eq. 18 H02) :
	 bS1=1d0+cb*S
	 bS1m=bS1**cm
	 F1=(K0/1.5d0)**3 * ca * (S**ck) * ((bS1m - 1d0)**cn)
	 dF1=F1* (ck/S + cn*cm*cb / (bS1m-1d0) * (bS1m/bS1) )
	 dF1=dF1* dS
c	 F2  (Eq. 20 H02) :
 	 cSp=cc*(S**cpp)
	 dSq1=1d0+cd*(S**cq)
	 eSr= ce*(S**cr)
	 fSt1=1d0+cf*(S**ct)
	 F2=1d0 + cSp/dSq1+eSr/fSt1
	 dF2=cpp* cSp/S/dSq1-cSp*cd*cq*(S**(cq-1d0))/(dSq1**2) 
	 dF2=cr * eSr/S/fSt1-eSr*cf*ct*(S**(ct-1d0))/(fSt1**2) + dF2
	 dF2=dF2* dS
c	 Phi  (Eq. 17 H02) :
	 phi=F1  * F2
	 dphi= F1*dF2+F2*dF1
	 f=grad+phi*(grad-gradad)-gradrad		! Eq. (63) CM
	 df=1.d0+dphi*(grad-gradad)+phi
	 corr=f/df
c	 print*,iter
c	 write(6,2000)corr,abs(phi-phip)/phi
	 grad=grad-corr
	 conv=abs(phi-phip)/phi .le. 1.d-10
	enddo
c	pause
	sig=a*(grad-gradad)	! Eq.(6) H02
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
	 
	 sig=a*(grad-gradad)	! Eq.(6) H02
	 dsig=a
	 dsighp=sig/a*dahp
	 dsigk= sig/a*dak
	 dsigcp=sig/a*dacp
	 dsigro=sig/a*daro	
	 dsiggr=sig/a*dagr
	 dsigad=-a

c	 F1  (Eq. 18 H02) :
	 bS1=1d0+cb*S
	 bS1m=bS1**cm
	 F1=(K0/1.5d0)**3 * ca * (S**ck) * ((bS1m - 1d0)**cn)
	 dF1=F1* (ck/S + cn*cm*cb / (bS1m-1d0) * (bS1m/bS1) )
				! ici dF1 = dF1/dS	 
c	 F2  (Eq. 20 H02) :
 	 cSp=cc*(S**cpp)
	 dSq1=1d0+cd*(S**cq)
	 eSr= ce*(S**cr)
	 fSt1=1d0+cf*(S**ct)
	 F2=1d0 + cSp/dSq1+eSr/fSt1
	 dF2=cpp* cSp/S/dSq1-cSp*cd*cq*(S**(cq-1d0))/(dSq1**2) 
	 dF2=cr * eSr/S/fSt1-eSr*cf*ct*(S**(ct-1d0))/(fSt1**2) + dF2
				! ici dF2 = dF2/dS
c	 Phi  (Eq. 17 H02) :
	 phi=F1  * F2
	 dphi= F1*dF2+F2*dF1
				! ici dphi = dphi/dS

	 dphihp=dphi* dsighp * S_sig
	 dphik =dphi* dsigk * S_sig  
	 dphicp=dphi* dsigcp * S_sig 
	 dphiro=dphi* dsigro * S_sig
	 dphigr=dphi* dsiggr * S_sig
	 dphiad=dphi* dsigad * S_sig 
	 dphi=dphi*dsig*  S_sig    ! dphi = dphi/dgrad

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

