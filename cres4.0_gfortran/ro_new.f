c*************************************************************************
 
	subroutine ro_new(p,t6,x,ztab,ro,conv)
 
c	package OPAL_EOS
c	recherche de ro par un newton, remplace la routine rhoofp
c	du package OPAL_EOS qui ne converge pas dans certains cas solaires
 
c	CESAM3
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM3
 
c entrees
c	p: pression
c	t6: temperature/1.d6
c	x: X hydrogene
c	ztab: Z
 
c sorties:
c	ro: densite
 
	implicit none
	
	include'ctephy.common'
	
	integer	mv
	parameter (mv=10)
	
	integer n,nmax
		
	real*4	p,t6,x,ro,lr,lp,granr,mum1,ztab,epsi,
     1	pp,tp,rop,xp
c	data	pp,tp,xp/-1.,-1.,-1./
 
	logical init,conv
c	data init/.true./
 
	real*4	esact,eos(mv)	
	common/e/esact,eos
	
	data	pp,tp,xp/-1.,-1.,-1./
	data init/.true./
 
	save
 
2000	format(1x,1p8d10.3)
 
	if(init)then
	 init=.false.
	 granr=kbol/amu*1.e6*1.e-12	!1.e6 pour T6, 1.e-12 pour p
	 nmax=20
	 epsi=1.e-5		!precision	
	endif
	
c	print*,index
 
c	initialisation de ro
 
	if(max( abs(pp-p)/p,abs(tp-t6)/t6,abs(xp-x) ) .lt. 0.1)then
	 ro=rop
	else
	 mum1=2.*x/ah+3.*(1.-x-ztab)/ahe4	!totalement ionise
	 ro=p/granr/mum1/t6	!ro estime
	endif
	
	lp=log(p)
	lr=log(ro)
	
	n=0
	conv=.false.
	do while((.not.conv .and. n .lt. nmax) .or. n .lt. 1)
	 call esac_dc(x,ztab,t6,ro,2,1,conv)	!calcul
	 if(.not.conv)return
	
	 lr=lr-(log(eos(1))-lp)/eos(2)
	 n=n+1
	
	 conv=abs(eos(1)-p)/p .lt. epsi	
	
	 ro=exp(lr)
	
c	 write(6,2000)p,ro,(eos(1)-p)/p	
c	 print*,conv,n
c	 pause	
	enddo
	
	if(.not.conv)then
	 write(6,10)t6*1.e6,p*1.e12,x
10	 format(1x,/,'!ATTENTION pas de convergence dans ro_new, t=',1pe10.3,
     1	' p=',1pe10.3,' X=',1pe10.3,' ro=',1pe10.3)
	 pp=-1.
	 tp=-1.
	 xp=-1.
	else
	 pp=p
	 tp=t6
	 xp=x
	 rop=ro
	endif	
	
c	write(6,2000)lr
	
	return
	
	end
	
