 
c************************************************************************
 
	subroutine etat_opal(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene)
	
c	interface de l'equation d'etat OPAL avec CESAM3	
c	package OPAL_EOS
c	la recherche de ro par un newton, remplace la routine rhoofp
c	du package OPAL_EOS qui ne converge pas dans certains cas solaires
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM3
	
	implicit none
	
	include 'cesam_3.parametres'
	include 'evol_chim_3.common'
        include 'modele_3.common'
 
	integer	mx,mv,nr,nt
	parameter (mx=5,mv=10,nr=121,nt=167)
	
c	iorder=10  gives all 1st and 2nd order data. See instructions in esac.
c	irad=0 does not add radiation  corrections
 
	real*8 p,t,ro,drop,drot,drox,drott,drotp,drotx,xchim(pnelem),
     1	u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene,
     2	stor,stor0,dstor
	
	real*4	x,ztab,t6,r,p4,unpdx,dx
c	data dx/0.001/
 
	character*40 fich(8)
	common/fich_etat/fich
	
	real*4	esact,eos(mv)	
	common/e/esact,eos
	
	character*3 fixedTP
c	data fixedTP/'yes'/  ! 'yes' for fixed T6,P; 'no' for fixed T6,rho
 
	logical deriv,init,ok
c	data	init/.true./
 
	data dx/0.001/
	data fixedTP/'yes'/  ! 'yes' for fixed T6,P; 'no' for fixed T6,rho
	data	init/.true./
 
2000	format(1x,1p8d10.3)
 
c            T6 is the temperature in units of 10**6 K
c            rho is the density in grams/cc
c            R=Rho/T6**3
 
c	on exploite les definitions de OPAL, initialement on a:
 
c            eos(1) is the pressure in megabars (10**12dyne/cm**2)
c            eos(2) is energy in 10**12 ergs/gm. Zero is zero T6
c            eos(3) is the entropy in units of energy/T6
c            eos(4) is dE/dRho at constant T6
c            eos(5) is the specific heat, dE/dT6 at constant V.
c            eos(6) is dlogP/dlogRho at constant T6.
c                   Cox and Guil1 eq 9.82
c            eos(7) is dlogP/dlogT6 at conxtant Rho.
c                   Cox and Guil1 eq 9.81
c            eos(8) is gamma1. Eqs. 9.88 Cox and Guili.
c            eos(9) is gamma2/(gaamma2-1). Eqs. 9.88 Cox and Guili
c            eos(10) is gamma3-1. Eqs 9.88 Cox and Guili
 
c            iorder sets maximum index for eos(i);i.e., iorder=1
c                   gives just the pressure
c
c            irad  if =0 no radiation correction; if =1 adds radiation
 
c            index(i),i=1,10  sets order in which the equation of state
c            variables are stored in eos(i).  Above order corresponds
c            to block data statement:
c                 data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/.
c            If you, for example, only want to return gamma1: set iorder=1
c            and set: data (index(i),i=1,10)/8,2,3,4,5,6,7,1,9,10/
 
c	dans le blockdata blk_eos
c	au lieu de: data(index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/
c	on a mis: data(index(i),i=1,10)/1,4,10,5,6,2,3,7,8,9/
 
c	on aura donc
c            eos(1) is the pressure in megabars (10**12dyne/cm**2)
c            eos(2) is dlogP/dlogRho at constant T6.
c            eos(3) is dlogP/dlogT6 at conxtant Rho.
c            eos(4) is energy in 10**12 ergs/gm. Zero is zero T6
c            eos(5) is dE/dRho at constant T6
c            eos(6) is the specific heat, dE/dT6 at constant V.
c            eos(7) is gamma1. Eqs. 9.88 Cox and Guili.
c            eos(8) is gamma2/(gaamma2-1). Eqs. 9.88 Cox and Guili
c            eos(9) is gamma3-1. Eqs 9.88 Cox and Guili
c            eos(10) is the entropy in units of energy/T6
c	avec, par exemple, iorder=6 on ne calculera que les 6 premieres valeurs
c	d'ou des economies de calculs eventuels
 
	if(init)then
	 init=.false.
	 unpdx=1.+dx
	 write(2,1)
	 write(6,1)
1	 format(//,5x,'-------------Equation d''etat-----------------',//,
     1	5x,'Equation d''etat OPAL si dans [ 5010 , 1.e8 ], EFF sinon',/,
     2	'der. num.; (p,T)-->(ro,T) calcule par ro_new et non rhoofp',/)
	 print*
	endif
	
c	write(6,2000)p,t
 
	p4=p*1.d-12
	x=xchim(1)
	ztab=Z0		!EOS OPAL n'utilise qu'une table en Z
	t6=t*1.d-6
	
	call ro_new(p4,t6,x,ztab,r,ok)
	if(.not.ok)goto10
		       	
	call esac(x,ztab,t6,r,6,1,ok) ! calc EOS; use r from rhoofp
	if(.not.ok)goto10
		
	ro=r
	u=eos(4)*1.e12
	
	drop=r/p/eos(2)		!1 / dp/dro
	drot=-r/t*eos(3)/eos(2)	!- dp/dT / dp/dro
	
	eos(5)=eos(5)*1.e12
	dup=eos(5)*drop			!du/dT dro/dp
	dut=eos(5)*drot+eos(6)*1.e6	!du/dro  dro/dT + du/dt
	
c	derivees / X	
 
	stor0=x
	stor=stor0*unpdx
	if(stor .lt. dx)stor=dx
	dstor=stor-stor0
	x=stor
	call ro_new(p4,t6,x,ztab,r,ok)
	if(.not.ok)goto10
		
        drox=(r-ro)/dstor
	call esac(x,ztab,t6,r,6,1,ok)
	if(.not.ok)goto10
		
	dux=(eos(4)*1.e12-u)/dstor
	x=stor0
	
	if(deriv)then
	 drotx=drot
	 drotx=(-r/t*eos(3)/eos(2)-drotx)/dstor
	
	 dutx=dut	
	 dutx=(-eos(5)*1.e12*r/t*eos(3)/eos(2)+eos(6)*1.e6-dutx)/dstor
 
c	 derivees par rapport a T	
	
	 stor0=t
	 stor=stor0*unpdx
	 dstor=stor-stor0
	 t6=stor*1.e-6
	
	 call ro_new(p4,t6,x,ztab,r,ok)
	 if(.not.ok)goto10 	
	 call esac(x,ztab,t6,r,6,1,ok)
	 if(.not.ok)goto10 	
	 drotp=drop
	 drotp=(r/p/eos(2)-drotp)/dstor			!1 / dp/dro
	
	 drott=drot	
	 drott=(-r/t/unpdx*eos(3)/eos(2)-drott)/dstor	!- dp/dT / dp/dro
	
	 dutp=dup
	 eos(5)=eos(5)*1.e12		
	 dutp=(eos(5)*r/p/eos(2)-dutp)/dstor
	 dutt=dut
	 dutt=(-eos(5)*r/t*eos(3)/eos(2)+eos(6)*1.e6-dutt)/dstor
	 t=stor0
	
	endif
	
	dh1=-100	!taux d'ionisation et degenerescence inconnus
	dhe1=-100	
	dhe2=-100
	degene=5.	
		
	return
	
c	en cas de Pb avec opal appel a EFF	
	
10	call etat_eff_ps_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene)
	
	return
	
	end
 
