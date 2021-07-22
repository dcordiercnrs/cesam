 
c***********************************************************************
 
	blockdata blk_eos
	
c	adaptation a CESAM3 du blockdata  du package OPAL_EOS
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
	implicit none
	
	integer mx,mv,nr,nt
	parameter (mx=5,mv=10,nr=121,nt=167)
	
	integer	i
 
	real*4	q(4),h(4),xxh
	common/aa/q,h,xxh
 
	integer	m,mf
	real*4	xz(mx,mv,nt,nr),t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),
     1	esk2(nt,nr),dfsx(mx),dfs(nt),dfsr(nr),xa(mx)
	common/a/xz,t6list,rho,t6a,esk,esk2,dfsx,dfs,dfsr,m,mf,xa
	
	integer iri(10),index(10),nta(nr)
	real*4	 zz(mx)
	common/b/iri,index,nta,zz
 
	data (xa(i),i=1,mx)/0.0,0.2,0.4,0.6,0.8/
	
	data (nta(i),i=1,nr)/73*167,14*159,
     1	2*99,91,87,79,75,67,67,63,59,55,53,51,51,49,47,45,45,41,41,39,
     2	37,35,33,31,29,29,25,25,21,21,17,17,13/
 
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
 
c            index(i),i=1,10  sets order in which the equation of state
c            variables are stored in eos(i).  Above order corresponds
c            to block data statement:
c                 data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/.
c            If you, for example, only want to return gamma1: set iorder=1
c            and set: data (index(i),i=1,10)/8,2,3,4,5,6,7,1,9,10/
c	c'est a dire que la 1iere variable est en 8ieme position
c	et la 8ieme variable est en premiere position
 
c	on deplace l'entropie de l'indice 3 a l'indice 10
c	et mis les derivees de ln P en rang 2 et 3
c	au lieu de: data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/
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
 
c	data(index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/
	data index /1,4,10,5,6,2,3,7,8,9/	
	
	end
	
