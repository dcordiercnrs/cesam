c****************************************************************
 
	subroutine lim_teff_3(list,l_rac,r_rac,xchim,p_rac,t_rac,
     1	dpsdl,dpsdr,dtsdl,dtsdr,t_eff,r_star,m_rac,dmsdl,dmsdr,
     2	p_atm,t_atm,m_atm,tau,r_atm,mstar,
     3	tdetau,etat,opa)
 
c	calcul de l'atmosphere
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM Version 3
 
c entree :
c	list=.true. : calcul reduit pour une liste
c	r_rac : rayon au raccord
c	l_rac : luminosite
c	xchim composition chimique par gramme
c	mstar: masse avec perte de masse
 
c sortie :
c	p_rac : pression au raccord
c	t_rac : temperature
c	m_rac : masse au raccord avec l'enveloppe
c	d*r, d*l : derivees /r /l
c	t_eff : temperature effective
c	r_star : rayon de l'etoile
c	p_atm, t_atm, r_atm, m_atm : pression, temperature, rayon, masse
c	tau : profondeur optique
 
c routines externes :
c	etat : equation d'etat
c	tdetau : loi t(tau,Teff,g)
c	opa : opacite
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
 
	integer pjac1
	parameter (pjac1=pne_atm*(pm_qs*(pn_atm-1)+1))
 
	integer knot,l,i
c	integer j
 
	integer*4 long

	external long

	real*8 r_rac,l_rac,dpsdr,dpsdl,dtsdr,dtsdl,dmsdr,dmsdl,bid,
     1	xchim(pnelem),p_rac,t_rac,m_rac,dr,dl,unpdx,t_eff,
     2	bp(pjac1),p_atm(1),r_atm(1),mstar,mstarp,cte10,
     3	m_atm(1),r_star,cte1,cte2,grav,tau(1),t_atm(1),rd,ld,
     4	x(pn_atm),xt(pn_atm*pm_qs+2),fx(pne_atm),dfxdx(pne_atm)
c	data unpdx/1.00001d0/,mstarp/-1.d0/
	logical ini,list
c	data ini/.true./
 
	character*70 modele_atm
 
	external tdetau,etat,opa
        data ini/.true./
	data unpdx/1.00001d0/,mstarp/-1.d0/
 
	save ini,mstarp,cte1,cte10,cte2
 
2000	format((1x,1p8d10.3))
 
c	write(6,*)'lim_ext,tau_max,n_atm,m_qs,n23',
c	1	tau_max,n_atm,m_qs,n23
 
	if(ini)then
	 ini=.false.
	
	 write(6,*)'------------------ATMOSPHERE-----------------------------'
	 write(2,*)'------------------ATMOSPHERE-----------------------------'
 
	 call tdetau(1.d0,6.d3,1.d4,bid,bid,bid,bid,bid,bid,bid)  !pour tau_ext
	 n_atm=min(n_atm,pn_atm)
 
	 write(2,*)' '
	 write(2,1)tau_max,tau_min,n_atm,n23
1	 format(1x,'limite exterieure avec atmosphere',/,
     1	1x,'tau au fond de l''atmosphere:',1pd10.3,
     2	' rayon et masse pris a tau* tel que T(tau*,Teff)=Teff',/
     3	1x,'tau exterieur:',1pd10.3,', nombre de couches: ',i2,
     4	', R* a la couche: ',i2)
	 write(2,*)' '
	 write(6,*)' '
	 write(6,1)tau_max,tau_min,n_atm,n23
	 write(6,*)' '
	
	 cte10=g*msol/rsol**2
	 cte2=lsol/pi/rsol**2/aradia/clight
 
	endif
 
	if(mstar .ne. mstarp)then
	 mstarp=mstar
	 cte1=cte10*mstar
	endif
 
	write(6,*)' '
	write(6,*)'------- L''atmosphere (debut) ------'
 
c	write(6,*)'avant colatm,tau_max,n_atm,m_qs,n23',
c	1	tau_max,n_atm,m_qs,n23
 
	write(6,10)r_rac,l_rac
10	format(1x,/,1x,'au raccord: R/Rsol=',1pd10.3,', L/Lsol =',1pd10.3,/)
 
	call colatm_3(r_rac,l_rac,xchim,bp,x,xt,knot,mstar,tdetau,etat,opa)
 
	r_star=bp(4)				!bp(q,4) est R_star, tous q
	t_eff=(cte2*l_rac/r_star**2)**(.25)
	grav=cte1/r_star**2
	p_rac=exp(bp(1))
	t_rac=exp(bp(2))
	m_rac=bp(5)
c	write(6,*)'p_rac,t_rac,m_rac'
c	write(6,2000)p_rac,t_rac,m_rac
c	pause'apres colatm'
 
	if(list)then
c	if(.true.)then
	 do i=1,n_atm
	  call sbsp1dn(nea,bp,x,xt,n_atm,mpr,knot,.true.,dfloat(i),
     1	l,fx,dfxdx)
c	  write(6,2000)(exp(fx(j)),j=1,2),(fx(j),j=3,5),(exp(fx(j)),j=6,7)
	  p_atm(i)=exp(fx(1))
	  t_atm(i)=exp(fx(2))
	  r_atm(i)=fx(3)
	  m_atm(i)=fx(5)
	  tau(i)=exp(fx(7))	
c	  write(6,2000)tau(i),p_atm(i),t_atm(i),r_atm(i),m_atm(i)
	 enddo
	 open(unit=31,form='unformatted',status='unknown',
     1	file=nom_fich2(:long(nom_fich2))//'_B.atm')
	 modele_atm=modele(:long(modele))//'(atmosphere)'
	 write(31)bp,x,xt,t_eff,grav,l_rac,xchim,nea,n_atm,
     1	knot,n23,modele_atm
	 close(unit=31)
 
	else		!derivees / r et l
	 write(6,*)' '
	 write(6,*)'derivees partielles numeriques / r et l'
	 rd=r_rac*unpdx
	 dr=rd-r_rac
	 call colatm_3(rd,l_rac,xchim,bp,x,xt,knot,mstar,tdetau,etat,opa)
	 dpsdr=(exp(bp(1))-p_rac)/dr
	 dtsdr=(exp(bp(2))-t_rac)/dr
	 dmsdr=(bp(5)-m_rac)/dr
 
	 ld=l_rac				!derivee /l
	 ld=ld*unpdx
	 dl=ld-l_rac
	 call colatm_3(r_rac,ld,xchim,bp,x,xt,knot,mstar,tdetau,etat,opa)
	 dpsdl=(exp(bp(1))-p_rac)/dl
	 dtsdl=(exp(bp(2))-t_rac)/dl
	 dmsdl=(bp(5)-m_rac)/dl
	endif
 
	write(6,*)' '
	write(6,*)'------- L''atmosphere (fin) ------'
 
	return
 
	end
