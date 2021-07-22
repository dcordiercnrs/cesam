 
c***********************************************************************
 
	subroutine colatm_3(ray,lum,xchim,bp,x,xt,knot,mstar,tdetau,etat,opa)
 
c	resolution du systeme d'equations differentielles non lineaires
c	de l'atmosphere par developpement sur B-splines
c	par iteration Newton avec derivees analytiques
c	adapte de eqdonld
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM: version 3
 
c entree
c	ray : rayon		!au point tau_f
c	lum : luminosite	!au point tau_f
c	xchim : composition chimique par gramme
 
c sortie
c	press: pression	!au point tau_f
c	temp: temperature	!au point tau_f
c	mass: masse		!au point tau_f
c	teff: temperature effective
c	rstar: rayon a tau=2/3
c	bp, x, xt: solution spline
 
c routines externes
c	tdetau, etat, opa
 
c	le reseau de points de raccord x, possede p points
c	la solution, les coefficients des splines d'ordre m+r sur le
c	reseau de points de table xt, sont dans b
 
c	les initialisations, les residus et les derivees sont
c	calcules dans eqatm
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
 
	integer pjac1,pjac2,pxx,pd
	parameter (pjac1=pne_atm*((pn_atm-1)*pm_qs+1), pjac2=(1+pm_qs)*pne_atm,
     1	pxx=(pn_atm-1)*pm_qs, pd=pm_qs+1)
 
	integer ni(pn_atm),knot,spi,ip,i,li,der,spl,indpc(pjac1),
     1	indcol,eq,nl,var,ligne,indice,ipe,vare,compt,cx,rp1
	
	real*8 x(1),xt(1),xchim(1),bp(1),a(pjac1*pjac2),b(pjac1),
     1	xl(pne_atm),xx(pxx),derxx(2*pd*pxx),ray,lum,mstar,
     2	lderxx(2*pd*pne_atm),ae(pne_atm*pne_atm*2),be(pne_atm),
     3	y(pne_atm*2),corr,err,er,dx,to_min,to_max
 
c	data to_min,to_max/0.2d0,0.7d0/
c	real*8 fx(7),dfxdx(7)
 
	logical init,lu
c	data init/.true./
c	data to_min,to_max/0.2d0,0.7d0/
 
	external tdetau,etat,opa
 
	data to_min,to_max/0.2d0,0.7d0/
	data init/.true./
 
	save init,xl,ni,rp1,nl
 
c	write(6,*)'entree dans colatm, tau_min,tau_max,n_atm,m_qs,n23,nea',
c	1	tau_min,tau_max,n_atm,m_qs,n23,nea
 
	if(init)then
	 init=.false.
	 to_min=log(to_min)		!bornes pour tau*
	 to_max=log(to_max)
 
c	 initialisation de la solution et des B-splines
 
	 call eqatm_3(1,knot,li,x,xt,xx,xl,mstar,derxx,lderxx,
     1	bp,xchim,cx,y,be,ae,ray,lum,tdetau,etat,opa,.false.)
 
c	 abscisses des conditions aux limites par ordre croissant
 
	 do li=2,nea
	  if(xl(li) .lt. xl(li-1))then
	   write(6,*)'dans colatm_3 les conditions aux limites',
     1 ' ne sont pas dans l''ordre croissant'
	   stop
	  endif
	 enddo	!li
 
c	 localisation des limites dans la grille de points de raccord
 
	 do ip=1,n_atm
	  ni(ip)=0	!ip
	 enddo
	 do li=1,nea
	  do ip=1,n_atm-1
	   if(x(ip) .le. xl(li) .and. x(ip+1) .gt. xl(li))then
	    ni(ip)=ni(ip)+1
	    go to 10
	   endif
	  enddo	!ip
	  if(x(n_atm) .eq. xl(li))ni(n_atm-1)=ni(n_atm-1)+1
 
c	  write(6,*)'n_atm dans colatm 2',n_atm
 
10	 enddo	!li
 
	 rp1=2
	 nl=nea*((n_atm-1)*m_qs+1)
 
	 dx=min(1.d-7,precix)	!limite inf pour precision relative
c	 dx=1.d-2
 
	endif
 
c"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'
 
c	integration
 
	compt=0
	lu=.false.
 
c	write(6,*)n_atm,nea,m_qs,knot
 
c	do i=1,n_atm
c	 call sbsp1dn(nea,bp,x,xt,n_atm,mpr,knot,.true.,x(i),ip,fx,dfxdx)
c	 write(6,2000)(exp(fx(cx)),cx=1,2),(fx(cx),cx=3,5),(exp(fx(cx)),cx=6,7)
c	enddo
 
20	li=0	!pour les limites
	ligne=0
	do ip=1,n_atm-1
	 spi=(ip-1)*m_qs	!indice-1 de la 1-iere B-spline
	 indcol=nea*spi+1
	 do i=1,m_qs
	  cx=(n_atm-1)*(i-1)+ip
 
c	calcul des coefficients des equations et des seconds membres
c	au point xx, a cause de la non linearite il faut
c	interpoler la solution provisoire pour obtenir la fonction et ses
c	derivees y(variable,derivee) au point xx ( ou xl(li) )
 
c	ae(ne*(ne*(der-1)+var-1)+eq)=coefficient de la (der-1)-ieme derivee de
c	la var-ieme variable de la eq-ieme equation
c	y(variable,derivee)=y(ne*(derivee-1)+variable)
 
	  do var=1,nea
	   do der=1,rp1
	    y(nea*(der-1)+var)=0.
	    do spl=1,mpr
	     y(nea*(der-1)+var)=y(nea*(der-1)+var)+
     1    derxx(rp1*(mpr*(spi+i-1)+spl-1)+der)*bp(nea*(spi+spl-1)+var)
	    enddo	!spl
	   enddo	!der
	  enddo	!var
	  do der=1,nea*nea*rp1	!mise a 0 des derivees
	   ae(der)=0.
	  enddo	!der
 
	  call eqatm_3(2,knot,li,x,xt,xx,xl,mstar,derxx,lderxx,
     1	bp,xchim,cx,y,be,ae,ray,lum,tdetau,etat,opa,lu)
 
	  do eq=1,nea		!eq-ieme equation
	   ligne=ligne+1
	   indpc(ligne)=indcol
	   b(ligne)=be(eq)
	   if(.not.lu)then
	    do spl=1,mpr	!spl - ieme spline
	     do var=1,nea
	      indice=nl*(nea*(spl-1)+var-1)+ligne
	      a(indice)=0.
	      do der=1,rp1 !der-1-eme derivee
	       a(indice)=a(indice)+
     1 ae(nea*(nea*(der-1)+var-1)+eq)*derxx(rp1*
     2 (mpr*(spi+i-1)+spl-1)+der)
	      enddo	!der
	     enddo	!var
	    enddo	!spl
	   endif	!lu
	  enddo		!eq
	 enddo		!i
 
	 if(ni(ip) .eq. 0)goto 30
	 do i=1,ni(ip)
	  li=li+1
	  do var=1,nea
	   do der=1,rp1
	    y(nea*(der-1)+var)=0.
	    do spl=1,mpr
	     y(nea*(der-1)+var)=y(nea*(der-1)+var)+
     1    lderxx(rp1*(mpr*(li-1)+spl-1)+der)*bp(nea*(spi+spl-1)+var)
	    enddo	!spl
	   enddo	!der
	  enddo		!var
	  do der=1,nea*nea*rp1	!mise a 0 des derivees
	   ae(der)=0.
	  enddo	!der
	  call eqatm_3(3,knot,li,x,xt,xx,xl,mstar,derxx,lderxx,
     1	bp,xchim,cx,y,be,ae,ray,lum,tdetau,etat,opa,lu)
 
	  ligne=ligne+1
	  indpc(ligne)=indcol
	  b(ligne)=be(1)
	  if(.not.lu)then
	   do spl=1,mpr
	    do var=1,nea
	     indice=nl*(nea*(spl-1)+var-1)+ligne
	     a(indice)=0.
	     do der=1,rp1    !derivee <r pour cond.lim.
	      a(indice)=a(indice)+
     1    ae(nea*(nea*(der-1)+var-1)+1)*
     2    lderxx(rp1*(mpr*(li-1)+spl-1)+der)
	     enddo	!der
	    enddo	!var
	   enddo	!spl
	  endif		!lu
	 enddo		!i
30	enddo		!ip
 
c	do ligne=1,nl,8
c	 write(6,1000)(b(ligne+i),i=0,7)
c	enddo	!ip
 
c	pause'apres les b'
 
c	do ligne=1,nl
c	 do i=1,mpr
c	  write(6,1000)(a(nl*(nea*(i-1)+var-1)+ligne),var=1,nea)
c	 enddo	!i
c	enddo	!ligne
1000	format((1x,1p10e10.3))
 
	call gausdn(a,b,indpc,nl,nea*mpr,1)
 
c	do spl=1,(p-1)*m_qs+1
c	 write(6,1000)(bp(nea*(spl-1)+var),b(nea*(spl-1)+var),var=1,nea)
c	enddo	!spl
 
c	limitation des corrections
 
	corr=1.
	do spl=1,(n_atm-1)*m_qs+1
	 do var=1,nea
	  indice=nea*(spl-1)+var
c	  write(6,*)indice,corr
40	  if(abs(bp(indice)) .gt. dx .and.
     1 	corr*abs(b(indice)) .gt. .6*abs(bp(indice)))then
	   corr=corr/2.
	   goto40
	  endif
	 enddo	!var
	enddo	!spl
 
c	corrections et b --> bp
 
	err=0.
	do spl=1,(n_atm-1)*m_qs+1
	 do var=1,nea
	  er=0.
	  indice=nea*(spl-1)+var
	  if(abs(bp(indice)) .gt. dx)er=abs(b(indice)/bp(indice))
	  bp(indice)=bp(indice)-b(indice)*corr
	  if(var .eq. 6)bp(indice)=max(to_min,min(to_max,bp(indice)))
	  err=max(er,err)
	  if(er .eq. err)then
	   vare=var
	   ipe=spl/m_qs+1
	  endif
	 enddo	!var
	enddo	!ip
 
c	do i=1,n_atm
c	 call sbsp1dn(nea,bp,x,xt,n_atm,mpr,knot,.true.,x(i),ip,fx,dfxdx)
c	 write(6,2000)(exp(fx(cx)),cx=1,2),(fx(cx),cx=3,5),(exp(fx(cx)),cx=6,7)
c	enddo
c	pause
 
	compt=compt+1
	write(6,*)' '
	write(6,100)compt,err,vare,ipe,corr
100	format(1x,'atmosphere iter.',i3,
     1	' err. max.',1pd8.1,' variable',i2,' couche',i3,' corr',1pd8.1)
c	 write(6,2000)bp(4),exp(bp(6))
c	 call sbsp1dn(nea,bp,x,xt,n_atm,mpr,knot,.true.,x(ipe),ip,fx,dfxdx)
c	 write(6,2000)(exp(fx(cx)),cx=1,2),(fx(cx),cx=3,5),(exp(fx(cx)),cx=6,7)
c	pause
	if(err .lt. dx)then
	 return
	elseif(compt .lt. 30 .or. err .lt. 1.)then
	 lu=err .lt. 5.d-3
	 lu=.false.		!lu existe pour des raisons historiques
	 if(compt .lt. 60)goto 20
	else
	 write(6,*)'pas de convergence dans colatm apres 30 iterations'
	 stop
	endif
 
2000	format((1x,1p8d10.3))
 
	end
 
