 
c***********************************************************************
 
	subroutine diffus_3(bp,q,qt,n,knot,chim,mc,mct,nc,knotc,
     1	r2,m23,mstar,chim_t,mc_t,mct_t,nc_t,knotc_t,
     2	dt,convd,convf,nzc,tot_rad,ok,new,	
     3	old_m23,new_m23,nm,new_m23t,knotm,
     4	etat,opa,conv,nuc,coeff_diff)	
	
c	integration de la composition chimique avec diffusion
c	cette routine gere les splines:
c	initialisation: a l'age age=0 la composition chimique est homogene
c	methode des elements finis Petrov-Galerkin
c	on suppose que la perte de masse pour l'element X est X Mdot
c	en faisant porter l'integrale de Petrov-Galerkin sur [0, M(t+dt)]
c	cette perte de masse est automatiquement prise en compte
c	pour une autre valeur de la perte de masse il faudrait integrer
c	sur [0, M(t)] en tenant compte du fait que sur [M(t+dt),M(t)]
c	les abondances sont nulles... donc modifier les produits scalaires.
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM version 3
 
c entrees
c	dt: pas temporel
c	r2,m23,q,qt,knot,n: pour interp. en t+dt
c	mc_t,mct_t,nc_t,knotc_t,chim_t: int. comp. chim. temps t
c	convd,convf,nzc: limites, nombre de ZM
c	mstar: masse au temps t+dt, avec perte de masse
c	tot_rad=.true.: modele totalement convectif
c	new=.true.: initialisation d'une nouvelle solution pour comp.chim
 
c entrees/sorties
c	mc,mct,nc,knotc,chim: int. comp. chim. temps t+dt
 
c sortie
c	ok=.true.: il y a eu convergence
 
c external
c	nuc: reactions thermonucleaires
c	coeff_diff: calcul du coefficient de diffusion
c	etat: equation d'etat
c	opa: opacite
c	conv: routine de convection
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer pae, pa
	parameter (pae=pnelem**2*2, pa=pchim*2*(pm_ch-1)*pnelem)
 
	integer i,j,k,id,nl,iv,ligne,indice,spi,l,ndis,idis(0:pnzc*2+1),
     1	dimf,ie,ig,indpc(pchim),colonne,col,compt,sple,vare,ipe,
     2	mi(pchimt),knot,n,knot_mix,nmix,ld,
     3	nc_t,knotc_t,nc,knotc,bloc,n_fi,knot_fi,
     4	nzc,convd(1),convf(0:*),izc,mfi(pchimt),nm,knotm
 
	real*8	d(pm_ch**2),bes(pnelem),bed(pnelem),mstar,old_m23(1),
     1	qg(pm_ch),wg(pm_ch),corr,err,er,mixd(2*pnzc),new_m23t(1),
     2	dt,mc_t(1),mct_t(1),chim_t(1),m_fi(pchimt),new_m23(1),
     3	dfi(pm_ch**2),mc(1),mct(1),chim(1),bp(1),q(1),qt(1),
     4	mix(pchimt),a_mixt(pchimt),a_mix(pchimt),y(pnelem*2),
     5	aes(pae),aed(pae),mt_fi(pchimt),r2(1),m23(1),a(pa),b(pchim)
	
	logical init,inversible,ok,tot_rad,new,newd
c	data init/.true./
 
	external coeff_diff,etat,opa,conv,nuc
 
	data init/.true./
 
	save init,bloc
 
2000	format(1x,1p8d10.3)
2002	format(1p13d8.1)
 
c	print*,tot_rad,nzc
c	print*,'convd',(convd(i),i=1,nzc+1)
c	print*,'convf',(convf(i),i=0,nzc)
 
c	definitions initialisation
 
	if(init)then
	 init=.false.
 
	 bloc=nbelem*(2*m_ch-1)		!longueur d'un bloc
 
c	 write(6,*)bloc,nbelem,nc,m_ch,tot_rad,nzc
 
	 write(6,*)
	 write(6,*)'Diffusion par methode de Petrov-Galerkin'
 
	 write(2,*)'Diffusion par methode de Petrov-Galerkin'
	 write(6,*)'avec discontinuite de la derivee des abondances aux LMR'
	 write(2,*)'avec discontinuite de la derivee des abondances aux LMR'
	 write(6,*)'ordre des splines:',m_ch
	 write(2,*)'ordre des splines:',m_ch
	 write(2,*)
	 write(6,*)
	 if(mdot .ne. 0.d0)then
	  write(6,21)mdot
	  write(2,21)mdot
21	  format(1x,'Avec un taux de variation externe de masse de',1pd10.3,
     1	' Msol/an')
	 else
	  write(6,*)'sans variation externe de masse'
	  write(2,*)'sans variation externe de masse'
	 endif
	endif		!initialisation
	
c	definition des bases
 
c	a_mix, a_mixt, mix, nmix, knot_mix : vecteur de melange
c base continue Fi :
c	mfi(multiplicite), m_fi(abscisses), mt_fi, n_fi, knot_fi (n_fi .ge. nc)
c base discontinue fi :
c	mi(multiplicite), mc(abscisses), mtc, nc, chim, knotc
c	pour simplifier le calcul des produits scalaires,
c	les points nodaux mtc sont necessairement des points nodaux de Fi
 
	nmix=nc			!dimension du vecteur de melange
	n_fi=nc			!dimension base Fi	
	do i=1,nc
	 a_mix(i)=mc(i)		!vecteur de melange
 
	 m_fi(i)=mc(i)		!abscisses, base Fi
	 mfi(i)=1		!multiplicite, base Fi
	 mi(i)=1		!multiplicite, base fi (peut etre modifiee)
	enddo
	
c	si le modele est tot. radiatif, les deux bases sont identiques
c	s'il est tot. convectif le melange est classique et l'integration
c	est faite dans evol_3
	
c	print*,tot_rad,nzc
c	print*,'convd',(convd(i),i=1,nzc+1)
c	print*,'convf',(convf(i),i=0,nzc)
		
	if(tot_rad)then
	 do i=1,nc
	  mix(i)=1.	!vecteur de melange +1 radiatif, -1 convectif
	 enddo
	 call sbsp1dn(1,mix,a_mix,a_mixt,nmix,2,knot_mix,.false.,a_mix(1),
     1	l,bes,bed)
 
c	 il y a des ZM	
c	 convf(0)=1 si convd(1) .ne. 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc si convf(nzc) .ne. nc, sinon convd(nzc+1)=-10000
	
	else			!tot_rad
	 ndis=0		!nb de discontinuites de la derivee 1iere dans fi
	 do izc=1,nzc+1
	  do i=convf(izc-1),convd(izc)		!dans ZR
	   if(i .eq. convf(izc-1) .and. i .gt. 1)then	!sur une limite ZR/ZC
	    ndis=ndis+1
	    idis(ndis)=i
	    mixd(ndis)=1	!a droite de la limite dans la ZR
	    mi(i)=m_ch-1	!der. 1iere discontinue de fi sur la limite
c	    print*,'ndis,i',ndis,i
c	    pause'debut de ZR a droite de ZM'
	
c	    si m_ch .gt. 2 on ajoute m_ch-2 points a Fi de part et
c	    d'autre des limites ZR/ZM
 
	    if(m_ch .gt. 2)then
	     n_fi=n_fi+1		!+1 pt. dans ZR, (reclasse ensuite)
	     m_fi(n_fi)=(mc(i)+mc(i+1))/2.
	     mfi(n_fi)=1		!multiplicite Fi
	     if(m_ch .eq. 4)then	!2 .le. m_ch .le. 4
	      n_fi=n_fi+1		!+1 pt. dans ZM, (reclasse ensuite)
	      m_fi(n_fi)=(mc(i)+mc(i-1))/2.
	      mfi(n_fi)=1		!multiplicite Fi
 
	     endif
	    endif	!m_ch .gt. 2
	   else
	    mix(i)=1	!dans la ZR
	    mi(i)=1	!multiplicite de fi dans ZR
	   endif	!sur une limite ZR/ZC
	  enddo	!dans ZR
	
	  if(izc .le. nzc)then	!pour les ZM
	   do i=convd(izc),convf(izc)
	    if(i .eq. convd(izc) .and. i .gt. 1)then	!sur une limite
	     ndis=ndis+1
	     idis(ndis)=i
	     mixd(ndis)=-1	!a droite de la limite dans la ZM
	     mi(i)=m_ch-1	!multiplicite de fi sur la limite
c	     print*,'ndis,i,izc',ndis,i,izc
c	     pause'debut de ZM a droite de ZR'
	
c	     si m_ch .gt. 2 on ajoute m_ch-2 points a Fi de part et
c	     d'autre des limites ZR/ZM
 
	     if(m_ch .gt. 2)then
	      n_fi=n_fi+1	!+1 pt. dans ZR, (reclasse ensuite)
	      m_fi(n_fi)=(mc(i)+mc(i-1))/2.
	      mfi(n_fi)=1
	      if(m_ch .eq. 4)then
	       n_fi=n_fi+1	!+1 pt. dans ZM, (reclasse ensuite)
	       m_fi(n_fi)=(mc(i)+mc(i+1))/2.
	       mfi(n_fi)=1
	      endif
	     endif	!m_ch .gt. 2
	    else
	     mix(i)=-1	!dans la ZM
	     mi(i)=1	!multiplicite de fi dans ZM
	    endif	!sur limite ZR/ZC
	   enddo	!dans ZM
	  endif		!pour les ZM
	 enddo	!izc
c	 print*,'ndis,idis,nmix',ndis,(idis(i),i=1,ndis),nmix
	 call sbsp_dis(1,a_mix,mix,ndis,idis,mixd,nmix,2,a_mixt,knot_mix)
	 if(m_ch .gt. 2)call shell(n_fi,m_fi)	!on reordonne Fi
	endif	!tot_rad
	
c	do i=1,nc
c	 call sbsp1dn(1,mix,a_mix,a_mixt,nmix,2,knot_mix,.true.,mc(i),j,bes,bed)
c	 print*,i,mc(i),bes(1)
c	enddo	
	
c	les extremites des vecteurs de multiplicite des bases Fi et fi
 
	mi(1)=m_ch		!pour fi
	mi(nc)=m_ch
	mfi(1)=m_ch		!pour Fi
	mfi(n_fi)=m_ch
 
 
c	print*,'mfi',n_fi
c	print*,(mfi(i),i=1,n_fi)
c	print*,'mi',nc
c	print*,(mi(i),i=1,nc)
 
 
c	print*,'mc',nc
c	write(6,2000)(mc(i),i=1,nc)
c	print*,'m_fi',n_fi
c	write(6,2000)(m_fi(i),i=1,n_fi)
 
c	vecteur nodal pour Fi
 
	call noeud(m_fi,mt_fi,mfi,n_fi,knot_fi)
c	print*,'base Fi',knot_fi,knot_fi-m_ch
c	write(6,2000)(mt_fi(i),i=1,knot_fi)
	
c	vecteur nodal pour fi
	
	call noeud(mc,mct,mi,nc,knotc)
	dimf=knotc-m_ch		!dimension des bases Fi et fi
	nl=nbelem*dimf		!nombre de lignes
c	print*,'base fi',knotc,dimf,knot_fi-m_ch
c	write(6,2000)(mct(i),i=1,knotc)
	j=0
 
	do i=1,nc-1
	 j=j+mi(i)
	enddo
	k=0
	do i=1,n_fi-1
	 k=k+mfi(i)
	enddo
c	print*,'dimension fi, Fi',j,k
c	pause'bases'
	
c	print*,'bloc,nbelem,nchim,nc,m_ch,dimf,knotc,nl'
c	print*, bloc,nbelem,nchim,nc,m_ch,dimf,knotc,nl
c	pause'dimensions'
c	do i=1,n
 
c	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),l,a,b)
c	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
c	1	knotc,.true.,a(5),j,bes,bed)
c	 write(6,*)bes(1),q(i),a(5),l,j
 
c	enddo
c	print*,(mct(i),i=1,knotc)
c	print*,ne,n,mpr,knot,nbelem,nc,m_ch,knotc
c	pause'avant integration'
 
	write(6,*)' '
	write(6,*)'---- integration de l''equation de diffusion (debut) -----'
 
c--------------------iterations Newton Raphson----------------------
 
	compt=0
	newd=new		!pour la premiere initialisation
20	do i=1,nl		!mises a 0
	 b(i)=0.d0
	enddo
	do i=1,nl*bloc		!matrice compressee
	 a(i)=0.d0
	enddo
 
c_____________ le point courant ______________________
 
c	aes(nbelem*(nbelem*(k-1)+j-1)+i)=
c	coefficient de la (k-1)-ieme derivee de la
c	j-ieme variable de la i-ieme equation <.> Ni
c	aed(nbelem*(nbelem*(k-1)+j-1)+i)=
c	coefficient de la (k-1)-ieme derivee de la
c	j-ieme variable de la i-ieme equation <.> d/dx Ni
c	a(<spline 1 . equation 1: variable 1,2,..,nbelem>
c	  <spline 1 . equation 2: variable 1,2,..,nbelem>
c		:	:	:	:	:	
c	  <spline 1 . equation ne: variable 1,2,..,nbelem>
c	  <spline 2 . equation ne: variable 1,2,..,nbelem>
c	  <spline 2 . equation ne: variable 1,2,..,nbelem>
c		:	:	:	:	:
c	  <spline 2 . equation ne: variable 1,2,..,nbelem>
c	  <spline 3 . equation ne: variable 1,2,..,nbelem>
c		:	:	:	:	:
c	  <spline dimf . equation ne: variable 1,2,..,nbelem>)
 
c	bes(i)=second membre de la i-ieme equation <.> Ni
c	bed(i)=second membre de la i-ieme equation <.> d/dx Ni
 
c	b(<spline 1 . equation 1, 2,.., nbelem>
c	  <spline 2 . equation 1, 2,.., nbelem>
c		:	:	:	:
c	  <spline dims . equation 1, 2,.., nbelem>)
 
c	indices de premiere colonne pour les produits scalaires
 
	do i=1,dimf		!pour chaque spline Ni
	 k=max(1,i-m_ch+1)	!indice s'il n'y avait qu'une variable
c	 write(6,*)'i,k',i,k
	 do ie=1,nbelem		!pour chaque equation
	  indpc(nbelem*(i-1)+ie)=nbelem*(k-1)+1
	 enddo
	enddo
c	write(6,*)'indices de premiere colonne'
c	write(6,*)(indpc(i),i=1,nbelem*dimf)
 
c	formation des produits scalaires inspire de l'algorithme 5.22 p. 203
c	de Schumaker
 
	do k=m_ch,dimf	!le vecteur nodal des Fi n'a que des points simples
	 spi=k-m_ch		!indice -1 de la premiere spline non nulle
	 call intgauss(mt_fi(k),mt_fi(k+1),qg,wg,m_ch) !Xi,Wi car int.Gauss
	 do ig=1,m_ch
	  call slinf(qg(ig),mct,knotc,l)
	  call bvald(qg(ig),mct,m_ch,l,2,d)	!d(der,spline)
	  if(newd)then			!initialisation
	   call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,knotc_t,.true.,
     1	qg(ig),ld,bes,bed)
	   do iv=1,nbelem
	    y(iv)=bes(iv)
	    y(iv+nbelem)=bed(iv)
	   enddo
 
	  else	
	   do iv=1,nbelem
	    do id=1,2
	     y(nbelem*(id-1)+iv)=0.
	     do j=1,m_ch
	      y(nbelem*(id-1)+iv)=y(nbelem*(id-1)+iv)+d(m_ch*(j-1)+id)*
     1     chim(nbelem*(spi+j-1)+iv)	!chim(variable,spline)
	     enddo	!j
	    enddo	!id
	   enddo	!iv
	  endif
 
c	  appel a l'eq. de diffusion
 
	  do i=1,pae	!mises a 0
	   aes(i)=0.d0
	   aed(i)=0.d0
	  enddo
	  do i=1,nbelem
	   bes(i)=0.d0
	   bed(i)=0.d0
	  enddo
c	  print*,k,pae,nbelem,y(1)
	  call eq_diffus_3(y,qg(ig),aes,aed,bes,bed,1,mstar,
     1	bp,q,qt,knot,n,r2,m23,mc_t,mct_t,nc_t,knotc_t,chim_t,dt,
     2	a_mix,nmix,mix,a_mixt,knot_mix,
     3	old_m23,new_m23,nm,new_m23t,knotm,
     4	etat,opa,conv,nuc,coeff_diff)
	  call bvald(qg(ig),mt_fi,m_ch,k,2,dfi)	!dd(der,spline)
	  do ie=1,nbelem		!pour chaque equation
	   do i=1,m_ch			!pour chaque spline i
	    ligne=nbelem*(spi+i-1)+ie
	    b(ligne)=b(ligne)+wg(ig)*(dfi(m_ch*(i-1)+1)*bes(ie)+
     1		dfi(m_ch*(i-1)+2)*bed(ie))
c	    write(6,2003)ligne,dfi(m_ch*(i-1)+1),dfi(m_ch*(i-1)+2),bes(ie),
c	1	bed(ie),y(nbelem*(1-1)+1),y(nbelem*(2-1)+1),qg(ig)
2003	    format(1x,i3,1p9d8.1)
	    do j=1,m_ch			!pour chaque spline j
	     do iv=1,nbelem		!pour chaque variable
	      colonne=nbelem*(spi+j-1)+iv
	      col=colonne-indpc(ligne)+1!colonne dans matrice compressee
c	      indice=nl*(colonne-1)+ligne	!indice matrice normale
	      indice=nl*(col-1)+ligne		!indice matrice compressee
c	      write(6,*)'ligne,colonne,col',ligne,colonne,col
	      do id=1,2	!pour 1: NiNj, 2: NiN'j, 3: NiN"j etc...
	       a(indice)=a(indice)+wg(ig)*(
     1		aes(nbelem*(nbelem*(id-1)+iv-1)+ie)*
     2		dfi(m_ch*(i-1)+1)*d(m_ch*(j-1)+id)+
     3		aed(nbelem*(nbelem*(id-1)+iv-1)+ie)*
     4	 	dfi(m_ch*(i-1)+2)*d(m_ch*(j-1)+id))
	      enddo	!id
	     enddo	!iv variable
	    enddo	!j
	   enddo	!i
	  enddo	!ie equation
	 enddo 	!ig
 
c	 appel a l'eq. de diffusion pour la perte de masse
c	 supprime car l'integration Petrov-Galerkin porte sur [0 M(t+dt)]
 
	 if(mdot .gt. 0.d0 .and. k .eq. dimf .and. .false.)then !perte de masse
	  call slinf(mc(nc),mct,knotc,l)
	  call bvald(mc(nc),mct,m_ch,l,2,d)	!dd(der,spline)
	  do iv=1,nbelem		!on ajoute X_i Mdot dt aux derniers
	   y(nbelem*(1-1)+iv)=0.	!produits scalaires
	   do j=1,m_ch
	    y(nbelem*(1-1)+iv)=y(nbelem*(1-1)+iv)+d(m_ch*(j-1)+1)*
     1   chim(nbelem*(spi+j-1)+iv)	!chim(variable,spline)
	   enddo	!j
	  enddo	!iv
	  do i=1,pae	!mises a 0
	   aes(i)=0.d0
	  enddo
	  do i=1,nbelem
	   bes(i)=0.d0
	  enddo
c	  print*,k,pae,nbelem,y(1)
 
	  call eq_diffus_3(y,qg(ig),aes,aed,bes,bed,2,mstar,
     1	bp,q,qt,knot,n,r2,m23,mc_t,mct_t,nc_t,knotc_t,chim_t,dt,
     2	a_mix,nmix,mix,a_mixt,knot_mix,
     3	old_m23,new_m23,nm,new_m23t,knotm,	
     4	etat,opa,conv,nuc,coeff_diff)
	
	  do ie=1,nbelem		!pour chaque equation
	   ligne=nbelem*(spi+m_ch-1)+ie	!spline m_ch = 1 seule, la derniere
	   b(ligne)=b(ligne)+bes(ie)
	   colonne=nbelem*(spi+m_ch-1)+ie	!spline m_ch = 1
	   col=colonne-indpc(ligne)+1	!colonne dans matrice compressee
c	   indice=nl*(colonne-1)+ligne	!indice matrice normale
	   indice=nl*(col-1)+ligne		!indice matrice compressee
c	   write(6,*)'ligne,colonne,col',ligne,colonne,col
	   a(indice)=a(indice)+aes(nbelem*(nbelem*(1-1)+ie-1)+ie)
	  enddo	!ie equation
	 endif		!k=dimf
	enddo		!k
 
c	do i=1,dimf
c	 write(6,1000)(a(nl*(j-1)+i),j=1,bloc),b(i)
c	enddo
c	pause' a et b'
 
c	write(6,*)nbelem,dimf,nl,ligne,nbelem*dimf,bloc
c	print*,(indpc(i),i=1,nl)
c	pause'indpc0'
c	do li=1,2		!nl
c	 print*,li
c	 write(6,1000)(a(nl*(i-1)+li),i=1,bloc),b(li)
c	 pause
c	enddo
1000	format(1p10d8.1)
c	pause'a,b'
 
c____________ solution _____________________________
 
	call gausdp_g(a,b,indpc,nl,nbelem*dimf,bloc,inversible)
c	print*,(indpc(i),i=1,nl)
c	pause'indpc1'
	if(.not.inversible)stop
 
c	write(6,*)'inversible/les b',inversible
c	write(6,2000)(b(i),i=1,dimf*nbelem)
c	write(6,*)(indpc(i),i=1,nl)
c	write(6,2000)(ab_min(i),i=1,nbelem)
c	pause'indpc2'
 
c	limitation des corrections
 
	corr=1.
	do j=1,dimf		!pour chaque spline
	 do iv=1,nbelem
	  indice=nbelem*(j-1)+iv
c	  write(6,*)indice,corr
	  do while(abs(chim(indice)) .gt. ab_min(iv) .and.
     1 	corr*abs(b(indice)) .gt. .6*abs(chim(indice))
     2	.and. (iv .eq. 1 .or. iv .eq. ihe4))
	   corr=corr/2.
	  enddo
	 enddo	!iv
	enddo	!j
 
c	corrections et b --> chim
 
	err=0.
	do j=1,dimf
	 do iv=1,nbelem
	  er=0.
	  indice=nbelem*(j-1)+iv
	  if(chim(indice) .gt. ab_min(iv) .and.
     1	(iv .eq. 1 .or. iv .eq. ihe4))
     2	er=abs(b(indice)/chim(indice))
	  chim(indice)=chim(indice)-b(indice)*corr
	  chim(indice)=max(chim(indice),ab_min(iv)*1.d-5)	!ab>0
	  err=max(er,err)
	  if(er .eq. err)then
	   sple=j
	   vare=iv
	  endif
c	  print*,j,iv,er,b(indice),chim(indice)
	 enddo	!iv
	enddo	!j
c	pause
 
	if(iw .gt. 1)chim(iw)=0.	!MW=0 au centre
	
c	do i=1,nc		!nc
c	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,knotc,.true.,mc(i),k,bes,bed)
c	 print*,i,mc(i)
c	 write(6,2001)(bes(j),j=1,nbelem)
c	 write(6,2001)(bed(j),j=1,nbelem)
c	 print*,mc(i),bes(1)
2001	 format(1x,1p10d8.1)
c	enddo
c	pause'solution1'
 
	ipe=sple/m_ch
	compt=compt+1
	write(6,*)' '
	write(6,100)compt,err,vare,sple,ipe,corr
100	format(1x,'iteration',i3,' err. max.',1pd8.1,' variable',i2,
     1' B-spline',i4,' couche',i4,' corr',1pd8.1)
	
c	do i=1,n
c	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),l,a,b)
c	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
c	1	knotc,.true.,min(a(5),mc(nc)),j,bes,bed)
c	 write(6,*)mc(1),mc(nc),a(5),bes(1)
c	enddo
c	print*,ne,n,mpr,knot,nbelem,nc,m_ch,knotc
c	pause'solution2'
	
	newd=.false.		!on utilisera la nouvelle solution
	ok=err .lt. precix
	if(ok)then
	 write(6,*)' '
	 write(6,*)'----- integration de l''equation de diffusion (fin) -------'
	elseif(compt .le. 8)then
	 goto 20
	else
	 write(6,*)'pas de conv. dans diffus apres 9 iterations'
	endif
	
c	pause'fin de diffus_3'
 
	return
 
	end
 
