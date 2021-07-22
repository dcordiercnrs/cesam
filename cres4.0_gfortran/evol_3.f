 
c********************************************************************
 
	subroutine evol_3(new,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,r2_,m23_,dt_,
     4	bp__,q__,n__,qt__,knot__,chim__,mc__,nc__,mct__,knotc__,
     5	r2__,m23__,dt__,kmax,estim,g_max,
     6	old_m23,new_m23,nm,new_m23t,knotm,
     7	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     8	old_m23__,new_m23__,nm__,new_m23t__,knotm__,
     9	etat,opa,conv,nuc,coeff_diff)	
	
c	gestion de l'evolution temporelle de la composition chimique
c	les points d'integration en comp. chim. sont les points de raccord
 
c	les ZC sont melangees
c	pour le moment angulaire Mw, si iw > 1 il y a melange de Mw
c	dans les ZM; on integre w r**2 au temps t et on divise par
c	(somme dm au temps t+dt)
c	dans ZR on conserve Mw = w r**2
 
c	ensuite on reprend la solution avec ou sans diffusion;
c	dans ZM on la multipliera Mw (moyen) par
c	(somme dm temps t+dt) / (somme r**2 dm temps t+dt) et on multipliera
c	par r**2; on reconstruit le vecteur de composition chimique avec
c	discontinuite des xchim aux LMR.
c	ainsi w, en tout point, sera obtenu en divisant Mw interpole, par r**2
 
c	Auteur:
c	P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	23 09 96 suppression du mouvement des ZR/ZC avec diffusion
c	09 10 96 modif de la gestion de somme dm/somme r2 dm dans ZC
 
c	version 3
 
c entree
c	new=.true.: premiere integration t-->dt, .false. : ameliorations
c	age: age initial et final
c	n: nombre de points
c	lim: nombre de limites ZC/ZR temps t+dt
c	jlim: indices des limites ZC/ZR temps t+dt
c	mc_t,mct_t,nc_t,knotc_t,chim_t: pour int. comp.chim. temps t
c	bp_t,q_t,qt_t,knot_t,n_t  var. pples. temps t
c	m_zc, r_zc, r_ov: masses, rayons aux limites ZR/ZM en Msol, Rsol
c	lconv: true/false debut d'une ZC/ZR aux temps t+dt
c	n: nombre de points
c	dt: pas temporel (fixe dans UPDATE)
c	bp,q,qt,knot: solution temps t+dt
c	dt_,mc_,mct_,nc_,knotc_,chim_: pour int. comp.chim. temps t-1
c	mstar: masse avec perte de masse
c	old_m23,new_m23,nm,new_m23t,knotm : interpolation m(t+dt)-->m(t)
 
c sortie
c	g_max: il faut reduire le pas temporel a cause de TdS ou de non CV
c	estim : estimation de l'erreur max
c	dts : estimation du pas temporel a utiliser pour le pas suivant
c	mc,mct,nc,knotc,chim comp. chim. temps t+dt
c	kmax : indice de couche de l'erreur max
 
c external
c	etat, opa, conv, nuc
 
c	modele totalement convectif: lim=1,jlim(1)=n,lconv(1)=.false.
c	modele totalement radiatif: lim=0,jlim(i)=-100,lconv(i)=.false.
c	les m_zc sont en m/Mstar
c	chim(i,iw)= moment angulaire par unite de masse
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer	n,knot,nc,knotc,lim,jlim(1),idis(0:2*pnzc+1),ndis,
     1	n_t,knot_t,nc_t,knotc_t,lim_t,jlim_t(1),n_zr,n_zc,
     2	n_,knot_,nc_,knotc_,n__,knot__,nc__,knotc__,
     3	i,j,k,ll,kmax,kmaxi(pnelem),izc,nzc,
     4	convd(pnzc+1),convf(0:pnzc),nadd,jfin,it,
     5	nm,nm_,nm__,knotm,knotm_,knotm__
			
	real*8	bp(1),q(1),qt(1),chim(1),mc(1),mct(1),m_zc23(2*pnzc),
     1	r2(1),m23(1),dts,m_zc(1),r_zc(1),r_ov(1),mstar,age,
     2	bp_t(1),q_t(1),qt_t(1),chim_t(1),mc_t(1),mct_t(1),
     3	r2_t(1),m23_t(1),dt,m_zc_t(1),r_zc_t(1),r_ov_t(1),
     4	bp_(1),q_(1),qt_(1),chim_(1),mc_(1),mct_(1),
     5	r2_(1),m23_(1),dt_,bp__(1),q__(1),qt__(1),chim__(1),
     6	mc__(1),mct__(1),r2__(1),m23__(1),dt__,estim(1),
     7	old_m23(1),new_m23(1),new_m23t(1),nuc_m(pnelem),
     8	old_m23_(1),new_m23_(1),new_m23t_(1),	
     9	old_m23__(1),new_m23__(1),new_m23t__(1)	
	
	real*8	esti(pnelem),est,dtnew,dtn,dm(pnch),f(pne),dfdx(pne),
     1	compx(pnch*pnelem),compy(pnch*pnelem),drott,drotp,pas,mc_max,
     3	p,t(pnch),ro(pnch),drop,drot,drotx,dutt,dutp,dutx,mk,
     4	bidd(pnelem),drox,u,dup,dut,dux,p_t,t_t(pnch),bid,
     5	ro_t(pnch),nh1,nhe1,nhe2,compv(pchim),
     6	g(pnzc*pnelem),chimd(2*pnzc*pnelem),m_zc23t,mk_t,
     7	comp1(pnelem),comp2(pnelem),mc_maxi(pnelem),lt_inf,
     8	ab_max(pnelem),lamb,mass_zc,pas_min,mstar23,dmsr2dm(2*pnzc),
     9	p_,t_(pnch),ro_(pnch),p__,t__(pnch),ro__(pnch)
	
	logical init,new,g_max,ok,tot_conv,tot_rad,lconv(1),lconv_t(1),
     1	demi,core_reg,zc_aug,tot_conv_t
c	data init/.true./,demi/.true./
 
	external etat,opa,conv,nuc,coeff_diff
 
	data init/.true./,demi/.true./
 
	save
 
2000	format((1x,1p8d10.3))
2001	format(1x,1p1d10.3,(1p3d22.15))
 
c	write(6,*)'n,nbelem,nchim,nc,knotc',n,nbelem,nchim,nc,knotc
c	write(6,*)'dt,dt_,t_inf,lim,jlim,l_conv',dt,dt_,t_inf,lim,
c	1	(jlim(i),i=1,lim),(lconv(i),i=1,lim)
c	pause'entree evol_3'
	
	if(init)then
	 init=.false.
	 lt_inf=log(t_inf)
	 pas_min=1.d-6	!distance minimale entre deux couches
	 n_zr=10	!nb. de points minimum dans ZR ou ZC
	 n_zc=10
	endif
	mc_max=(1.d0-1.d-5)*mstar  !au dela on supprime les ZR/ZC	
 
	if(dt .lt. dtmin)then
	 write(6,*)'EVOL_3(1) dt trop petit abandon/dt,dtmin'
	 write(6,2000)dt,dtmin
	 !pause'ABANDON'
	 stop
	endif
	
c	masse externe	
	
	mstar23=min(mstar**(2.d0/3.d0),m23(n))
 
c	modele totalement convectif: lim=1,jlim(1)=n,lconv(1)=.false.
c	nzc=1, convd(1)=1, convf(1)=nc
c	modele totalement radiatif: lim=0,jlim(i)=-100,lconv(i)=.false.
c	nzc=0, convd(1)=nc, convf(0)=1
 
	tot_conv=lim .eq. 1 .and. (jlim(1) .eq. n) .and. .not.lconv(1)
	tot_rad=lim .eq. 0
	
c	print*,tot_conv,lim,jlim(1),.not.lconv(1)
c	pause'tot_conv'
 
c	elimination des limites ZR/ZC trop externes
 
	if(.not.tot_conv .and. .not.tot_rad)then
	 i=1
	 do while(i .le. lim)
c	  print*,i,m_zc(i),mc_max
	  if(m_zc(i) .gt. mc_max)then
	   print*	!,i,lim,lconv(i)
	   write(6,*)'Pour le melange convectif elimination des limites',
     1	' ZR/ZC trop externes #', (j,j=i,lim)
 
	   if(lim .eq. 1 .and. .not.lconv(1))then 	!fin de l'unique ZC
	    jlim(1)=n		!que l'on deplace a la fin du modele
	    write(6,*)'Modele completement melange'
	    tot_conv=.true.	!qui devient totalement convectif
	   else
	    lim=max(i-1,0)
c	    print*,'lim, jlim, lconv',lim,(jlim(i),lconv(i),i=1,lim)
c	    write(6,2000)(m_zc(i),i=1,lim)
	    tot_rad=lim .le. 0
	    if(tot_rad)then
	     jlim(i)=-100
	     lconv(i)=.false.
	     write(6,*)'Modele sans zone melangee'
	    endif
	   endif
	  endif
	  i=i+1
	 enddo
c	 pause'apres les eliminations de limites ZR/ZC'
	endif
 
c	les masses (m**2/3) pour la comp chim
c	on prend la repartition deduite de bp, puis on ajoute les limites
c	ZR/ZC, enfin on s'assure que les ZC et ZR retenues ont, au moins,
c	n_rad et n_zc points
	
	do i=1,lim
	 m_zc23(i)=m_zc(i)**(2.d0/3.d0)
	enddo
 
	do i=1,n
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),ll,f,dfdx)
	 mc(i)=f(5)
	enddo
	nc=n+1
	mc(nc)=mstar23

c -------------------------------------------------------
c	un pas temporel sur 2 on prend les mc au milieu
 
c	if(new)demi=.not.demi
c	if(demi)then
c	 do i=nc,2,-1
c	  mc(i)=(mc(i)+mc(i-1))/2.
c	 enddo
c	 nc=nc+1
c	 mc(nc)=mstar23
c	endif
c -------------------------------------------------------

c	on ajoute les m_zc au temps t+dt et du temps t	
	
	do i=1,lim
	 mc(i+nc)=m_zc23(i)
	enddo
	nc=nc+lim
	
	do i=1,lim_t
	 mc(i+nc)=m_zc_t(i)**(2.d0/3.d0)
	enddo
	nc=nc+lim_t
	
c	tri, suppression des doubles et limites interne et externe
	
	call shell(nc,mc)		!tri
c	print*,nc
c	print*,(mc(i),i=1,nc)
c	pause'apres le tri'
 
	k=1			!suppression des mc < 0
	do i=2,nc
	 if(mc(i) .gt. 0.)then
	  k=k+1
	  mc(k)=mc(i)
	 endif
	enddo
	nc=k
	mc(1)=0.
	
	k=1			!suppression des doubles
	do i=2,nc
	 if(mc(k) .ne. mc(i))then
	  k=k+1
	  mc(k)=mc(i)
	 endif
	enddo
	nc=k
	
	do i=k,1,-1		!limite externe a m23=mstar23
	 if(mc(i) .ge. mstar23)nc=i
	enddo
	mc(nc)=mstar23
	
c	elimination des couches plus serrees que pas_min
 
	do i=2,nc-1
	 if(mc(i)-mc(i-1) .lt. pas_min)then
	  do j=1,lim
	   if(mc(i) .eq. m_zc23(j))mc(i-1)=mc(i)	!on gardera les m_zc23
	  enddo
	  mc(i)=mc(i-1)	!on double
	 endif
	enddo	
	k=1			!suppression des doubles
	do i=2,nc
	 if(mc(k) .ne. mc(i))then
	  k=k+1
	  mc(k)=mc(i)
	 endif
	enddo
	nc=k
	
	print*
	print*,'nombre de couches pour la comp. chim.',nc
	if(demi)print*,'couches aux points semi entiers'
	
c	print*,nc
c	print*,(mc(i),i=1,nc)
c	pause'apres limite externe'
	
c	write(6,*)'EVOL_3 limites, n, lim',n,lim
c	do i=1,lim
c	 write(6,*)i,jlim(i),lconv(i),m_zc(i)
c	enddo
 
	est=0.			!initialisations precision
	do i=1,nbelem
	 estim(i)=0.d0
	enddo
 
c	delimitation des zones convectives
c	il y a nzc zones convectives chacune entre convd(.) et convf(.)
 
	do i=1,pnzc		!initialisations
	 convd(i)=-100
	 convf(i)=-100
	enddo
 
	if(tot_conv)then
	 nzc=1
	 convd(1)=1
	 convf(1)=nc
	elseif(tot_rad)then
	 nzc=0
	 convd(1)=nc
	 convf(0)=1
	else
	 nzc=0		!nombre de ZC
	 ll=1		!indice des limites
	 do i=1,nc
c	  print*,i,ll,nzc,lconv(ll),m_zc23(ll),mc(i)
	  if(m_zc23(ll) .eq. mc(i))then
	   if(lconv(ll))then		!debut de ZC
	    nzc=nzc+1			!une ZC en plus
	    convd(nzc)=i		!indice du debut de la ZC a droite
 
	   else				!fin de ZC
	    if(nzc .eq. 0)then		!la zone externe est ZC
	     nzc=1
	     convd(nzc)=1
	    endif
	    convf(nzc)=i		!indice de fin de la ZC a gauche
	   endif
	   ll=min(lim,ll+1)		!limite suivante
	  endif
	 enddo	!i, pas de limite en n
	 if(nzc .ge. 1 .and. convf(nzc) .lt. 0)convf(nzc)=nc
	endif
	
c-------------------------------------------------------------------------
 
	if(.false.)then
c	if(.true.)then
	 print*,lim,nzc,(convd(i),convf(i),i=1,nzc)
	 do i=1,nzc
	  write(6,*)i,convd(i),convf(i),mc(convd(i)),mc(convf(i))
	 enddo	!i
	 write(6,*)' '
	 j=1
	 do i=1,nc
	  if(i .eq. convd(j))then	
	   if(i .ne. 1)write(6,*)mc(i),i
	   write(6,*)'debut de ZC',j,i
	   write(6,*)mc(i),i
	  elseif(i .eq. convf(j))then
	   write(6,*)mc(i),i
	   write(6,*)'fin de ZC',j,i
	   if(i .ne. nc)write(6,*)mc(i),i	
	   j=min(j+1,nzc)
	  else
	   write(6,*)mc(i),i
	  endif
	 enddo
	 pause'avant ajustement des ZR et ZC'
	endif
c------------------------------------------------------------------------
 
c	les ZR doivent avoir au moins m_ch+1 couches
 
c	print*,(mc(i),i=1,nc)
c	pause'repartition initiale'	
 
	if(nzc .gt. 0)then
	 if(convd(1) .gt. 1)then
	  nadd=n_zr+1-convd(1)	!pour une ZR centrale
	  if(nadd .gt. 0)then
	   print*
	   write(6,*)'adjonction de',nadd,' couches pour la ZR centrale'
	   if(nadd .gt. n_zr+1)then
	    print*,'Dans evol_3 (1) : erreur nadd > n_zr+1'
	    !pause'ABANDON'
	    stop
	   endif
	   do i=nc,convd(1),-1		!decalages
	    mc(i+nadd)=mc(i)
	   enddo
	   do i=1,nzc
	    convd(i)=convd(i)+nadd
	    convf(i)=convf(i)+nadd
	   enddo
	   nc=nc+nadd	!adjonction de masses
	   pas=(mc(convd(1))-mc(1))/(n_zr+1)
	   do i=2,n_zr+1
	    mc(i)=mc(1)+pas*i
	   enddo
	  endif
	 endif
 
	 do k=1,nzc-1	!pour les ZR internes
	  nadd=n_zr+1+convf(k)-convd(k+1)
	  if(nadd .gt. 0)then
	   print*
	   write(6,*)' + ',nadd,' couches pour la ZR qui suit la ZC',k
	   if(nadd .gt. n_zr+1)then
	    print*,'Dans evol_3 (2) erreur nadd > n_zr+1'
	    !pause'ABANDON'
	    stop
	   endif
	   do i=nc,convd(k+1),-1	!decalages
	    mc(i+nadd)=mc(i)
	   enddo
	   do i=k+1,nzc
	    convd(i)=convd(i)+nadd
	    convf(i)=convf(i)+nadd
	   enddo
	   nc=nc+nadd	!adjonction de masses
	   pas=(mc(convd(k+1))-mc(convf(k)))/(n_zr+1)
	   do i=1,n_zr
	    mc(convf(k)+i)=mc(convf(k))+pas*i
	   enddo
	  endif
	 enddo
 
	 nadd=0
	 if(convf(nzc) .ne. nc)nadd=n_zr+1+convf(nzc)-nc
	 if(nadd .gt. 0)then	!pour une ZR externe
	  print*
	  write(6,*)'adjonction de',nadd,' couches pour une ZR externe'
	  if(nadd .gt. n_zr+1)then
	   print*,'Dans evol_3 (3) erreur nadd > n_zr+1'
	   !pause'ABANDON'
	   stop
	  endif
c	  print*,'avant',nc,k,mc(convf(nzc)),convf(nzc)
	  pas=(mc(nc)-mc(convf(nzc)))/(n_zr+1)
	  nc=nc+nadd
	  mc(nc)=mstar23	
c	  print*,'apres',nc
	  do i=1,n_zr
	   mc(convf(nzc)+i)=mc(convf(nzc))+pas*i
c	   print*,convf(nzc)+i,i,mc(convf(nzc)+i),mc(convf(nzc)),pas*i
	  enddo
c	  print*,(mc(i),i=convf(nzc),nc)
c	  pause'nouvelle repartition'
	 endif
	
c	 les ZC doivent avoir aussi plus de n_zc+1 couches
 
	 do izc=1,nzc	!pour les ZC
	  nadd=n_zc+1+convd(izc)-convf(izc)
	  if(nadd .gt. 0)then
	   print*
	   print*,'adjonction de',nadd,' couches dans la ZC',izc
	   if(nadd .gt. n_zc+1)then
	    print*,'Dans evol_3 (4) erreur nadd > n_zc+1'
	    !pause'ABANDON'
	    stop
	   endif
c	   print*,'avant',nc,izc,mc(convd(izc)),convd(izc),
c	1	mc(convf(izc)),convf(izc)
	   do i=nc,convf(izc),-1	!decalages
	    mc(i+nadd)=mc(i)
	   enddo
	   convf(izc)=convf(izc)+nadd
	   do i=izc+1,nzc
	    convd(i)=convd(i)+nadd
	    convf(i)=convf(i)+nadd
	   enddo
	   nc=nc+nadd	!adjonction de masses	
c	   print*,'apres',nc,izc,mc(convd(izc)),convd(izc),
c	1	mc(convf(izc)),convf(izc)	
	   pas=(mc(convf(izc))-mc(convd(izc)))/(n_zc+1)
	   do i=1,n_zc
	    mc(convd(izc)+i)=mc(convd(izc))+pas*i
c	    print*,convd(izc)+i,i,izc,mc(convd(izc)+i),mc(convd(izc)),pas*i
	   enddo
c	   print*,convd(izc),convf(izc)
c	   print*,(mc(i),i=1,nc)
c	   pause'nouvelle repartition'	
	  endif
	 enddo
	endif		!sur nzc
 
c	limites fictives pour faciliter les algorithmes
	
	if(convd(1) .ne. 1)then
	 convf(0)=1
	else
	 convf(0)=10000
	endif
	if(convf(nzc) .ne. nc)then
	 convd(nzc+1)=nc
	else
	 convd(nzc+1)=-10000
	endif
	
c-----------------------------------------------------------------
 
c	if(.true.)then		!test
	if(.false.)then
	print*,'nzc',nzc,nc
	 do i=0,nzc+1
	  write(6,*)i,convd(i),convf(i),mc(convd(i)),mc(convf(i))
	 enddo	!i
	 write(6,*)' '
	 j=1
	 do i=1,nc
	  if(i .eq. convd(j))then
	   if(i .ne. 1)write(6,*)mc(i),i,mc(i)**(3.d0/2.d0)	
	   write(6,*)'debut de ZC',j,i
	   write(6,*)mc(i),i,mc(i)**(3./2.)
	  elseif(i .eq. convf(j))then
	   write(6,*)mc(i),i,mc(i)**(3./2.)
	   write(6,*)'fin de ZC',j,i
	   if(i .ne. nc)write(6,*)mc(i),i,mc(i)**(3.d0/2.d0)	
	   j=min(j+1,nzc)
	  else
	   write(6,*)mc(i),i,mc(i)**(3.d0/2.d0)
	  endif
	 enddo
	 pause'apres ajustement des ZR et des ZC'
	
	endif
	
c-------------------------------------------------------------------------
 
c	dans ZC
c	g(i,nbelem) on met la moyenne des Xi sur la ZC
c	g(i,iw) vitesse angulaire constante de la rotation solide de la ZC
c	dans ZR
c	dans compv(i,nbelem) c'est la comp. chim.
c	dans compv(i,iw) moment angulaire
 
c-----------------------  diffusion --------------------------------------
 
	if(diffusion .and. .not.tot_conv)then !melange classique si tot_conv
 
c	 do i=1,nc
c	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
c	1	knotc,.true.,mc(i),ll,comp1,bidd)
c	  write(6,2002)mc(i),(comp1(j),j=1,nbelem)
2002	  format(1p13d8.1)
c	 enddo
 
	 call diffus_3(bp,q,qt,n,knot,chim,mc,mct,nc,knotc,
     1	r2,m23,mstar,chim_t,mc_t,mct_t,nc_t,knotc_t,
     2	dt,convd,convf,nzc,tot_rad,ok,new,	
     3	old_m23,new_m23,nm,new_m23t,knotm,
     4	etat,opa,conv,nuc,coeff_diff)	
 
c	 estimation de la precision
 
c	 write(6,2002)(mct_t(i),i=1,knotc_t)
c	 pause'estimation de la precision'
 
	 do i=1,nc
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(i),ll,comp1,bidd)
	  bid=min(mc(i),new_m23(nm))
	  call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	bid,ll,f,dfdx)	!masse au temps t mk_t
	  mk_t=min(max(mc_t(1),f(1)),mc_t(nc_t))
	  call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,	
     1	knotc_t,.true.,mk_t,ll,comp2,bidd)
	
c	  estimation de la precision
 
	  do j=1,nbelem
	   esti(j)=abs(comp1(j)-comp2(j))/max(abs(comp1(j)),abs(comp2(j)))
	   if(abs(comp1(j)) .gt. ab_min(j))then
	    if(j .eq. 1 .or. j .eq. ihe4)est=max(est,esti(j))
	    if(est .eq. esti(j))kmax=nc
	    if(estim(j) .le. esti(j))then		!precision par element
c	     write(6,*)j,estim(j),esti(j),est
	     estim(j)=esti(j)
	     kmaxi(j)=i
	     ab_max(j)=comp1(j)
	     mc_maxi(j)=mc(i)
	    endif
	   endif
	  enddo
	 enddo
	
c	 pause'apres diffusion'
 
c----------------------------- diffusion ------------------------------------
 
	else	!integration et melange classiques
 
c	 write(6,*)'knotc',knotc,nbelem
c	 write(6,2000)dt
 
c	 write(6,*)'avant integration : nc_t,knotc_t,new',nc_t,knotc_t,new
c	 write(6,*)'mc_t',mc_t(nc_t)
c	 write(6,2000)(mc_t(i),i=1,nc)
c	 write(6,*)'mct_t',mct_t(knotc_t)
c	 write(6,2000)(mct_t(i),i=1,knotc_t)
c	 print*,'xchim'
c	 do i=1,nc_t
c	  call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,
c	1	knotc_t,.true.,mc_t(i),ll,compx,compy)
c	  write(6,2000)mc_t(i),(compx(j)*nucleo(j),j=1,min(7,nbelem))
c	  write(6,2002)mc_t(i),(compx(j)*nucleo(j),j=1,nbelem)
c	 enddo
c	 pause'avant integration'
 
c	 on integre par rapport a t en melangeant les ZM
c	 d'abord les ZC puis les ZR
c	 on garde dans g(nbelem,nzc) les valeurs des ZM
c	 la variable en masse est mc=m**2/3
 
c	 convf(0)=1 si convd(1) .ne. 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc si convf(nzc) .ne. nc, sinon convd(nzc+1)=-10000
 
c	 pour le moment angulaire Mw, si iw > 1 il y a melange de Mw
c	 dans les ZM; on integre w r**2 au temps t et on divise par
c	 (somme dm au temps t+dt)
c	 dans ZR on conserve Mw = w r**2
 
c	 ensuite on reprend la solution avec ou sans diffusion;
c	 dans ZM on la multipliera Mw (moyen) par
c	 (somme dm temps t+dt) / (somme r**2 dm temps t+dt) et on multipliera
c	 par r**2; on reconstruit le vecteur de composition chimique avec
c	 discontinuite des xchim aux LMR.
c	 ainsi w, en tout point, sera obtenu en divisant Mw interpole, par r**2
 
c	 les ZC d'abord
 
	 do izc=1,nzc
	  mass_zc=0.
	  do i=convd(izc),convf(izc)-1
	   k=i-convd(izc)+1
	   mk=(mc(i)+mc(i+1))/2.d0		!point milieu
	   dm(k)=(mc(i+1)-mc(i))*sqrt(mk)
	   mass_zc=mass_zc+dm(k)		!masse de la ZC
c	   write(6,2000)mass_zc,dm(k),mc(i),mc(i+1),sqrt(mk),mk
 
c	   la composition chimique en mk_t au temps t, et en mk au temps t+dt
 
	   bid=min(mk,new_m23(nm))
	   call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	bid,ll,f,dfdx)	!masse au temps t mk_t
	   mk_t=min(max(mc_t(1),f(1)),mc_t(nc_t))
	   call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,
     1	knotc_t,.true.,mk_t,ll,comp1,bidd)		!perte de Mw
	   if(iw .gt. 1 .and. mk_t .gt. 0.)comp1(iw)=comp1(iw)*mk/mk_t
	
c	   write(6,*)'apres sbsp1dn 1, izc,i,k',izc,i,k
c	   write(6,2000)mk,(comp1(j),j=1,min(7,nbelem))
c	   pause'apres sbsp1dn 1'
 
	   do j=1,nbelem
	    comp1(j)=max(comp1(j),1.d-100)
	    compx(nbelem*(k-1)+j)=comp1(j)
	   enddo
	   if(new)then		!comp. chim. en t+dt
	    do j=1,nbelem
	     comp2(j)=comp1(j)
	     compy(nbelem*(k-1)+j)=comp2(j)
	    enddo
	   else
	    call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(mk,mct(knotc)),ll,comp2,bidd)
c	    write(6,2000)mk,(comp2(j),j=1,min(7,nbelem))
c	    pause'apres sbsp1dn 2'
	    do j=1,nbelem	
	     comp2(j)=abs(comp2(j))
	     compy(nbelem*(k-1)+j)=comp2(j)
	    enddo	!j
	   endif		!new
	   call chim_gram_3(comp1,bidd,nuc_m)
	   call chim_gram_3(comp2,bidd,nuc_m)	
	
c	   les T et ro en t+dt, t, t-dt_, t-dt_-dt__, par interpolation
c	   inverse en m23 par inter_3
 
c	   write(6,*)'avant, izc,i,k',izc,i,k
 
	   call inter_3('mu',bp,q,qt,n,knot,min(mk,m23(n)),f,dfdx,r2,m23)
	   p=exp(f(1))
	   t(k)=exp(f(2))
	   call etat(p,t(k),comp2,.false.,			!les ro
     1	ro(k),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	   write(6,2000)p,t(k),ro(k),mk
	
	   call inter_3('mu',bp_t,q_t,qt_t,n_t,knot_t,min(mk_t,m23_t(n_t)),
     1	f,dfdx,r2_t,m23_t)
	   		
	   p_t=exp(f(1))
	   t_t(k)=exp(f(2))
	   call etat(p_t,t_t(k),comp1,.false.,
     1	ro_t(k),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	   write(6,2000)p_t,t_t(k),ro_t(k),mk
 
	   if(ordre .gt. 2)then		!sans perte de masse necessairement
	    if(dt_ .gt. 0.d0)then
	
	     call sbsp1dn(1,old_m23_,new_m23_,new_m23t_,nm_,2,knotm_,.true.,
     1	min(mk_t,new_m23_(nm_)),ll,f,dfdx)		
	     mk_t=f(1)
	     call inter_3('mu',bp_,q_,qt_,n_,knot_,min(mk_t,m23_(n_)),
     1	f,dfdx,r2_,m23_)
	     p_=exp(f(1))
	     t_(k)=exp(f(2))
	     call sbsp1dn(nbelem,chim_,mc_,mct_,nc_,m_ch,
     1	knotc_,.true.,min(mk,mc_(nc_)),ll,comp1,bidd)
	     call chim_gram_3(comp1,bidd,nuc_m)		
	     call etat(p_,t_(k),comp1,.false.,
     1	ro_(k),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	     write(6,2000)p_,t_(k),ro_(k)
	    endif
	    if(dt__ .gt. 0.d0)then
	     call sbsp1dn(1,old_m23__,new_m23__,new_m23t__,nm__,2,knotm__,
     1	.true.,min(mk_t,new_m23__(nm__)),ll,f,dfdx)		
	     mk_t=f(1)
	     call inter_3('mu',bp__,q__,qt__,n__,knot__,min(mk_t,m23__(n__)),
     1	f,dfdx,r2_,m23_)
	     p__=exp(f(1))
	     t__(k)=exp(f(2))
	     call sbsp1dn(nbelem,chim__,mc__,mct__,nc__,m_ch,
     1	knotc__,.true.,min(mk,mc__(nc__)),ll,comp1,bidd)
	     call chim_gram_3(comp1,bidd,nuc_m)		
	     call etat(p__,t__(k),comp1,.false.,
     1	ro__(k),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	     write(6,2000)p__,t__(k),ro__(k)
	    endif
	   endif
	  enddo		!i dans chaque ZC
 
c	  print*,k
c	  write(6,2000)(compx(nbelem*(i-1)+1),i=1,k)
c	  write(6,2000)(compy(nbelem*(i-1)+1),i=1,k)	  	
c	  write(6,2000)(dm(i),i=1,k)
c	  write(6,2000)(t(i),i=1,k)
c	  write(6,2000)(ro(i),i=1,k)	  	
c	  write(6,2000)mass_zc
c	  pause'apres les dm(k)'
 
c	  integration temporelle
 
	  do i=1,k
	   dm(i)=dm(i)/mass_zc
	  enddo
 
	  call rk_imps_3(t__,ro__,dt__,t_,ro_,dt_,
     1	t_t,ro_t,compx,t,ro,compy,dt,esti,ok,nuc,k,dm)
	
c	  print*,ok
c	  write(6,2000)(compx(i),i=1,nbelem)	
c	  write(6,2000)(compy(i),i=1,nbelem)
c	  pause'apres rk_imps_3'
	
c	  on garde, provisoirement, les abondances dans g
	
	  if(ok)then		!il y a eu convergence
	   do j=1,nbelem		
	    g(nbelem*(izc-1)+j)=compy(j) !g: comp. chim. de la ZC
	
	    if(estim(j) .le. esti(j))then
c	     write(6,*)j,estim(j),esti(j),est
	     estim(j)=esti(j)!precision par element
	     kmaxi(j)=convd(izc)
	     ab_max(j)=compy(j)
	     mc_maxi(j)=mc(convd(izc))
	    endif
	   enddo		!j=1,nbelem
	   bid=abs(compy(1)-compx(1))
	   est=max(est,bid)
	   if(est .eq. bid)kmax=convd(izc)
	   if(ihe4 .gt. 1)then
	    bid=abs(compy(ihe4)-compx(ihe4))
	    est=max(est,bid)
	    if(est .eq. bid)kmax=convd(izc)
	   endif
	  else			!ok
	   write(6,*)'probleme a la ZC',izc
	   g_max=.true.			!reinitialisation et dt-->dt/2
	   if(dt .gt. dtmin)return	!pas de CV dans rk_imps_3
	   write(6,*)'EVOL_3(2) dt trop petit abandon/dt,dtmin'
	   write(6,2000)dt,dtmin
	   !pause'ABANDON'
	   stop
	  endif
c	  print*,'zone melangee',izc
c	  write(6,2000)mc(convd(izc)),mc(convf(izc))
c	  write(6,2002)(g(nbelem*(izc-1)+j),j=1,nbelem)
c	  pause'ZM'
	 enddo	!izc
	
c 	 ensuite les ZR
 
	 dm(1)=1.d0 !fictif
	 do izc=1,nzc+1
c	  print*,'ZR',izc,convf(izc-1),convd(izc)
	  do i=convf(izc-1),convd(izc)
 
c	   la composition chimique en mk_t au temps t et en mc(i) au temps t+dt
c	   les mc peuvent differer des mct
 
	   bid=min(mc(i),new_m23(nm))
	   call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	bid,ll,f,dfdx)	!masse au temps t mk_t
	   mk_t=min(max(mc_t(1),f(1)),mc_t(nc_t))
	   call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,	
     1	knotc_t,.true.,mk_t,ll,comp1,bidd)	!perte de Mw
	   if(iw .gt. 1 .and. mc(i) .gt. 0.)comp1(iw)=comp1(iw)*mc(i)/mk_t
	
c	   write(6,*)'apres sbsp1dn 1, izc,i',izc,i
c	   write(6,2000)mc(i),(comp1(j),j=1,min(7,nbelem))
c	   pause'apres sbsp1dn 1 radiatif'
 
	   do j=1,nbelem
	    comp1(j)=max(comp1(j),1.d-100)
	    compx(j)=comp1(j)	
	   enddo
	   if(new)then		!comp. chim. en t+dt
	    do j=1,nbelem
	     comp2(j)=comp1(j)
	     compy(j)=comp2(j)
	    enddo
	   else
	    call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(mc(i),mct(knotc)),ll,comp2,bidd)
c	    write(6,2000)(comp2(j),j=1,min(7,nbelem))
c	    pause'apres sbsp1dn 2 radiatif'
	    do j=1,nbelem	
	     comp2(j)=abs(comp2(j))
	     compy(j)=comp2(j)
	    enddo	!j
	   endif		!new
	   call chim_gram_3(comp1,bidd,nuc_m)		
	   call chim_gram_3(comp2,bidd,nuc_m)		
	
c	   les T et ro en t+dt, t, t-dt_, t-dt_-dt__, par interpolation
c	   inverse en m23 par inter_3
 
c	   write(6,*)'avant, izc,i',izc,i
 
	   call inter_3('mu',bp,q,qt,n,knot,min(mc(i),m23(n)),f,dfdx,r2,m23)
	   p=exp(f(1))
	   t(1)=exp(f(2))
	   call etat(p,t(1),comp2,.false.,			!les ro
     1	ro(1),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	   write(6,2000)p,t(1),ro(1),mc(i)
 
	   call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	min(mc(i),new_m23(nm)),ll,f,dfdx)
	   mk_t=f(1)	
	   call inter_3('mu',bp_t,q_t,qt_t,n_t,knot_t,min(mk_t,m23_t(n_t)),
     1	f,dfdx,r2_t,m23_t)
	   p_t=exp(f(1))
	   t_t(1)=exp(f(2))
	   call etat(p_t,t_t(1),comp1,.false.,
     1	ro_t(1),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	   write(6,2000)p_t,t_t(1),ro_t(1),mc(i)
 
	   if(ordre .gt. 2)then
	    if(dt_ .gt. 0.d0)then	!sans perte de masse necessairement
	     call sbsp1dn(1,old_m23_,new_m23_,new_m23t_,nm_,2,knotm_,.true.,
     1	min(mk_t,new_m23_(nm_)),ll,f,dfdx)		
	     mk_t=f(1)
	     call inter_3('mu',bp_,q_,qt_,n_,knot_,min(mk_t,m23_(n_)),
     1	f,dfdx,r2_,m23_)
	     p_=exp(f(1))
	     t_(1)=exp(f(2))
	     call sbsp1dn(nbelem,chim_,mc_,mct_,nc_,m_ch,
     1	knotc_,.true.,min(mc(i),mc_(nc_)),ll,comp1,bidd)
	     call chim_gram_3(comp1,bidd,nuc_m)		
	     call etat(p_,t_(1),comp1,.false.,
     1	ro_(1),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	     write(6,2000)p_,t_(1),ro_(1)
	    endif
	    if(dt__ .gt. 0.d0)then
	     call sbsp1dn(1,old_m23__,new_m23__,new_m23t__,nm__,2,knotm__,
     1	.true.,min(mk_t,new_m23__(nm__)),ll,f,dfdx)		
	     mk_t=f(1)
	     call inter_3('mu',bp__,q__,qt__,n__,knot__,min(mk_t,m23__(n__)),
     1	f,dfdx,r2_,m23_)
	     p__=exp(f(1))
	     t__(1)=exp(f(2))
	     call sbsp1dn(nbelem,chim__,mc__,mct__,nc__,m_ch,
     1	knotc__,.true.,min(mc(i),mc__(nc__)),ll,comp1,bidd)
	     call chim_gram_3(comp1,bidd,nuc_m)	
	     call etat(p__,t__(1),comp1,.false.,
     1	ro__(1),drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c	     write(6,2000)p__,t__(1),ro__(1)
	    endif
	   endif
 
c	   integration temporelle
 
	   call rk_imps_3(t__,ro__,dt__,t_,ro_,dt_,
     1	t_t,ro_t,compx,t,ro,compy,dt,esti,ok,nuc,1,dm)
c	   print*,'dans ZR',i
c	   write(6,2000)mc(i),t(1),ro(1),(compy(j),j=1,min(4,nbelem))
	
c	   on place, provisoirement les abondances dans compv
	
	   if(ok)then		!il y a eu convergence
	    do j=1,nbelem
	     compv(nbelem*(i-1)+j)=compy(j)
	     if(abs(compy(j)) .gt. ab_min(j))then	
	      if(estim(j) .le. esti(j))then
c	       write(6,*)j,estim(j),esti(j),est
	       estim(j)=esti(j)	!precision par element
	       kmaxi(j)=i
	       ab_max(j)=compy(j)
	       mc_maxi(j)=mc(i)
	      endif
	     endif
	    enddo		!j=1,nbelem
	    bid=abs(compy(1)-compx(1))
	    est=max(est,bid)
	    if(est .eq. bid)kmax=i
	    if(ihe4 .gt. 1)then
	     bid=abs(compy(ihe4)-compx(ihe4))
	     est=max(est,bid)
	     if(est .eq. bid)kmax=i
	    endif
	
	   else			!ok
	    write(6,*)'probleme a la couche',i
	    g_max=.true.			!reinitialisation et dt-->dt/2
	    if(dt .gt. dtmin)return	!pas de CV dans rk_imps_3
	    write(6,*)'EVOL_3(2) dt trop petit abandon/dt,dtmin'
	    write(6,2000)dt,dtmin
	    !pause'ABANDON'
	    stop
	   endif
	  enddo	!i
	 enddo		!izc
	
c	 cas sans diffusion, mouvements des limites ZR/ZC	
 
	 if(.not.tot_conv .or. .not.tot_rad)then
	  tot_conv_t=lim_t .eq. 1 .and. jlim_t(1) .eq. n_t
     1		.and. .not.lconv_t(1)
	
c	  cas de la regression d'un coeur convectif :
c	  on interpole lineairement la comp.chim. sur la zone de regression
 
	  core_reg=.not.lconv(1) .and. .not.lconv_t(1)
     1	.and. m_zc(1) .lt. m_zc_t(1) .and. m_zc(1) .lt. mstar
	
c	  core_reg=.false.
	
c	  print*,core_reg,lconv(1),lconv_t(1),m_zc(1),m_zc_t(1)
	  if(core_reg)then
	   do j=1,nbelem !droite de lim. on met la comp. chim. de la ZC
	    compv(nbelem*(convf(1)-1)+j)=g(nbelem*(1-1)+j)
	   enddo	
	   m_zc23t=m_zc_t(1)**(2.d0/3.d0)
	   call slinf(m_zc23t,mc,nc,ll)
c	   print*,convf(1),ll,mc(convf(1)),m_zc23t
	
c	   la ZC a regresse de m_zc23t a mc(convf(1)) et on a :
c	   mc(ll+1) .gt. m_zc23t .ge. mc(ll) .ge. mc(convf(1))
c	   interp. lin. de la comp. chim. pour mc entre
c	   mc(ll+1) et mc(convf(1))
	
	   do i=convf(1)+1,ll
c	    print*,i,mc(i)
	    bid=(mc(i)-mc(convf(1)))/(mc(ll+1)-mc(convf(1)))
	    do j=1,nbelem
	     compv(nbelem*(i-1)+j)=compv(nbelem*(convf(1)-1)+j)+
     1	(compv(nbelem*ll+j)-compv(nbelem*(convf(1)-1)+j))*bid
	    enddo
	   enddo	
	  endif
	
c	  cas de l'augmentation de l'abscisse lagrangienne (masse) de la
c	  base d'une ZC suffisamment interne
c	  on interpole lineairement la comp.chim. sur la zone de regression
 
	  do izc=1,nzc
	   if(mc(convd(izc)) .le. mc_max)then
	    if(tot_conv_t .and. convd(izc) .gt. 1)then
	     zc_aug=.true.
	     m_zc23t=0
	    elseif(lim_t .eq. 1)then		!au temps t : une seule ZC
	     m_zc23t=m_zc_t(1)**(2.d0/3.d0)
	     zc_aug=m_zc23t .lt. mc(convd(izc)) .and. lconv_t(1)
	    else		!plusieurs ZC au temps t
	     jfin=lim_t-1
	     do while(it .le. jfin)
	      m_zc23t=m_zc_t(it)**(2.d0/3.d0)
	      zc_aug=m_zc23t .lt. mc(convd(izc)) .and. lconv_t(it)
     1	 .and. m_zc_t(it+1)**(2.d0/3.d0) .gt. mc(convd(izc))
	      if(zc_aug)then
	       jfin=-4	!il y a eu encadrement
	      else
	       it=it+1
	      endif
	     enddo
	    endif
	
c	    print*,zc_aug,izc,convd(izc),it,tot_conv_t,lim_t,
c	1	jlim_t(1),n_t,lconv_t(1)
c	    write(6,2000)m_zc23t,mc(convd(izc)),m_zc_t(it+1)**(2.d0/3.d0)
 
c	    zc_aug=.false.		!suppression
 
	    if(zc_aug)then
	     call slinf(m_zc23t,mc,nc,ll)
c	     print*,ll
	 	
c	     la limite de la ZC est passee de m_zc23t=mc(ll) a mc(convd(iz))
c	     interp. lin. de la comp. chim. pour mc entre mc(ll)
c	     et mc(convd(izc))
c	     pour iw on interpole lineairement le moment angulaire
 
	     do i=ll+1,convd(izc)
c	      print*,i,mc(i),mc(ll)
	      bid=(mc(i)-mc(convd(izc)))/(mc(ll)-mc(convd(izc)))
	      do j=1,nbelem
	       compv(nbelem*(i-1)+j)=g(nbelem*(izc-1)+j)+
     1	(compv(nbelem*(ll-1)+j)-g(nbelem*(izc-1)+j))*bid
	      enddo
	     enddo
c	     print*,convd(izc)
c	     write(6,2000)(compv(nbelem*(j-1)+iw),j=1,convd(izc)-1),
c	1	compv(nbelem*(ll-1)+iw)
	    endif
	   endif
	  enddo
	 endif
c	 pause'apres mouvements des ZR/ZC'
 
c	 on modifie le vecteur de composition chimique
 
c	 convf(0)=1 si convd(1) .ne. 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc si convf(nzc) .ne. nc, sinon convd(nzc+1)=-10000
 
c	 print*,nzc
c	 print*,(convf(i),i=0,nzc)
c	 print*,(convd(i),i=1,nzc+1)
c	 pause'nzc'
	 	
	 ndis=0		!nombre de discontinuites
	 do izc=1,nzc+1
	  do i=convf(izc-1),convd(izc)	!zone radiative  ZR
	   if(i .eq. convf(izc-1) .and. i .gt. 1)then
	    ndis=ndis+1
	    idis(ndis)=i
	    do j=1,nbelem
	     chimd(nbelem*(ndis-1)+j)=compv(nbelem*(i-1)+j)
	    enddo	!j
c	    print*,'ndis,i',ndis,i
c	    write(6,2000)compv(nbelem*(i-1)+1),chimd(nbelem*(ndis-1)+1)
c	    write(6,2000)chimd(nbelem*(ndis-1)+iw)
c	    print*,mc(i)
c	    pause'debut de ZR a droite de ZM'
	
	   else
	    do j=1,nbelem
	     chim(nbelem*(i-1)+j)=compv(nbelem*(i-1)+j)
	    enddo
	   endif
	  enddo
	
	  if(izc .le. nzc)then		!pour ZM	
	   do i=convd(izc),convf(izc)
	    if(i .eq. convd(izc) .and. i .gt. 1)then
	     ndis=ndis+1
	     idis(ndis)=i
	     do j=1,nbelem	
	      chimd(nbelem*(ndis-1)+j)=g(nbelem*(izc-1)+j)
	     enddo
c	     print*,'ndis,i,izc',ndis,i,izc
c	     write(6,2000)g(nbelem*(izc-1)+1),chimd(nbelem*(ndis-1)+1)
c	     pause'debut de ZM a droite de ZR'
 
	    else
	     do j=1,nbelem	
	      chim(nbelem*(i-1)+j)=g(nbelem*(izc-1)+j)
	     enddo
	    endif	!i .eq. convd
	   enddo	!i=convd , convf
	  endif	!izc le. nzc
	 enddo	!nzc
	
c	 print*,nc,nzc,ndis,(idis(i),i=1,ndis)
c	 ll=1
c	 do i=1,nc
c	  write(6,2000)mc(i),(chim(nbelem*(i-1)+j),j=1,min(4,nbelem))
c	  write(6,2000)mc(i),chim(nbelem*(i-1)+iw)	
c	  if(i .eq. idis(ll))then
c	   print*,idis(ll)	
c	   write(6,2000)mc(i),(chimd(nbelem*(ll-1)+j),j=1,min(4,nbelem))
c	   write(6,2000)mc(i),chimd(nbelem*(ll-1)+iw)
c	   print*,mc(i)	
c	   print*
c	   ll=min(ll+1,ndis)
c	  endif
c	 enddo
c	 do i=1,ndis
c	  write(6,2000)(chimd(nbelem*(i-1)+j),j=1,min(4,nbelem))
c	  write(6,2000)chimd(nbelem*(i-1)+iw)	
c	 enddo
c	 print*,ndis,(idis(i),i=1,ndis)	
c	 write(6,2000)(chim(nbelem*(i-1)+iw),i=1,nc)
c	 write(6,2000)(chimd(nbelem*(ll-1)+iw),ll=1,ndis)
c	 pause'avant sbsp_dis'
 
	 if(tot_rad .or. tot_conv)then
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.false.,mc(1),ll,comp2,bidd)	
	 else
	  call sbsp_dis(nbelem,mc,chim,ndis,idis,chimd,nc,m_ch,mct,knotc)
	 endif
 
c	 write(6,*)'apres integration : nc,knotc',nc,knotc
c	 write(6,*)'mc',mc(nc)
c	 write(6,2000)(mc(i),i=1,nc)
c	 write(6,*)'mct',mct(knotc)
c	 write(6,2000)(mct(i),i=1,knotc)
c	 print*,'xchim'
c	 do i=1,nc
c	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
c	1	knotc,.true.,mc(i),ll,compx,compy)
c	  write(6,2000)mc(i),(compx(j)*nucleo(j),j=1,min(7,nbelem))
c	  write(6,2002)mc(i),(compx(j)*nucleo(j),j=1,nbelem)
c	  write(6,2002)mc(i),(compx(j),j=1,nbelem)	
c	 enddo
c	 pause'apres integration'
c
c	 pause'test rota'
c	 do i=1,nc
c	  call inter_3('mu',bp,q,qt,n,knot,min(m23(i),m23(n)),f,dfdx,r2,m23)
c	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
c	1	knotc,.true.,m23(i),ll,compx,compy)
c	  if(r2(i) .gt. 0.)print*,m23(i),r2(i),compx(iw)/r2(i),compx(1)	
c	  if(f(3) .gt. 0.)then
c	   print*,m23(i),f(3),compx(iw)/f(3)
c	  else
c	   print*,m23(i),f(3),compy(iw)/dfdx(3)*dfdx(5)
c	  endif
c	  write(6,2000)compx(iw)
c	 enddo
c	 pause'apres test'
 
	endif	!diffusion ou integration et melange classiques
 
c-----------debut de rotation non solide------------------------
 
c	if(.false.)then
	if(iw .gt. 1 .and. .not.tot_rad)then
 
c	 avec rotation non solide (iw > 1) on determine la valeur du moment
c	 angulaire dans les ZC
 
c	 calcul de somme r2 dm et de somme dm dans chaque ZC
	
	 do i=1,nzc
	  mass_zc=0
	  bid=0
	  do j=convd(i),convf(i)-1
	   mk=(mc(j)+mc(j+1))/2.		!point milieu
	   dm(j)=(mc(j+1)-mc(j))*sqrt(mk)
	   mass_zc=mass_zc+dm(j)	
	   call inter_3('mu',bp,q,qt,n,knot,mk,f,dfdx,r2,m23)
	   bid=bid+dm(j)*f(3)
	  enddo
c	  write(6,2000)mass_zc,bid
	  dmsr2dm(i)=mass_zc/bid
	 enddo			!f,dfdx: vt
c	 pause'les integrales'
	
c	 convf(0)=1 si convd(1) .ne. 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc si convf(nzc) .ne. nc, sinon convd(nzc+1)=-10000
 
c	 print*,nzc
c	 print*,(convf(i),i=0,nzc)
c	 print*,(convd(i),i=1,nzc+1)
c	 pause'nzc'
	 	
	 ndis=0		!nombre de discontinuites
	 do izc=1,nzc+1
	  do i=convf(izc-1),convd(izc)	!zone radiative  ZR
	   if(i .eq. convf(izc-1) .and. i .gt. 1)then
	    call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(i)+1.d-14,ll,compx,compy)	
	    ndis=ndis+1
	    idis(ndis)=i
	    do j=1,nbelem
	     chimd(nbelem*(ndis-1)+j)=compx(j)
	    enddo	!j
c	    print*,'debut de ZR a droite de ZC ndis,i',ndis,i,mc(i)+1.d-14	
c	    write(6,2000)chimd(nbelem*(ndis-1)+iw),compx(iw)
 
	   elseif(i .eq. convd(izc) .and. i .lt. nc)then
	    call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(i)-1.d-14,ll,compx,compy)  !petit decalage
	    do j=1,nbelem
	     compv(nbelem*(i-1)+j)=compx(j)
	    enddo	!j
c	    print*,'fin de ZR,i, mc(i)',i,mc(i)-1.d-14	
c	    write(6,2000)compv(nbelem*(i-1)+iw)
c	    print*	
 
	   else	
	    call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(i),ll,compx,compy)
	    do j=1,nbelem
	     compv(nbelem*(i-1)+j)=compx(j)
	    enddo
c	    print*,'dans ZR,i',i
c	    write(6,2000)compv(nbelem*(i-1)+iw)
 
	   endif
	  enddo		!pour la ZR
	
	  if(izc .le. nzc)then		!pour ZM	
	   do i=convd(izc),convf(izc)
	    if(i .eq. convd(izc) .and. i .gt. 1)then
	     call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(i)+1.d-14,ll,compx,compy)
	     ndis=ndis+1
	     idis(ndis)=i
	     do j=1,nbelem		
	      if(j .eq. iw)then
	       call inter_3('mu',bp,q,qt,n,knot,min(mc(i),m23(n)),
     1	f,dfdx,r2,m23)
	       chimd(nbelem*(ndis-1)+iw)=compx(iw)*dmsr2dm(izc)*f(3)	
	      else
	       chimd(nbelem*(ndis-1)+j)=compx(j)
	      endif
	     enddo
c	     print*,'debut de ZC a droite de ZR, ndis,i,izc',ndis,i,izc,f(3)	
c	     write(6,2000)chimd(nbelem*(ndis-1)+iw),compx(iw),dmsr2dm(izc)
c	     print*,'debut de ZC a droite de ZR'
c	     write(6,2000)compx(iw)*dmsr2dm(izc),compx(iw)/w_rot,1./dmsr2dm(izc)
	
	    elseif(i .eq. convf(izc) .and. i .lt. nc)then
	     call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(i)-1.d-14,ll,compx,compy)
	     do j=1,nbelem	
	      if(j .eq. iw)then
	       call inter_3('mu',bp,q,qt,n,knot,min(mc(i),m23(n)),
     1	f,dfdx,r2,m23)
	       compv(nbelem*(i-1)+iw)=compx(iw)*dmsr2dm(izc)*f(3)	
	      else
	       compv(nbelem*(i-1)+j)=compx(j)
	      endif
	     enddo
c	     print*,'fin de ZC,i',i,f(3)	
c	     write(6,2000)compv(nbelem*(i-1)+iw),compx(iw),dmsr2dm(izc)
c	     print*,'fin de ZC'
c	     write(6,2000)compx(iw)*dmsr2dm(izc),compx(iw)/w_rot,1./dmsr2dm(izc)
c	     print*
	
	    else
	     call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(i),ll,compx,compy)
	     do j=1,nbelem	
	      if(j .eq. iw)then
	       call inter_3('mu',bp,q,qt,n,knot,min(mc(i),m23(n)),
     1	f,dfdx,r2,m23)
	       compv(nbelem*(i-1)+iw)=compx(iw)*dmsr2dm(izc)*f(3)	
	      else
	       compv(nbelem*(i-1)+j)=compx(j)
	      endif	!iw
	     enddo	!j
c	     print*,'dans ZC,i',i,f(3)
c	     write(6,2000)compv(nbelem*(i-1)+iw),compx(iw),dmsr2dm(izc)
c	     print*,'dans ZC'
c	     write(6,2000)compx(iw)*dmsr2dm(izc),compx(iw)/w_rot,1./dmsr2dm(izc)
 
	    endif	!i .eq. convd
	
	   enddo	!i=convd , convf
	  endif	!izc le. nzc
	 enddo	!nzc
	
	 do i=1,pchim		!compv ---> chim
	  chim(i)=compv(i)
	 enddo
	
c	 print*,nc,ndis,(idis(i),i=1,ndis)
c	 write(6,2000)(chim(nbelem*(i-1)+iw),i=1,nc)
c	 write(6,2000)(chimd(nbelem*(ll-1)+iw),ll=1,ndis)	
c	 pause'fin de rota'
 
	 if(ndis .le. 0)then
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.false.,mc(1),ll,comp2,bidd)
	 else
	  call sbsp_dis(nbelem,mc,chim,ndis,idis,chimd,nc,m_ch,mct,knotc)
	 endif
	endif		!iw > 1
	
c--------------fin de rotation non solide------------------------
	
c	estimation du pas temporel suivant
 
	if(est .ne. 0.d0)then
	 dtnew=.9*dt*(precit/est)**(1.d0/(ordre+1))
c	 dtnew=.95*dt*precit/est
	else
	 dtnew=1.2*dt
	endif
	dtn=max(.8d0*dt,min(dtnew,1.2d0*dt))
 
	do i=1,nbelem		!les masses avec erreur maximale
	 mc_maxi(i)=sqrt(mc_maxi(i))**3
	enddo
	
c	ecritures	
 
	write(6,111)kmax,est,dt,dtn
	write(6,112)(ab_max(i)*nucleo(i),i=1,nchim)
	write(6,112)(estim(i),i=1,nchim)
	write(6,112)(mc_maxi(i),i=1,nchim)
	write(6,113)(kmaxi(i),i=1,nchim)
111	format(1x,/,1x,'EVOL_3 couche',i4,' variation maximale pour H ou He4:'
     1	,1pd10.3,/,1x,'dt=',1pd10.3,' dt optimal=',1pd10.3,
     2	' abondance // variation/scale // masse // couche:')
112	format(1x,1p12d8.1)
113	format(1x,12(2x,i4,2x))
 
	dts=min(dtmax,dtn,agemax-age-dt)
	g_max=.false.
 
c	write(6,*)'fin evol dt,dts,dtn,dtmax,agemax,age,agemax-age-dt'
c	write(6,2000)dt,dts,dtn,dtmax,agemax,age,agemax-age-dt
 
	return
 
	end
