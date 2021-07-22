c*************************************************************************
 
	subroutine dess_hr_DC(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,teff,
     1	r2,m23,age,lim,lconv,m_zc,mstar)
	
c	dessin de variables en cours d'evolution
c	Auteur: P. Morel,  D. Cordier
c       Adaptation du HR en fonction de la masse de l'étoile : Nov. 2000
c	CESAM_3
 
c entrees, significations evidentes
c	p,t,m,l,r,ro,grad_ad,grad_mj,alfa,delta,kap,cp,teff,age
c	epsilon(5,i): energie thermo et gravifique, voir cesam_3
c	n: nombre de couches
c	chim,mc,mct,nc,knotc: pour interpolation de la comp. chim.
c	lim: nombre de limites ZR/ZC
c	lconv: .true. si debut de ZC
c	m_zc: masse aux limites ZC/ZR
c	mstar: masse totale au temps du dessin
	
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'ctephy.common'
	
	integer pndes
	parameter(pndes=1000)
 
	integer ih,ihe3,ili7,ic12,in14,i,ndes,nc,knotc,lx,
     1	lim,nxsub,nysub,long,n,knot

	external long

	real*8	chim(1),mc(1),mct(1),xdbl,teff,age,dc,dq,f(pne),dfdx(pne),
     1	xchim(pnelem),dxchim(pnelem),m_zc(1),qx(pndes),
     2	bp(1),q(1),qt(1),nuc_m(pnelem),mstar,r2(1),m23(1)
	
	real*4	x(pndes),y1(pndes),y2(pndes),y3(pndes),y4(pndes),y5(pndes),
     1	xmin,xmax,ymin,ymax,inch,xleft,ybot,xmin1,xmax1,ymin1,ymax1,
     2	xtick,ytick,disp,coord,fjust,h,dh,ld,dl,y6(pndes),y7(pndes),pas,
     3	xmax3,xmin3,ymax3,ymin3,xmin4,xmax4,lteff,llext,lteffp,llextp,
     4	ymin4,ymax4,x4(2),qd(pndes),pmax,tmax,rmax,lmax,mmax,yp(pndes),
     5	yt(pndes),yr(pndes),yl(pndes),ym(pndes),ym2,ym4,ym5,ym6,ym7
c	data inch/2.54/,pas/-1./	
 
	character*10 device,xopt,yopt,side
	character*50 text
 
	logical init,lconv(1),efface
 
	data init/.true./,efface/.false./	
	data inch/2.54/,pas/-1./	
 
2000	format(1x,1p8d10.3)	
 
	if(init)then
	 init=.false.
	
c	 indice de He4
 
	 ih=1
	 ihe3=0
	 ili7=0
	 ic12=0
	 in14=0
	 do i=1,nchim
	  if(nom_elem(i) .eq. 'He3')ihe3=i
	  if(nom_elem(i) .eq. 'Li7')ili7=i
	  if(nom_elem(i) .eq. 'C12')ic12=i
	  if(nom_elem(i) .eq. 'N14')in14=i
	 enddo
c	 print*,ih,ihe3,ihe4,ili7,ic12,nchim
c	 pause
 
	 device=' '
	 print*,'entrer le nom du device: /XSERVE, /PS ou ?'
	 print*,'(RETURN = /XSERVE)'
	 read(5,'(a)')device
	 if ( device .eq. ' ' ) then
	    device='/XSERVE'
	 end if
c	 print*,device(:long(device))
c	 pause
	 if(device(:long(device)) .eq. '/PS')then
	  device=nom_fich2(:long(nom_fich2))//'.ps/PS'
	  print*,'nom du fichier du plot: ',nom_fich2(:long(nom_fich2))//'.ps'
	 endif
	 call pgbegin(0,device,1,1)
	 call pgscf(2)		!roman font
c	 call pgsch(1.2)	!hauteur des caracteres
	 call pgslw(3)		!epaisseur du trait	
		
	 ndes=pndes
	 dc=1./ndes			!pour la chimie
	 do i=1,ndes
	  x(i)=dc*i
	 enddo
 
	 h=7./inch	!hauteur du cadre
	 dh=1.5/inch	!saut entre cadres
 
	 ld=10./inch	!largeur du cadre
	 dl=2./inch	!espace en largeur
 
	 xleft=1.8/inch
	 ybot=1.8/inch
 
c	 premier cadre: H, He3, He4
 
	 xmin=-0.01
	 xmax=1.01
	 ymax=1.1
	 ymin=-.01
	
	 call pgvsize(xleft,xleft+ld,ybot,ybot+h)
	 call pgwindow(xmin,xmax,ymin,ymax)
 
	 xopt='bcnst'
	 xtick=0.
	 nxsub=0
 
	 yopt='bcnst'
	 ytick=0.
 
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
 
	 side='b'	!abscisses
 
	 disp=2
	 coord=.5
 
	 fjust=.5
	 text='M/M\\d\\(0856)\\u'
 
	 call pgmtext(side,disp,coord,fjust,text)
	 side='l'	!ordonnees
	 disp=2
 
	 coord=.5
	 fjust=.5
         call pgsci(2)
	 text='X, '
	 call pgmtext(side,disp,coord,fjust,text)
         call pgsci(1)
         text='   Y, '
         call pgmtext(side,disp,coord,fjust,text)
         call pgsci(7)
         text='      Z'
         call pgmtext(side,disp,coord,fjust,text)
         call pgsci(1)
 
c	 second cadre: P, T, R, L en fonction de m/Mstar
	 xmin1=xmin
	 xmax1=xmax
         ymin1=ymin
	 ymax1=ymax
 
	 call pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
	 call pgwindow(xmin1,xmax1,ymin1,ymax1)
	
	 xopt='bcnst'
 
	 xtick=0.
	 nxsub=0
 
	 yopt='bcnst'
	 ytick=0.
 
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
 
	 side='b'	!abscisses
	 disp=2
 
	 coord=.5
	 fjust=.5
 
	 text='M/M\\(0856)\\u'
 
	 call pgmtext(side,disp,coord,fjust,text)
	
	 side='l'	!ordonnees
 
	 disp=2
	 coord=.5
 
	 fjust=.5
	 text='P, '
         call pgsci(2)
	 call pgmtext(side,disp,coord,fjust,text)
         text='    T, '
         call pgsci(3)
         call pgmtext(side,disp,coord,fjust,text)
         call pgsci(5)
         text='       R,'
         call pgmtext(side,disp,coord,fjust,text)
         call pgsci(6)
         text='          L'
         call pgmtext(side,disp,coord,fjust,text)
 
         call pgsci(1)
 
c	 pause
	
	 side='t'	!sommet
	 disp=1
 
	 coord=.5
	 fjust=.5
 
	 call pgmtext(side,disp,coord,fjust,nom_fich2)
			
c	 troisieme cadre: diagramme HR
 
	 lteffp=log10(teff)
	 lx=n
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(n),lx,f,dfdx)	
	 llextp=log10(f(4))*1.5
c        adaptation en fonction de la masse totale ---------------
         if ( mtot .le. 3.d0 ) then
	  xmin3=4.5
          xmax3=3.55
          ymin3=1.0
          ymax3=3.2
         end if
         if ( (mtot .gt. 3.d0) .AND. (mtot .lt. 7.d0) ) then
	  xmin3=4.5
          xmax3=3.55
          ymin3=2.0
          ymax3=4.2
         end if
         if ( mtot .ge. 7.d0 ) then
	  xmin3=4.5
          xmax3=3.55
          ymin3=2.8
          ymax3=5.2
         end if
c        --------------------------------------------------------

	 call pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot,ybot+h)
	 call pgwindow(xmin3,xmax3,ymin3,ymax3)
	
	 xopt='bcnst'
 
	 xtick=0.
	 nxsub=0
 
	 yopt='bcnst'
	 ytick=0.
 
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
 
	 side='b'	!abscisses
	 disp=2
 
	 coord=.5
	 fjust=.5
 
	 text='Log\\d10\\u T\\deff\\u'
	 call pgmtext(side,disp,coord,fjust,text)
 
	 side='l'	!ordonnees
	 disp=2
 
	 coord=.5
	 fjust=.5
 
	 text='Log\\d10\\u L/L\\d\\(2281)\\u'
	 call pgmtext(side,disp,coord,fjust,text)
 
c	 quatrieme cadre: ZC
 
	 xmin4=-3.
	 xmax4=600
	 ymin4=-.05
	 ymax4=1.05
	 call pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot+h+dh,ybot+h+dh+h)
	 call pgwindow(xmin4,xmax4,ymin4,ymax4)
	
	 xopt='bcnst'
 
	 xtick=0.
	 nxsub=0
 
	 yopt='bcnst'
	 ytick=0.
 
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	
	 side='b'	!abscisses
 
	 disp=2
	 coord=.5
	 fjust=.5
 
	 text='pas temporel'
	 call pgmtext(side,disp,coord,fjust,text)
	
	 side='l'	!ordonnees
 
	 disp=2
	 coord=.5
 
	 fjust=.5
	 text='ZC ---- M/M\\d\\(0856)\\u'
 
	 call pgmtext(side,disp,coord,fjust,text)	
c	 pause
 
	 side='t'	!sommet
	 disp=1
 
	 coord=.5
	 fjust=.5
 
	 call pgmtext(side,disp,coord,fjust,nom_fich2)
	
	 call pgslw(1)		!epaisseur du trait
	 	
	endif
c-------------------------------------------------------------
c	les abondances
	if(efface)then
	 call pgvsize(xleft,xleft+ld,ybot,ybot+h)	
	 call pgwindow(xmin,xmax,ymin,ymax)
	 call pgsls(1)
	 call pgsci(0)
	 call pgline(ndes,x,y1)
	 if(ihe3 .gt. 1)call pgline(ndes,x,y2)
	 if(ihe4 .gt. 1)call pgline(ndes,x,y3)
	 if(ili7 .gt. 1)call pgline(ndes,x,y4)	
	 if(ic12 .gt. 1)call pgline(ndes,x,y5)
	 if(in14 .gt. 1)call pgline(ndes,x,y6)
	 if(iw .gt. 1)  call pgline(ndes,x,y7)	
	 call pgsci(1)
	endif
 
	do i=1,ndes
	 xdbl=(x(i)*mstar)**(2./3.)
	 xdbl=max(mc(1),min(xdbl,mc(nc)))
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,xdbl,lx,xchim,dxchim)
	 call chim_gram_3(xchim,dxchim,nuc_m)
	 y1(i)=xchim(ih)
	 if(ihe3 .gt. 1)y2(i)=xchim(ihe3)
	 if(ihe4 .gt. 1)y3(i)=xchim(ihe4)
	 if(ili7 .gt. 1)y4(i)=xchim(ili7)
	 if(ic12 .gt. 1)y5(i)=xchim(ic12)
	 if(in14 .gt. 1)y6(i)=xchim(in14)
	 if(iw .gt. 1)then
	  call inter_3('mu',bp,q,qt,n,knot,xdbl,f,dfdx,r2,m23)
	  if(f(3) .gt. 0.)y7(i)=xchim(iw)/f(3)
	 endif
	 y7(1)=y7(2)
c	 write(6,2000)x(i),y1(i),y2(i),y3(i),y4(i),y5(i),y(6),y(7)
	enddo
 
        do i= 1, ndes	! Calcul de Y
           y2(i)=y6(i)+y3(i)
        end dO
	do i= 1, ndes	! Calcul de z
           y3(i)=1.d0-y1(i)-y2(i)
        end do
 
	ym2=0
	ym4=0
	ym5=0	
	ym6=0
	ym7=0	
	do i=1,ndes
	 ym2=max(ym2,y2(i))
	 ym4=max(ym4,y4(i))	
	 ym5=max(ym5,y5(i))
	 ym6=max(ym6,y6(i))
	 ym7=max(ym7,y7(i))	 	 	
	enddo
c	do i=1,ndes
c	 if(ym2 .ne. 0.)y2(i)=y2(i)/ym2
c	 if(ym4 .ne. 0.)y4(i)=y4(i)/ym4	
c	 if(ym5 .ne. 0.)y5(i)=y5(i)/ym5	
c	 if(ym6 .ne. 0.)y6(i)=y6(i)/ym6
c	 if(ym7 .ne. 0.)y7(i)=y7(i)/ym7	
c	 write(6,2000)x(i),y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y7(i)	  	
C	ENDdo
c	pause
 
	lteff=log10(teff)
	lx=n
	call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(n),lx,f,dfdx)	
	llext=log10(f(4))*1.5
	
	call pgvsize(xleft,xleft+ld,ybot,ybot+h)	
	call pgwindow(xmin,xmax,ymin,ymax)
	call pgsci(2)   ! X
	call pgline(ndes,x,y1)
	call pgsci(1)	! Y
	call pgline(ndes,x,y2)
	call pgsci(7)   ! Z
        call pgline(ndes,x,y3)
 
	call pgsci(1)	
	efface=.true.
 
c	pause
 
c	P, T, R, L en fonction de m/Mstar
 
	if(efface)then
	 call pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
	 call pgwindow(xmin1,xmax1,ymin1,ymax1)
	 call pgsci(0)
	 call pgline(ndes,ym,yp)	
	 call pgline(ndes,ym,yt)	
	 call pgline(ndes,ym,yr)
	 call pgline(ndes,ym,yl)
	endif
	
	dq=(n-1.d0)/dble(ndes-1)	!pour les variables	
	do i=1,ndes
	 qx(i)=dq*(i-1)+1.
	 qd(i)=qx(i)
	enddo
	
	call pgsci(1)
	
	do i=1,ndes
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,qx(i),lx,f,dfdx)
	 yp(i)=f(1)
	 yt(i)=f(2)	
	 yr(i)=f(3)		
	 yl(i)=f(4)
	 ym(i)=f(5)	 		
	enddo
	pmax=0
	tmax=0
	rmax=0
	lmax=0
	mmax=0
	do i=1,ndes
	 pmax=max(pmax,yp(i))
	 tmax=max(tmax,yt(i))	
	 rmax=max(rmax,yr(i))	
	 lmax=max(lmax,yl(i))
	 mmax=max(mmax,ym(i))
	enddo	
	do i=1,ndes
	 yp(i)=yp(i)/pmax	
	 yt(i)=yt(i)/tmax	
	 yr(i)=yr(i)/rmax	
	 yl(i)=yl(i)/lmax
	 ym(i)=ym(i)/mmax
	enddo
	 	
	call pgsci(2)
	call pgline(ndes,ym,yp)
	call pgsci(3)
	call pgline(ndes,ym,yt)
	call pgsci(5)
	call pgline(ndes,ym,yr)
	call pgsci(6)
	call pgline(ndes,ym,yl)
	call pgsci(1)	
	
c	diagramme HR
 
 
	call pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot,ybot+h)
	call pgwindow(xmin3,xmax3,ymin3,ymax3)
 
	call pgmove(lteffp,llextp)
	call pgdraw(lteff,llext)
 
	lteffp=lteff
	llextp=llext
	
c	zones radiatives, convectives
 
 
c	print*,lim,(lconv(i),i=1,lim)
c	print*,(m_zc(i),i=1,lim),age
 
	pas=pas+1.
	x4(1)=pas
	x4(2)=pas
 
	if(lim .gt. 0)then	!lim=0: modele totalement radiatif pas de dessin	
	 call pgvsize(xleft+ld+dl,xleft+2.*ld+dl,ybot+h+dh,ybot+2.*h+dh)
	
c	 print*,lim,(lconv(i),i=1,lim)
c	 write(6,2000)(m_zc(i),i=1,lim),mstar,mtot
 
	 call pgwindow(xmin4,xmax4,ymin4,ymax4)
	 if(lim .eq. 1 .and. .not.lconv(1) .and.
     1	m_zc(1)/mstar .ge. .99)then!modele totalement convectif
	  y4(1)=0.
	  y4(2)=1.
c	  call pgsls(4)		!dotted
	  call pgline(2,x4,y4)
	 else
	  i=1
	  do while(i .le. lim)
	   if(lconv(i))then	!debut de ZC
	    y4(1)=m_zc(i)/mstar
	    if(i .eq. lim)then		!ZC externe
	     y4(2)=1.
	    else
	     y4(2)=m_zc(i+1)/mstar
	     i=i+1
	    endif	
c	    call pgsls(4)		!dotted
	    call pgline(2,x4,y4)
	   else
	    if(i .eq. 1)then			!fin de ZC centrale
	     y4(1)=0
	     y4(2)=m_zc(1)/mstar
c	     call pgsls(4)		!dotted
	     call pgline(2,x4,y4)
	    else
	     print*,'probleme'
	     call pgend
	     stop
	    endif
 
	   endif
	   i=i+1
	  enddo	!while
 
	 endif
	endif
 
	
	return
 
	
	end
 
