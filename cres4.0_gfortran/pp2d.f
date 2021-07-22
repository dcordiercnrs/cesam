 
c************************************************************************
 
	subroutine pp2d(x,xt,xr,nx,nrx,mx,xx,knotx,lx,dfxydx,f,s,fxy,init,
     1	        y,yt,yr,ny,nry,my,yy,knoty,ly,dfxydy)
 
c	interpolation 2D equivalente a sbsp2d mais on calcule avec les
c	polynomes par morceaux d'ou un gain de temps justifie a partir
c	de quelques nx*ny appels
 
c	interpolation dans le tableau f(nx,ny) par produit tensoriel
c	de splines polynomiales d'ordres mx{y}>1, au point
c	(xx,yy) : x(1).le.xx.le.x(nx), y(1).le.yy.le.y(ny), on aura aussi
c	 xt(lx).le.xx.lt.xt(lx+1)	yt(ly).le.yy.lt.yt(lt+1)
c	en entree prendre lx{y}=mx{y}
c	au premier appel, init=.false. il y a initialisation :
c	calcul de knotx=2mx+nx, formation du tableau xt(knotx). idem pour y
c	!! ATTENTION le tableau f(nx,ny) initial est remplace par le
c	tableau des coefficients des B-splines en xy!!
c	ce calcul etant fait selon de Boor p.342
 
c	ax et ay doivent etre dimensionnes au moins a 	max(nx*nx,ny*ny)
c	l et m	"	"	"	"	"	max(nx,ny)
 
c	les tables de travail xr, yr sont conservees car on peut faire appel
c	sequentiellement a des tables diverses
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer dimab,dimxy,dimaxy,dimqxy,dimdxy
	parameter (dimab=64, dimxy=6)
	parameter (dimaxy=dimab**2, dimqxy=dimxy+1, dimdxy=dimxy**2)
 
	integer nx,ny,knotx,knoty,lx,ly,mx,my,l(dimab),m(dimab),i,j,k,
     1	llx,lly,ix,iy,sx,sy,nrx,nry
 
	real*8 x(1),y(1),f(1),xt(1),yt(1),s(1),xx,yy,fxy,dfxydx,dfxydy,
     1	ax(dimaxy),ay(dimaxy),qx(dimqxy),qy(dimqxy),dx(dimdxy),
     2	dy(dimdxy),pivot,d2fdxy,deriv,fac(dimxy),xx1,yy1,
     3	xr(1),yr(1)
 
	logical init
 
	if(init)goto 100
 
c	initialisations
 
	if(nx .gt. dimab .or. mx .gt. dimxy)goto 200
	if(ny .gt. dimab .or. my .gt. dimxy)goto 200
	
	do i=1,dimaxy
	 ax(i)=0.
	 ay(i)=0.
	end do
 
	lx=mx
	ly=my
 
c	formation des tableaux xt, yt et points de raccord
 
	call snoein(x,xt,nx,mx,knotx)
	call snoein(y,yt,ny,my,knoty)
	nrx=nx-mx+2
	nry=ny-my+2
	do ix=1,nrx
	 xr(ix)=xt(ix+mx-1)
	enddo
	do iy=1,nry
	 yr(iy)=yt(iy+my-1)
	enddo
 
c	ay <- matrice 1-D des interpolations en y
 
	do j=1,ny
	call slinf(y(j),yt,knoty,ly)
	call bval(y(j),yt,my,ly,qy)
	 do k=1,my
	  ay(ny*(ly-my+k-1)+j)=qy(k)
	 enddo	!k
	enddo	!j
 
c	ax <- transposee de ay, inversion de ax et ay=f*ax
 
	call strans(ay,ny,ax)
 
	call minv(ax,ny,pivot,l,m)
	
	if(pivot .eq. 0.)goto 220
	call smatpr(f,nx,ny,ax,ny,ay)
 
c	ax <- matrice 1-D des interpolations en x
 
	do i=1,dimaxy
	 ax(i)=0.
	end do
	do j=1,nx
	call slinf(x(j),xt,knotx,lx)
	call bval(x(j),xt,mx,lx,qx)
	 do k=1,mx
	  ax(nx*(lx-mx+k-1)+j)=qx(k)
	 enddo	!k
	enddo	!j
 
c	ax=inverse de ax et f=ax*ay
 
	call minv(ax,nx,pivot,l,m)
	if(pivot .eq. 0.)goto 230
	call smatpr(ax,nx,nx,ay,ny,f)
 
c	calcul des coefficients des pp. : 12-21 de schumaker
 
c	de facon a optimiser l'algoritme de horner avec derivee
c	dans s, les coefficients du polynome (derives/fac(.)) sont classes
c	par puissance decroissante :
c	fy3x3, fy3x2, fy3x, fy3; fy2x3, fy2x2, fy2x, fy2;
c	fyx3, fyx2, fyx,fy; fx3, fx2, fx, f
 
c	s(ix,iy,lx,ly)=s(mx*(my*((nrx-1)*(ly-1)+lx-1)+iy-1)+ix)
c	il faut ajouter 1 a l'ordre de la derivation, ainsi
c	pour fy2x on aura ix=mx-2, iy=my-3 et, pour l'ordre de derivation i,j
c	correspondant aux indices k=i+1, l=j+1 pour bvald, on aura
c	ix=mx-k+1, iy=my-l+1
 
	fac(1)=1.	!calcul de (i-1)!
	fac(2)=1.
	do i=3,dimxy
	 fac(i)=fac(i-1)*float(i-1)
	enddo
 
	llx=mx
	lly=my
	do lx=1,nrx-1
	 call slinf(xr(lx),xt,knotx,llx)
	 call bvald(xr(lx),xt,mx,llx,mx,dx)
	 do ly=1,nry-1
	  call slinf(yr(ly),yt,knoty,lly)
	  call bvald(yr(ly),yt,my,lly,my,dy)
	  do ix=1,mx
	   do iy=1,my
	    deriv=0.
	    do sx=1,mx
	     do sy=1,my
	      deriv=deriv+f(nx*(lly-my+sy-1)+llx-mx+sx)*
     1	dx(mx*(sx-1)+ix)*dy(my*(sy-1)+iy)
	     enddo	!sy
	    enddo	!sx
	    s(mx*(my*((nrx-1)*(ly-1)+lx-1)+my-iy)+mx-ix+1)=
     1	deriv/fac(ix)/fac(iy)
	   enddo	!iy
	  enddo	!ix
	 enddo	!ly
	enddo	!lx
 
c	localisation de (xx,yy)
 
 100	if(xx .lt. x(1))then
	 xx1=xr(1)
	 lx=1
c	 write(6,*)'dans pp2d extrapolation pour x=',xx,' .lt. x(1)=',x(1)
	elseif(xx .gt. x(nx))then
	 xx1=xr(nrx)
	 lx=nrx-1
c	 write(6,*)'dans pp2d extrapolation pour x=',xx,' .gt. x(nx)=',x(nx)
	else
	 call slinf(xx,xr,nrx,lx)
	 xx1=xx
	endif
 
	if(yy .lt. y(1))then
	 yy1=yr(1)
	 ly=1
c	 write(6,*)'dans pp2d extrapolation pour y=',yy,' .lt. y(1)=',y(1)
	elseif(yy .gt. y(ny))then
	 yy1=yr(nry)
	 ly=nry-1
c	 write(6,*)'dans pp2d extrapolation pour y=',yy,' .gt. y(ny)=',y(ny)
	else
	 call slinf(yy,yr,nry,ly)
	 yy1=yy
	endif
 
c	interpolation
 
	do iy=1,my	!lly premier indice pour y**(my-iy)
	 lly=mx*(my*((nrx-1)*(ly-1)+lx-1)+iy-1)+1
	 call horder(mx,lly,s,xx1-xr(lx),ax(iy),ay(iy))
	enddo	!iy
 
	call horder(my,1,ax,yy1-yr(ly),fxy,dfxydy)
	call horder(my,1,ay,yy1-yr(ly),dfxydx,d2fdxy)
 
	return
 
 200	write(6,*)'dans pp2d regarder le commentaire sur les dimensions'
	write(6,201)nx,ny,mx,my
 201	format(1x,'nx=',i4,' ou ny=',i3,' ou mx=',i3,' ou my=', i3,
     1	' est trop grand')
	stop
 
220	write(6,*)'dans l inversion 1 de ax'
	stop
 
230	write(6,*)'dans l inversion 2 de ax'
	stop
 
	end
