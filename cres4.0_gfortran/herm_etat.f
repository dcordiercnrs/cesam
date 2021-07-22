 
c*******************************************************************
 
	subroutine herm_etat(f,dfx,lx,fz,dfzx,d2fzxy,herx,hery,
     1	n,tx,ly,dfy,dfxy,dfzy,d2fzy2,deriv,ind)
 
c	interpolation avec une base de hermite adapte a l'equation d'etat
 
c	Auteur: P. Morel Observatoire de NICE
c	version du 22 11 91
 
 
c entree
c	f : fonction en (x,y)
c	dfx, dfy : derivee
c	dfxy : d2f/dxdy
c	lx,ly : localisation
c	deriv=.true. : calcul des derivees secondes
c	n : nombre de tables
c	herx,hery: bases de hermite en X et Y
c	ind : indices en y
 
c sortie
c	fz, dfzx, dfzy, d2fzxy, d2fzy2 : fonction et derivees en z
 
	implicit none
 
	integer lx,ly,i,j,kx,ky,n,l,indice,ind(2),tx
 
	real*8  f(1),dfx(1),fz(1),dfzx(1),herx(0:2,2,2),
     1	d2fzxy(1),dfy(1),dfzy(1),d2fzy2(1),hery(0:2,2,2),dfxy(1)
 
	logical deriv,init
	data init/.true./
 
2000	format((1x,1p8d10.3))
 
	do l=1,n
	 fz(l)=0.
	 dfzx(l)=0.
	 dfzy(l)=0.
	 d2fzxy(l)=0.
	 d2fzy2(l)=0.
	 do i=1,2
	  kx=lx+i-1
	  do j=1,2
	   ky=ind(i)+j-1
	   indice=n*(tx*(ky-1)+kx-1)+l
c	   write(6,*)'l,i,j,kx,ky,indice',l,i,j,kx,ky,indice
c	   write(6,2000)f(indice),dfx(indice),dfy(indice),dfxy(indice)
	   fz(l)=fz(l)+f(indice)*herx(0,1,i)*hery(0,1,j)+
     1	     dfx(indice)*herx(0,2,i)*hery(0,1,j)+
     2	     dfy(indice)*herx(0,1,i)*hery(0,2,j)+
     3	    dfxy(indice)*herx(0,2,i)*hery(0,2,j)
 
	   dfzx(l)=dfzx(l)+f(indice)*herx(1,1,i)*hery(0,1,j)+
     1		 dfx(indice)*herx(1,2,i)*hery(0,1,j)+
     2		 dfy(indice)*herx(1,1,i)*hery(0,2,j)+
     3		dfxy(indice)*herx(1,2,i)*hery(0,2,j)
 
	   dfzy(l)=dfzy(l)+f(indice)*herx(0,1,i)*hery(1,1,j)+
     1		 dfx(indice)*herx(0,2,i)*hery(1,1,j)+
     2		 dfy(indice)*herx(0,1,i)*hery(1,2,j)+
     3		dfxy(indice)*herx(0,2,i)*hery(1,2,j)
 
	  if(deriv)then
 
	   d2fzxy(l)=d2fzxy(l)+f(indice)*herx(1,1,i)*hery(1,1,j)+
     1		     dfx(indice)*herx(1,2,i)*hery(1,1,j)+
     2		     dfy(indice)*herx(1,1,i)*hery(1,2,j)+
     3		    dfxy(indice)*herx(1,2,i)*hery(1,2,j)
 
	   d2fzy2(l)=d2fzy2(l)+f(indice)*herx(0,1,i)*hery(2,1,j)+
     1		     dfx(indice)*herx(0,2,i)*hery(2,1,j)+
     2		     dfy(indice)*herx(0,1,i)*hery(2,2,j)+
     3		    dfxy(indice)*herx(0,2,i)*hery(2,2,j)
	  endif		!deriv
	  enddo		!j
	 enddo		!i
	enddo 		!l
c	write(6,*)'fz,dfzdx,dfzdy,d2fzxy,d2fzy2'
c	write(6,2000)(fz(i),i=1,n)
c	write(6,2000)(dfzx(i),i=1,n)
c	write(6,2000)(dfzy(i),i=1,n)
c	write(6,2000)(d2fzxy(i),i=1,n)
c	write(6,2000)(d2fzy2(i),i=1,n)
c	write(6,*)' '
 
	return
 
	end
