 
c******************************************************************
 
	subroutine noeu_dis(x,n,nd,id,m,xt,knot)
 
c	determine la sequence de noeuds de raccord xt(knot)
c	pour une interpolation "optimale"
c	en tenant compte de discontinuites
c	identique a snoein s'il n'y a pas de discontinuites i.e. nd=0
 
c	par B-splines d'ordre m sur la suite strictement
c	croissante x de n points de donnee, cf. de Boor p.219 formule (10)
c	aux limites le polynome d'interpolation s'appuie sur m points de table
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 15 09 92
 
c entrees
c	x : abscisses strictement croissantes
c	n : nombre de points
c	id : indice des discontinuites avec, en sortie id(0)=1, id(nd+1)=n
c	nd : nombre de discontimuites
c	m : ordre des splines
 
c sorties
c	xt : points de table
c	knot : nombre de points de table
 
	implicit none
 
	integer m,knot,n,i,j,ij,nd,id(0:*)
 
	real*8 x(1),xt(1),mm1,eps
	data eps/1.d-14/
 
2000	format((1x,1p8d10.3))
 
c	write(6,*)'NOEU_DIS n,nd,m,id/x',n,nd,m,(id(i),i=1,nd)
c	write(6,2000)(x(i),i=1,n)
 
c	verification de la stricte croissance de la suite des x
 
	do i=1,n-1
	 if(x(i) .ge. x(i+1))then
	  write(6,*)'dans noeu_dis la suite des abscisses n''est pas',
     1	' strictement croissante, en i=',i
	  write(6,*)'nombre de points: ',n
	  write(6,*)'abscisses: '
	  write(6,*)(x(knot),knot=1,i)
	  write(6,*)' '
	  write(6,*)(x(knot),knot=i+1,n)
	  stop
	 endif
	enddo
 
c	pour l'interpolation spline il faut n .ge. m
 
       if(n .ne. 2 .and. n .lt. m)then
          write(6,11)n,m
11        format(1x,'dans noeu_dis n=',i3,' .lt. ',i3,'=m')
          stop
       endif
 
       if(nd .gt. 0)then
         do i=2,nd
           if(id(i)-id(i-1) .lt. m)then
           write(6,*) 'dans noeu_dis on impose m+1 points'
           write(6,*) ' entre 2 discontinuites'
           write(6,*) 'nombre de discontinuites:',nd
           write(6,*) 'indices des discontinuites:',(id(j),j=1,nd)
           stop
           endif
         enddo
 
	 mm1=m-1
	 id(0)=1
	 id(nd+1)=n
	 knot=0
 
	 do i=1,m		!m points de table en 1
	  knot=knot+1
	  xt(knot)=x(1)-eps
	 enddo
 
	 do ij=1,nd+1			!entre chaque discontinuite
	  do i=id(ij-1),id(ij)-m		!entre 2 discontinuites
	   knot=knot+1
	   xt(knot)=0
	   do j=i+1,i+m-1		!moyenne de Schomberg
	    xt(knot)=xt(knot)+x(j)
	   enddo
	   xt(knot)=xt(knot)/mm1
	  enddo		!i
	  do j=1,m			!a chaque discontinuite
	   knot=knot+1
	   xt(knot)=x(id(ij))+eps
	  enddo	!j
	 enddo		!ij
	else
	 write(6,*)'erreur: dans noeu_dis, nb. de discontinuites nd=',nd
	 stop
	endif
 
	return
 
	end
 
