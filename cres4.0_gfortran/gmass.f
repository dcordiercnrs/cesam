 
 
c******************************************************************8
 
      function gmass(x,z,amoles,eground,fracz,frac)
 
c	adaptation a CESAM3 de la routine gmass du package OPAL_EOS
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
	implicit none
 
	integer	i
	
	real*4	x,z,amoles,eground,fracz,anum(6),frac(7),xc,xn,xo,xne,
     1	xc2,xn2,xo2,xh,xtot,gmass,anume,xhe,xne2,amas(7),eion(7)
	data (anum(i),i=1,6)/10.,8.,7.,6.,2.,1./
	data (eion(i),i=1,6)/-3388.637,-1970.918,-1431.717,-993.2303,-76.2315,-15.29409/
 
 
      xc=0.247137766
      xn=.0620782
      xo=.52837118
      xne=.1624188
      amas(7)=1.0079
      amas(6)=4.0026
      amas(5)=12.011
      amas(4)=14.0067
      amas(3)=15.9994
      amas(2)=20.179
      amas(1)=0.00054858
      fracz=z/(xc*amas(5)+xn*amas(4)+xo*amas(3)+xne*amas(2))
      xc2=fracz*xc
      xn2=fracz*xn
      xo2=fracz*xo
      xne2=fracz*xne
      xh=x/amas(7)
      xhe=(1.-x -z)/amas(6)
      xtot=xh+xhe+xc2+xn2+xo2+xne2
      frac(6)=xh/xtot
      frac(5)=xhe/xtot
      frac(4)=xc2/xtot
      frac(3)=xn2/xtot
      frac(2)=xo2/xtot
      frac(1)=xne2/xtot
      eground=0.0
      amoles=0.0
      do i=1,6
      eground=eground+eion(i)*frac(i)
      amoles=amoles+(1.+anum(i))*frac(i)
      enddo
      anume=amoles-1.
      gmass=anume*amas(1)
      do i=2,7
      gmass=gmass+amas(i)*frac(i-1)
      enddo
 
      return
 
      end
