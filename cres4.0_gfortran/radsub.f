 
c***********************************************************************
 
	subroutine radsub(t6,density,moles,tmass)
 
 
c	adaptation a CESAM3 de la routine radsub du package OPAL_EOS
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
 	implicit none
 	
 	real*4	t6,density,moles,tmass,k,molenak,Na,unitf,unitfold,c,sigma,
     1	sigmac,sigmacc,aprop
 
c	data Na/6.0221367e+23/,k/1.38065e-16/,unitf/0.9648575/,
c     1	unitfold/0.9652/,c/2.9979245e+10/,sigma/5.67051e-5/
c     2	sigmac/1.8914785e-15/,sigmacc/1.8914785e-3/,aprop/83.14511/
 
	real*4	rat,pr,er,sr,revise,fixerror,st,chir,chitt,pt,et,
     1	cvtt,gam3pt,gam1t,gam2pt,dedrhoat
	 	
	integer mx,mv,nr,nt
	parameter (mx=5,mv=10,nr=121,nt=167)
 
	real*4	esact,eos(mv)	
	common/e/esact,eos
 
	integer iri(10),index(10),nta(nr)
	real*4	 zz(mx)
	common/b/iri,index,nta,zz
 
	data Na/6.0221367e+23/,k/1.38065e-16/,unitf/0.9648575/,
     1	unitfold/0.9652/,c/2.9979245e+10/,sigma/5.67051e-5/
     2	sigmac/1.8914785e-15/,sigmacc/1.8914785e-3/,aprop/83.14511/
	save
 
c	Physical constants
c       Na=6.0221367e+23
c       k =1.380658e-16 !   erg/degree K
c       Na*k=6.0221367E+23*1.380658e-16 erg/degree K=8.314511E+7 erg/degree K
c           =8.314511e+7*11604.5 erg/eV=0.9648575E+12 erg/eV
c           Define unitf= Na*k/e+12=0.9648575
c           unitf=0.9648575  ! obtained with latest physical constants
c           unitfold=0.9652   ! old units- still used in the EOS code
c           In these units energy/density is in units of Mb-CC/gm
c           Pressure is in units of E+12 bars=Mb
c       sigma is the Stefan-Boltzmann constant
c       sigma=5.67051E-5 !   erg /(s*cm**2*K**4)
c       c=2.99792458E+10 !   cm/sec
 
c     rat=sigma/c    ! dyne/(cm**2*K**4)
 
c     rat=rat*1.e+24  !  Convert degrees K to units 10**6 K (T replaced with T6)
      rat=sigmacc
 
      pr=4./3.*rat*t6**4   ! Mb
      er=3.*pr/density   ! Mb-cc/gm
      sr=4./3.*er/t6   ! Mb-cc/(gm-unit T6)
      revise=unitf/unitfold
      eos(iri(1))=eos(iri(1))*revise
      eos(iri(2))=eos(iri(2))*revise
      eos(iri(3))=eos(iri(3))*revise
      eos(iri(4))=eos(iri(4))*revise
      eos(iri(5))=eos(iri(5))*revise
      pt=eos(iri(1))+pr
      et=eos(iri(2))+er
      fixerror=1./(.0114045*.0116045) !4/14/96 Corrects for earlier
c        multiplication by .0114045 where divide by .0116045 needed
c        Converts eV to T6
      st=eos(iri(3))*fixerror+sr
      chir=eos(iri(6))*eos(iri(1))/pt
      chitt=(eos(iri(1))*eos(iri(7))+4.*pr)/pt
c     gam1t(jcs,i)=(p(jcs,i)*gam1(jcs,i)+4./3.*pr)/pt(jcs,i)
c     gam2pt(jcs,i)=(gam2p(jcs,i)*p(jcs,i)+4.*pr)/pt(jcs,i)
c     gam3pt(jcs,i)=gam1t(jcs,i)/gam2pt(jcs,i)
      molenak=moles*aprop  ! Mb-cc/unit T6
      cvtt=(eos(iri(5))*molenak/tmass+4.*er/t6)
      gam3pt=pt*chitt/(cvtt*density*t6)
      gam1t=chir+chitt*gam3pt
      gam2pt=gam1t/gam3pt
 
c     normalize cvt to 3/2 when gas is ideal,non-degenerate,
c     fully-ionized, and has no radiation correction
c     cvt=(eos(5)*molenak/tmass+4.*er/t6)
c    x  /molenak
      dedrhoat=eos(iri(4))-er/density
      eos(iri(1))=pt
      eos(iri(2))=et
      eos(iri(3))=st
      eos(iri(4))=dedrhoat
      eos(iri(5))=cvtt
      eos(iri(6))=chir
      eos(iri(7))=chitt
      eos(iri(8))=gam1t
      eos(iri(9))=gam2pt
      eos(iri(10))=gam3pt
 
      return
 
      end
