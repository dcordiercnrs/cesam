c
c******************************************************************************
c
	subroutine hheion(eah, eahe, eahep, xih, xihe, xihep, secder)
c
c  given ratios between succesive states of ionization, and
c  derivatives, for H, He and He+ in eah(1-10), eahe(1-10),
c  and eahep(1-10), sets the corresponding degrees of ionization,
c  and derivatives, into xih(1-10), xihe(1-10) and xihep(1-10)
c
      implicit double precision (a-h,o-z)
	implicit integer(i-n)
      logical secder
      dimension eah(10), eahe(10), eahep(10), xih(10), xihe(10),
     *  xihep(10), eahs(10), eahes(10), eahepc(10)
c
      if(secder) then
        iders=10
      else
        iders=4
      end if
c
c  set combined ratio for he+
c
      eea=eahep(1)
      if(eea.gt.0) then
        eea1=eahe(1)
        eahepc(1)=eea*eea1
        ii=4
        if(secder) then
          do 15 l=2,4
          do 15 m=l,4
          ii=ii+1
   15     eahepc(ii)=eea*eahe(ii)+eahep(l)*eahe(m)+eahe(l)*eahep(m)+
     .      eea1*eahep(ii)
        end if
        do 20 i=2,4
   20   eahepc(i)=eahep(i)*eea1+eea*eahe(i)
c
      else
c
        do 22 i=1,iders
   22   eahepc(i)=0
c
      end if
c
c  set x h+, x he+, x he++ and derivatives
c
      dnm=1+eah(1)
      xih(1)=eah(1)/dnm
c
c  hydrogen fraction
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      do 25 i=2,iders
   25 eahs(i)=eah(i)/dnm
c
      ii=4
      do 30 l=1,3
      l1=l+1
      eal=eahs(l1)
      xih(l1)=eal/dnm
      if(secder) then
        do 28 m=l,3
        m1=m+1
        ii=ii+1
   28   xih(ii)=(eahs(ii)-2*eal*eahs(m1))/dnm
      end if
   30 continue
c
c  helium fractions
c
      dnm=1+eahe(1)+eahepc(1)
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      do 35 i=2,iders
      eahes(i)=eahe(i)/dnm
   35 eahepc(i)=eahepc(i)/dnm
c
      ii=4
      eeahe=eahe(1)
      eeahep=eahepc(1)
      xihe(1)=eeahe/dnm
      xihep(1)=eeahep/dnm
      do 40 l=2,4
      ealhe=eahes(l)
      ealhep=eahepc(l)
      anmhe=(1+eeahep)*ealhe-eeahe*ealhep
      anmhep=(1+eeahe)*ealhep-eeahep*ealhe
      xihe(l)=anmhe/dnm
      xihep(l)=anmhep/dnm
c
c  second derivatives
c
      if(secder) then
        do 38 m=l,4
        ii=ii+1
        eamhe=eahes(m)
        eamhep=eahepc(m)
c
c  for xi .lt. 1.e-10 zero second derivatives
c
        if(xihe(1).le.1.e-10) then
          xihe(ii)=0
        else
          eamhe=eahes(m)
          eamhep=eahepc(m)
          xihe(ii)=ealhe*eamhep-ealhep*eamhe+
     *      ((1+eeahep)*eahes(ii)-eeahe*eahepc(ii)
     *      -2*anmhe*(eamhe+eamhep))/dnm
        end if
c
        if(xihep(1).le.1.e-10) then
          xihep(ii)=0
        else
          xihep(ii)=ealhep*eamhe-ealhe*eamhep+
     *      ((1+eeahe)*eahepc(ii)-eeahep*eahes(ii)
     *      -2*anmhep*(eamhep+eamhe))/dnm
        end if
   38   continue
      end if
   40 continue
c
      return
      end
