c
c*****************************************************************************
c
	subroutine f4n(t,vol,snn,nspe)
c
c     derivatives with respect to number abundances.
c     calls f4 of MHD package and prepares quantities not provided by f4
c     (degeneracy-related stuff).
c
c === all units c.g.s. and degrees Kelvin
c
      implicit real*8 (a-h,o-z)
	implicit integer(i-n)
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
      dimension snn(nspe)
c
      common /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      common /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      common /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      common /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      common /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )
 
	logical ok
	common /marche/ok
 
c
      do 5 is=1,nspes
 5    sn(is) = snn(is)
c
      ckt   = ck * t
c
	call neweta(t,vol,sn(ise),ier)
	if(ier.ne.0) then
         write(6,*) 'error in s/r f4n. failure in neweta. f4 not called'
	 write(6,*)'appel a ETAT_EFF'
	 ok=.false.	!appel a etat_eff
	 return
	end if
c
      call f4mhd(t,vol)
c
      return
      end
