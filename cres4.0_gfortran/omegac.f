 
c*****************************************************************************
c
	function omegac(i,iz)
c  calculates statistical weight of i-th ionization stage of element
c  with number iz
      implicit double precision (a-h,o-z)
	implicit integer(i-n)
      common/hvomeg/ iom(26),iom1(20)
      common/hvomcl/ iomfll
c
      save
c
      if(i.le.1.and.iz.ge.19) go to 20
      if(i.eq.iz) go to 15
      omegac=iom(iz-i)
      return
c
c  statistical weight for fully ionized atom.
c  before 5/1/84 was always set to 15, for some bizarre reason.
c  for transition introduce flag iomfll so that iomfll = 0 corresponds
c  to old situation, and iomfll .ne. 0 gives correct value omegac = 1.
c
   15 omegac=1
      if(iomfll.eq.0) omegac=15
      return
   20 omegac=iom1(2*(iz-19)+i+1)
      return
      end
