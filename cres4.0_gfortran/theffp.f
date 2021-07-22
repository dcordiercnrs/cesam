 
c*******************************************************************
 
	subroutine theffp(p10,t10,chem,ychem,f)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
 
	implicit real*8 (a-h,o-z)
 
c	implicit real*16 (a-h,o-z)
	integer ihvz
      common/eqscnt/ anh0,ihvz
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     &  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/eqsaux/psi
c
c================= output =========================================
      common/eos/ald8,cp8,psi8,rhp8,rht8,vlro8,vmol8,
     *vmy8,vmyp8,vmyt8,vna8,
     *xio8(3),xmol8,e8,uint8,dpdrs8,scgs8,cv8
c==================================================================
c
      umod  = log(10.d0)
      xc=chem
      yc=ychem
      zc=1.-xc-yc
c
      call eqstp(p10,t10,xc,yc,zc,flnew)
c
      psi8   =psi
      vlro8  = log10(rho(1))
      rhp8   =rho(2)/pt(2)
      rht8   = -dlt(1)
      vna8   =dad(1)
      xio8(1)=xii1(1)
      xio8(2)=xii1(2)
      xio8(3)=xii1(3)
      xmol8  = 0.
c>>>>>>>>>>>>>>>>>> no molecules in eff
      rmu    = xc/ah+yc/ahe+zc/az
      vmol8  = 1./rmu
      cp88   =cp(1)
      rgas0  = 8.31434 e 7
      cp8    =cp88*vmol8/rgas0
      emax   =(xc/ah+2.*yc/ahe+zc*anh0)*av
      e8     =ane(1)/emax
      vmy8   =vmol8/(1.+e8)
      fact   =1.0 + e8
      xelp   =(ane(2)/pt(2))/emax
      xelt   =(-ane(2)*pt(3)/pt(2) + ane(3))/emax
      vmyp8  = - xelp/(fact*umod)
      vmyt8  = - xelt/(fact*umod)
      dpdrs8 = gm1
      f=10.d0**flnew
 
	return
 
	end
