 
c****************************************************************
 
c FICHIER CEFF.FOR
c Paquet de Jorgen Christensen-Dalsgaard sans programme principal
c ordre alphabetique
 
c*****************************************************************************
c
	subroutine bilinc(x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z)
c
c------- performs bilincear interpolation of the fonction f,
c------- given on three arbitray points (x0,y0),(x1,y1),(x2,y2),
c------- where the respective fonction values are z0,z1,z2.
c
      implicit double precision(a-h,o-z)
	implicit integer(i-n)
      x10=x1-x0
      x20=x2-x0
      y10=y1-y0
      y20=y2-y0
      z10=z1-z0
      z20=z2-z0
c
      det=x10*y20-y10*x20
c
      if (det.eq.0) goto 999
c
      dzdx=(z10*y20-z20*y10)/det
      dzdy=(z20*x10-z10*x20)/det
c
      z = z0 + (x-x0)*dzdx + (y-y0)*dzdy
c
      goto 1000
 999  print 8000,x10,x20,y10,y20
      stop
1000  return
8000  format(/' collinear points in s/r bilinc. error stop.',
     . ' x1-x0,x2-x0,y1-y0,y2-y0 = ',/1x,1p4g15.6/)
      end
