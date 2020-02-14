c...................
	function bnu(tr,wl)
	include 'param.dec'

	wlcm=wl*1.e-4
	xnu=cluz/wlcm
         dd=h*cluz/wlcm/bk/tr
            if (dd.gt.70.) then
              bnu=2.*h*cluz/(wlcm**3)*exp(-dd)
            elseif (dd.lt.1.e-10) then
              bnu=2.*h*cluz/(wlcm**3)/dd
            else
              bnu=2.*h*cluz/(wlcm**3)/
     *          (exp(dd)-1.)
            endif
	return
	end
c..............................................
	function fz(z)
	common/alt1/fz0,scale,rdcm,deltar
	
	fz1a=exp(-(z**2)/2.)
	fz1b=exp(-(z**2*(1.+deltar)**2/2.))
	fz1=fz1a+fz1b
	fz2=sqrt(z**2+(rdcm/scale)**2)
	fz3=fz1*fz2
	fz=1.0-fz3/fz0
	return
	end
c....................................................
	subroutine wcorr(t)
c	implicit double precision(a-h,o-z)
c	incluye longitudes de onda alrededor del maximo de
c	la function de Planck para asegurar que las integrales
c	estan bien
	parameter(nmax=4000)
        common /ww2/wl0(nmax),nwl0
        common /ww/wl(nmax),nwl

		nwl=nwl0
		do i=1,nwl
		wl(i)=wl0(i)
		enddo
	l=nwl+1
	wmax=4.8e3/t
	dw=wmax/20.
	dw=0.05
	xlm=alog10(wmax)
	wl(l)=wmax
	l=l+1
	i=2
	do while(wl(l-1).lt.wl0(nwl0))
	wl(l)=10.**(xlm+(i-1)*dw)
	l=l+1
	i=i+1
	enddo	
	i=2
	do while(wl(l-1).gt.1.e-2*wl0(1))
        wl(l)=10.**(xlm-(i-1)*dw)
        l=l+1
	i=i+1
        enddo
	nwl=l-1
	if(nwl.gt.nmax)stop 'wcorr'
	call xordr(wl,nwl)
	call una(wl,nwl)
	if(wl(1).lt.1.e-8) then
		do i=1,nwl-1
		wl(i)=wl(i+1)
		enddo
		nwl=nwl-1
	endif
c	write(6,*)nwl
c	write(6,100)(wl(i),i=1,nwl)
100     format(5f15.8)
	return
	end
c------------------------
       subroutine una(x,n)
c	implicit double precision(a-h,o-z)
c      elimina puntos iguales (hasta 7 digitos)
       dimension x(4000),x1(4000)
       x1(1)=x(1)
      l=2
                do i=2,n
                if(abs(x(i)-x(i-1)).gt.1.e-8)then
                x1(l)=x(i)
                l=l+1
                endif
                enddo
                n=l-1
              do i=1,n
               x(i)=x1(i)
               enddo
       return   
      end
c.................................
      SUBROUTINE XORDR (X,N)
c	implicit double precision (a-h,o-z)
      DIMENSION X(4000)
      DO 1 K=1,N
      XX=RMINIM(X,K,N,I)
      X(I)=X(K)
1     X(K)=XX
      RETURN
      END
c.................................
      FUNCTION RMINIM (A,I1,I2,I3)
c	implicit double precision (a-h,o-z)
      DIMENSION A(4000)
      RMINIM=A(I1)
      DO 3 J=I1,I2
      IF (A(J)-RMINIM) 4,4,3
    4 RMINIM=A(J)
      I3=J
    3 CONTINUE
      RETURN
      END
c........................................

c	::::::::::::::
c	zeroin.f
c	::::::::::::::
      function zeroin(ax, bx, f, tol)
c	implicit real*8 (a-h,o-z)
c	compute eps, the relative machine precision
	common/band/iwrr
      external f

	
      eps = 1.0
   10 eps = eps / 2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0) goto 10
    7 a = ax
      b = bx
      fa = f(a)
c	type *,'fa=',fa
      fb = f(b)
c	type *,'fb=',fb
      contr = fa * fb
      if (contr .gt. 0.) then
      write(unit=*, fmt=*) 'problemas en Zeroin'
      write(unit=*, fmt=*) k, 'a=', ax, 'b=', bx
      write(unit=*, fmt=*) 'fa=', fa, 'fb=', fb
	iwrr=1
      go to 100
      end if
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (abs(fc) .ge. abs(fb)) goto 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
   40 tol1 = ((2.0 * eps) * abs(b)) + (0.5 * tol)
      xmm = 0.5 * (c - b)
      if (abs(xmm) .le. tol1) goto 90
      if (fb .eq. 0.0) goto 90
      if (abs(e) .lt. tol1) goto 70
      if (abs(fa) .le. abs(fb)) goto 70
      if (a .ne. c) goto 50
      s = fb / fa
      p = (2.0 * xmm) * s
      q = 1.0 - s
      goto 60
   50 q = fa / fc
      r = fb / fc
      s = fb / fa
      p = s *((((2.0 * xmm) * q) * (q - r))-((b - a) * (r - 1.0)))
      q = ((q - 1.0) * (r - 1.0)) * (s - 1.0)
   60 if (p .gt. 0.0) q = - q
      p = abs(p)
      if ((2.0 * p).ge.(((3.0 * xmm) * q)-abs(tol1 * q))) goto 70
      if (p .ge. abs((0.5 * e) * q)) goto 70
      e = d
      d = p / q
      goto 80
   70 d = xmm
      e = d
   80 a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1,xmm)
      fb = f(b)

      if ((fb * (fc / abs(fc))) .gt. 0.0) goto 20
      goto 30
	
   90 zeroin = b
100      return 
      end
c....................................
c....................................
      real function linter(x0, x, y, i0,nstar)
        parameter (n=1206)
      dimension x(n), y(n)
	common/tstar1/tstar

      if (x(2) .lt. x(1)) then
c       x en orden decreciente
      do while (x0 .lt. x(i0))
      i0 = i0 + 1
      end do
      else
c       x en orden creciente
	if (x0.gt.x(nstar)) then
	go to 2
	endif
	if (x0.le.x(1)) then
	linter=y(1)-4.*(x0-x(1))-
     *  1.68*(3720./tstar)*((1./(10.**x0))-(1./(10.**x(1))))
	return
	endif
      do while (x0 .gt. x(i0))
      i0 = i0 + 1
      end do
      end if

      if (i0 .gt. 1) i0 = i0 - 1
c        if ((i0+1).gt.nstar) i0=nstar-1
2	continue
c	if (x0.ge.x(nstar)) then
c	i0=nstar
	 if (x0.ge.x(1202)) then
       i0=1202

c	linter = -3*(x0-x(i0))+y(i0)
	linter = -2.7*(x0-x(i0))+y(i0)
	else
      linter =
     &  ((y(i0) * (x(i0 + 1) - x0)) + (y(i0 + 1) * (x0 - x(i0))))
     & / (x(i0 + 1) - x(i0))
	endif
      return
      end

c........................................................................
        function trapz(x,y,nn)
        parameter(n=50)
        dimension x(n),y(n)

        trapz=0.
        if (nn.eq.1) then
        trapz=y(1)
        else
        do i=2,nn
        trapz=trapz+0.5*(x(i)-x(i-1))*(y(i)+y(i-1))
        enddo
        endif
	return
	end

c------------------------------

	function plank(nu,T)
c	implicit real*8 (a-h,o-z)
	implicit none
	real*4 dummy,dd
	real*4 plank
	real*4 nu,T
	real*4 h,k,c,Pi
	parameter (h=6.6260755d-27,  k=1.380658d-16	,  c=2
     1    .99792458d10,  Pi=3.14159 	)
	
         dd=h*nu/(k*T)
            if (dd.gt.70.) then
c              plank=((2.*h*(nu**3))/(c**2))*exp(-dd)
			   plank=(((2.*h*(nu**2))/(c**2))*nu)*exp(-dd)
c			  write(*,*)'left',(((2.*h*(nu**2))/(c**2))*nu)
c			  write(*,*)'plank',plank
c			  write(*,*)'nu**3',nu**3
c			  write(*,*)'nu**2',nu**2
c			  write(*,*)'h',h
c			  write(*,*)'h*(nu**3)',h*(nu**3)
            elseif (dd.lt.1.e-10) then
              plank=(((2.*h*nu**2)/(c**2))*nu)/dd
c			  write(*,*)((2.*h*nu**3)/(c**2))
            else
              plank=(((2.*h*nu**2)/(c**2))*nu)/
     *          (exp(dd)-1.)
c	 write(*,*)((2.*h*nu**3)/(c**2))
            endif
	return

	end function plank
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	function fscat(gg,mu)
	implicit none
	real*4 gg,mu
	real*4 fscat
	fscat=(1.-gg**2)/(1.+gg**2-2*gg*mu)**(1.5)
c	write(6,*)gg,mu,fscat,1.-gg**2,(1.+gg**2-2*gg*mu)**1.5
	return
	end 
	
c----------------
c........................................

c	::::::::::::::
c	zeroin.f
c	::::::::::::::
      function zeroind(ax, bx, f, tol)
 	implicit real*8 (a-h,o-z)
c	compute eps, the relative machine precision
	common/band/iwrr
      external f

	
      eps = 1.0
   10 eps = eps / 2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0) goto 10
    7 a = ax
      b = bx
      fa = f(a)
c	type *,'fa=',fa
      fb = f(b)
c	type *,'fb=',fb
      contr = fa * fb
      if (contr .gt. 0.) then
      write(unit=*, fmt=*) 'problemas en Zeroin'
      write(unit=*, fmt=*) k, 'a=', ax, 'b=', bx
      write(unit=*, fmt=*) 'fa=', fa, 'fb=', fb
	iwrr=1
      go to 100
      end if
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (abs(fc) .ge. abs(fb)) goto 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
   40 tol1 = ((2.0 * eps) * abs(b)) + (0.5 * tol)
      xmm = 0.5 * (c - b)
      if (abs(xmm) .le. tol1) goto 90
      if (fb .eq. 0.0) goto 90
      if (abs(e) .lt. tol1) goto 70
      if (abs(fa) .le. abs(fb)) goto 70
      if (a .ne. c) goto 50
      s = fb / fa
      p = (2.0 * xmm) * s
      q = 1.0 - s
      goto 60
   50 q = fa / fc
      r = fb / fc
      s = fb / fa
      p = s *((((2.0 * xmm) * q) * (q - r))-((b - a) * (r - 1.0)))
      q = ((q - 1.0) * (r - 1.0)) * (s - 1.0)
   60 if (p .gt. 0.0) q = - q
      p = abs(p)
      if ((2.0 * p).ge.(((3.0 * xmm) * q)-abs(tol1 * q))) goto 70
      if (p .ge. abs((0.5 * e) * q)) goto 70
      e = d
      d = p / q
      goto 80
   70 d = xmm
      e = d
   80 a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1,xmm)
      fb = f(b)

      if ((fb * (fc / abs(fc))) .gt. 0.0) goto 20
      goto 30
	
   90 zeroin = b
100      return 
      end
