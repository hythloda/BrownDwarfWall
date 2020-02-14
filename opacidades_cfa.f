	subroutine opa_prom_tot(t,ts,xkaps,xchis)
c.	basado en meanopa.f
c.	calcula opacidades promedio de rosseland y Planck
c.	de polvo, a partir de indices de refraccion

c	implicit double precision (a-h,o-z)
	parameter(kmax=100,nmol=16,nelem=30,nmax=4000)
	dimension ppesop(nmax),ppesopr(nmax)
     *		,pesop(nmax),pesopr(nmax)
        dimension g(nmax),gr(nmax)
	dimension xg(nmax),xgr(nmax),xgt(nmax)
	character*24 d_and_t
	character*14 disk

c.	de las longitudes de onda:
        common /ww2/wl0(nmax),nwl0
        common /ww/wl(nmax),nwl
	common/salva/dusto(nmax),alb(nmax)
	common/tsublim/torg,ttroi,tice
	common/salvasil/dustsil(nmax),albsil(nmax)
	common/salvaorg/dustorg(nmax),alborg(nmax)
        common/salvatroi/dusttroi(nmax),albtroi(nmax)
        common/salvaice/dustice(nmax),albice(nmax)


	include 'param.dec'
c	include 'blockdata'

	c=cluz

c    write(disk,'(a14)') "/Users/ccespa/"

c	call fdate(d_and_t)
c	write(*,*)'TIME begin readin: ', d_and_t(12:19)
        if (wl(1).eq.0.) then
        open (1,file=
     1   "/Users/ccespa/Desktop/ccespa/programs/pared/long_mar2004.ent",
     1   status="old")
        read(1,*)idummy
        read(1,*)nwl
	write(6,*)nwl
        read(1,*)idummy
        read(1,*)(wl(i),i=1,nwl)
        close(1)
	endif
c	call fdate(d_and_t)
c	write(*,*)'TIME end readin: ', d_and_t(12:19)


c	 nwl0=nwl
c	do il=1,nwl0
c	wl0(il)=wl(il)
c	enddo		

c	call wcorr(ts)

c	call fdate(d_and_t)
c	write(*,*)'TIME begin dust: ', d_and_t(12:19)

	if (dustsil(100).eq.0.) then
	write(*,*)'EMPIEZO A LEER'
                do i=1,nwl
c					if (i.eq.1) then	
c					call fdate(d_and_t)
c					write(6,*)'TIME begin opasil: ', d_and_t(12:19)
c					end if
c                write(*,*)'comienzo dustsil', i
				dustsil(i)=opasil(wl(i),albsil(i))
c					if (i.eq.1) then	
c					call fdate(d_and_t)
c					write(6,*)'TIME end opasil: ', d_and_t(12:19)
c					end if
					
c						if (i.eq.1) then	
c						call fdate(d_and_t)
c						write(6,*)'TIME begin opaorg: ', d_and_t(12:19)
c						end if
c                write(*,*)'comienzo dustorg', i
				dustorg(i)=opaorg(wl(i),alborg(i))
c							if (i.eq.1) then	
c							call fdate(d_and_t)
c							write(6,*)'TIME end opaorg: ', d_and_t(12:19)
c							end if
c							if (i.eq.1) then	
c							call fdate(d_and_t)
c							write(6,*)'TIME begin opatroi: ', d_and_t(12:19)
c							end if
c                write(*,*)'comienzo dusttroi', i
				dusttroi(i)=opatroi(wl(i),albtroi(i))
c							if (i.eq.1) then	
c							call fdate(d_and_t)
c							write(6,*)'TIME end opatroi: ', d_and_t(12:19)
c							end if				
c							if (i.eq.1) then	
c							call fdate(d_and_t)
c							write(6,*)'TIME begin opaice: ', d_and_t(12:19)
c							end if		
c                write(*,*)'comienzo dustice', i					
				dustice(i)=opaice(wl(i),albice(i))
c							if (i.eq.1) then	
c							call fdate(d_and_t)
c							write(6,*)'TIME end opaice: ', d_and_t(12:19)
c							end if				
				
				enddo
	write(*,*)'TERMINO DE LEER'
	endif

c	call fdate(d_and_t)
c	write(*,*)'TIME end dust: ', d_and_t(12:19)
	
          plancks=0.
          planckts=0.
          xnorpp=0.

c	call fdate(d_and_t)
c	write(*,*)'TIME begin opa: ', d_and_t(12:19)
	
	do i=1,nwl

	opa1=dustsil(i)
	alb1=albsil(i)

	if (t.lt.torg) then
        opa2=dustorg(i)
        alb2=alborg(i)
        else
        opa2=0.
        alb2=0.
        endif

	if (t.lt.ttroi) then
        opa3=dusttroi(i)
        alb3=albtroi(i)
        else
        opa3=0.
        alb3=0.
        endif

	if (t.lt.tice) then
        opa4=dustice(i)
        alb4=albice(i)
        else
        opa4=0.
        alb4=0.
        endif

c	dusto(i)=opa1
c	alb(i)=(opa1*alb1)/dusto(i)

	dusto(i)=opa1+opa2+opa3+opa4
        alb(i)=(opa1*alb1+opa2*alb2+opa3*alb3+opa4*alb4)/dusto(i)

	  wlcm=wl(i)*1.e-4
	  xnu=cluz/wlcm
	dd=h*cluz/wlcm/bk/ts
            if (dd.gt.70.) then
              bnu=2.*h*cluz/(wlcm**3)*exp(-dd)
            elseif (dd.lt.1.e-10) then
              bnu=2.*h*cluz/(wlcm**3)/dd
            else
              bnu=2.*h*cluz/(wlcm**3)/
     *          (exp(dd)-1.)
            endif

          ppesop(i)=bnu/wlcm/wlcm

	wlang=wlcm*1.e8

	  dust=(1.-alb(i))*dusto(i)
	    xkap=dust
	    xkaptot=dusto(i)

	  xg(i)=xkap*ppesop(i)
	  xgt(i)=xkaptot*ppesop(i)


	if (i.ge.2) then
	 xgg=.5*(xg(i-1)+xg(i))
	xggt=.5*(xgt(i-1)+xgt(i))

	ppesopm=.5*(ppesop(i-1)+ppesop(i))

	  del=1.e-4*(wl(i)-wl(i-1))

	  plancks=plancks+xgg*del
	  planckts=planckts+xggt*del
          xnorpp=xnorpp+ppesopm*del

	endif
	enddo

c	call fdate(d_and_t)
c	write(*,*)'TIME begin opa: ', d_and_t(12:19)

	osp=(plancks)/xnorpp
	ost=(planckts)/xnorpp

	xkaps=osp
	xchis=ost
1020	continue

	write(6,*)t,ts,xkaps,xchis

	return
1000	end


c.................................................
	function opasil(wl,alb)
	real MIMCUT
	character*24 d_and_t
        parameter(na=50)
        PARAMETER  ( MAXANG = 7, MOMDIM = 200 )
        COMPLEX  CREFIN, SFORW, SBACK, S1( MAXANG ), S2( MAXANG ),
     $          TFORW( 2 ), TBACK( 2 )

        LOGICAL  PERFCT,ANYANG, PRNT(2)
        dimension XMU(MAXANG)
	common/sizes/amin,amax,p
        dimension x(na),y(na),y2(na),y3(na),y4(na)


	data absil/0.0034/
	data rhosil/3.3/
	data pi/3.1416/
	
c	if (i.eq.1) then
c	call fdate(d_and_t)
c	write(6,*)'TIME begin opasil: ', d_and_t(12:19)
c	end if

c	write(*,*)'amin,amax,p',amin,amax,p
	PERFCT = .FALSE.
        ANYANG = .TRUE.
        MIMCUT = 1.E-6
        NUMANG=MAXANG

        DO 10 I = 1, NUMANG
         ANGLE = ( I - 1 )*180. / ( NUMANG - 1 )
         XMU( I ) = COS( PI / 180.*ANGLE )
   10 CONTINUE

        PRNT( 1 ) = .FALSE.
        PRNT( 2 ) = .FALSE.

c	call eps_draine(wl,xnr,xni)
	call eps_jager(wl,xnr,xni)
c	call eps_sil_jager2003(wl,xnr,xni)

	crefin=cmplx(xnr,-xni)
	deltaa=(alog10(amax)-alog10(amin))/(na-1)
	wlmic=wl

         do ia=1,na
         al=alog10(amin)+deltaa*(ia-1)
         amic=10.**al

         aa=amic*1.e-4
         xx=2.*pi*amic/wlmic
	
	if (xx.lt.1.e4) then
c	call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opasil: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 	call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opasil: ', d_and_t(12:19)
        else
        xx=1.e4
c			call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opasil: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 	call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opasil: ', d_and_t(12:19)
        endif


        x(ia)=alog(amic)

        y4(ia)=amic**4*amic**(-p)
        y(ia)=y4(ia)/aa*qext
        y2(ia)=y4(ia)/aa*qsca
        y3(ia)=y4(ia)/aa*(qext-qsca)
	enddo

        dnorm=trapz(x,y4,na)
        dopa=(
     *trapz(x,y,na)*3./4./rhosil*absil)/dnorm
        dscatt=(
     *trapz(x,y2,na)*3./4./rhosil*absil)/dnorm

	
	alb=dscatt/dopa
        opasil=dopa

	
	return
	end
c..................................................
	subroutine eps_draine(xl,xnr,xni)
c.	lee tabla de funciones dielectricas de draine 
c.	e interpola en lambda
c	implicit double precision (a-h,o-z)
	parameter (nwlmax=4000)
	character*40 nada
	common/retsil2/wl(nwlmax),rem(nwlmax),xim(nwlmax),nwl
	

	if (wl(1).eq.0.) then
	   open(41,file='draine_sil_suave.dat',
     *	status='old',form='formatted')

        read(41,'(a)')nada
        read(41,'(a)')nada
        read(41,'(a)')nada
        read(41,'(a)')nada
        read(41,'(a)')nada
        read(41,'(a)')nada
	read(41,'(a)')nada
        read(41,'(a)')nada
        read(41,'(a)')nada

        do il=1,nwlmax
        read(41,*,end=101)wl(il),eps1,eps2,rem(il),xim(il)
	rem(il)=rem(il)+1.
	
        if (wl(il).lt.1.e-2) go to 101
        enddo
101     nwl=il-1
	close(41)
	endif



c	if (xl.ge.wl(1)) then
c	extrapolamos logaritmicamente
	if (xl.ge.wl(20)) then
c.	tomamos los indices como constantes e iguales a los valores 
c.	en lambda=8.035E+02 (n=20)
	xnr=rem(20)
	xni=xim(20)
	elseif (xl.lt.wl(nwl)) then
c	extrapolamos logaritmicamente
	x=alog10(xl)
	x1=alog10(wl(nwl-1))
	x2=alog10(wl(nwl))
	y1=alog10(rem(nwl-1))
	y2=alog10(rem(nwl))
	xnr=y1+(y2-y1)*(x-x1)/(x2-x1)
	xnr=10.**xnr

	y1=alog10(xim(nwl-1))
	y2=alog10(xim(nwl))
	xni=y1+(y2-y1)*(x-x1)/(x2-x1)
	xni=10.**xni

	else
	do il=1,nwl
	if (xl.lt.wl(il)) then
	il1=il
	il2=il+1
	endif
	enddo
	
	x=alog10(xl)
	x1=alog10(wl(il1))
	x2=alog10(wl(il2))

	y1=rem(il1)
	y2=rem(il2)
	xnr=y1+(y2-y1)*(x-x1)/(x2-x1)

c	write(*,*)'real:'
c	write(*,*)'xl,x,x1,x2=',xl,x,x1,x2
c	write(*,*)'y1,y2,ep1=',y1,y2,ep1
c	write(*,*)'eps1(il1),eps1(il2)=',eps1(il1),eps1(il2)

	y1=(xim(il1))
	y2=(xim(il2))
	xni=y1+(y2-y1)*(x-x1)/(x2-x1)
c	write(*,*)'imag:'
c	write(*,*)'y1,y2,ep1=',y1,y2,ep1
c	write(*,*)'eps2(il1),eps2(il2)=',eps2(il1),eps2(il2)

	endif
c	ep1=10.**ep1
c	ep2=10.**ep2

	
c	xnimag=0.5*(sqrt(ep1**2+ep2**2)-ep1)
c          xnimag=sqrt(xnimag)
c        xnreal=0.5*(sqrt(ep1**2+ep2**2)+ep1)
c        xnreal=sqrt(xnreal)


cNOTA:
c.	xnr es la parte real del indice de refraccion
c.	xni es la parte imaginaria del indice de refraccion

c	ep1=xnr**2-xni**2
c        ep2=xnr*xni*2.

	return
	end
c................................................
	subroutine eps_jager(xl,xnr,xni)
c.	lee tabla de funciones dielectricas de draine 
c.	e interpola en lambda
c	implicit double precision (a-h,o-z)
	parameter (nwlmax=4000)
	character*40 nada
	character*120 filesil,filetroi,fileice1,fileice2
        common/dustfiles/filesil,filetroi,fileice1,fileice2

	common/retsil2/wl(nwlmax),rem(nwlmax),xim(nwlmax),nwl
	

c	if (wl(1).eq.0.) then
c	open(41,file='/home/ramiro/paola/COMUN/olmg50.lnk',
c     *	status='old',form='formatted')
c        do il=1,nwlmax
c        read(41,*,end=101)wl(il),rem(il),xim(il)
cc        if (wl(il).lt.1.e-2) go to 101
c        enddo
c101     nwl=il-1
c	close(41)
c	endif

	
       if (wl(1).eq.0.) then
c      open(41,file='olmg50.ric',
       open(41,file=filesil,
     * status='old',form='formatted')
	read(41,'(a)')nada
	read(41,'(a)')nada
	read(41,'(a)')nada
	read(41,'(a)')nada
	read(41,'(a)')nada
	read(41,'(a)')nada
        read(41,'(a)')nada
        read(41,'(a)')nada
        read(41,'(a)')nada
	
        do il=1,nwlmax
        read(41,*,end=101)wl(il),rem(il),xim(il)
	wl(il)=wl(il)*1.e4
c        if (wl(il).lt.1.e-2) go to 101
        enddo
101     nwl=il-1
       close(41)
       endif



	if (xl.le.wl(1)) then
	xnr=rem(1)
	xni=xim(1)
	elseif (xl.gt.wl(nwl)) then
c	extrapolamos logaritmicamente
	x=alog10(xl)
	x1=alog10(wl(nwl-1))
	x2=alog10(wl(nwl))
	y1=alog10(rem(nwl-1))
	y2=alog10(rem(nwl))
	xnr=y1+(y2-y1)*(x-x1)/(x2-x1)
	xnr=10.**xnr

	y1=alog10(xim(nwl-1))
	y2=alog10(xim(nwl))
	xni=y1+(y2-y1)*(x-x1)/(x2-x1)
	xni=10.**xni

	else
	do il=1,nwl
	if (xl.gt.wl(il)) then
	il1=il
	il2=il+1
	endif
	enddo
	
	x=alog10(xl)
	x1=alog10(wl(il1))
	x2=alog10(wl(il2))

	y1=rem(il1)
	y2=rem(il2)
	xnr=y1+(y2-y1)*(x-x1)/(x2-x1)

	y1=(xim(il1))
	y2=(xim(il2))
	xni=y1+(y2-y1)*(x-x1)/(x2-x1)
	endif
c	ep1=xnr**2-xni**2
c       ep2=xnr*xni*2.

	return
	end
c........................................

	function opaorg(wl,alb)
	real MIMCUT
	character*24 d_and_t
        parameter(na=50)
        PARAMETER  ( MAXANG = 7, MOMDIM = 200 )
        COMPLEX  CREFIN, SFORW, SBACK, S1( MAXANG ), S2( MAXANG ),
     $          TFORW( 2 ), TBACK( 2 )

        LOGICAL  PERFCT,ANYANG, PRNT(2)
	common/sizes/amin,amax,p
        dimension XMU(MAXANG)
        dimension x(na),y(na),y2(na),y3(na),y4(na)


c	data aborg/0.0041/
	data aborg/0.001/
	data rhoorg/1.5/
	data pi/3.1416/

	PERFCT = .FALSE.
        ANYANG = .TRUE.
        MIMCUT = 1.E-6
        NUMANG=MAXANG

        DO 10 I = 1, NUMANG
         ANGLE = ( I - 1 )*180. / ( NUMANG - 1 )
         XMU( I ) = COS( PI / 180.*ANGLE )
   10 CONTINUE

        PRNT( 1 ) = .FALSE.
        PRNT( 2 ) = .FALSE.

	call eps_org(wl,xnr,xni)

	crefin=cmplx(xnr,-xni)
	deltaa=(alog10(amax)-alog10(amin))/(na-1)
	wlmic=wl

         do ia=1,na
         al=alog10(amin)+deltaa*(ia-1)
         amic=10.**al

         aa=amic*1.e-4
         xx=2.*pi*amic/wlmic
	
	if (xx.lt.1.e4) then
c		call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opaorg: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opaorg: ', d_and_t(12:19)
        else
        xx=1.e4
c		call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opaorg: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opaorg: ', d_and_t(12:19)
        endif


        x(ia)=alog(amic)

        y4(ia)=amic**4*amic**(-p)
        y(ia)=y4(ia)/aa*qext
        y2(ia)=y4(ia)/aa*qsca
        y3(ia)=y4(ia)/aa*(qext-qsca)
	enddo

        dnorm=trapz(x,y4,na)
        dopa=(
     *trapz(x,y,na)*3./4./rhoorg*aborg)/dnorm
        dscatt=(
     *trapz(x,y2,na)*3./4./rhoorg*aborg)/dnorm

	
	alb=dscatt/dopa
        opaorg=dopa
	return
	end
c.........................................
	function opatroi(wl,alb)
	real MIMCUT
	character*24 d_and_t
        parameter(na=50)
        PARAMETER  ( MAXANG = 7, MOMDIM = 200 )
        COMPLEX  CREFIN, SFORW, SBACK, S1( MAXANG ), S2( MAXANG ),
     $          TFORW( 2 ), TBACK( 2 )

        LOGICAL  PERFCT,ANYANG, PRNT(2)
	common/sizes/amin,amax,p
        dimension XMU(MAXANG)
        dimension x(na),y(na),y2(na),y3(na),y4(na)

	data abtroi/7.68e-4/
	data rhotroi/4.83/
	data pi/3.1416/

	PERFCT = .FALSE.
        ANYANG = .TRUE.
        MIMCUT = 1.E-6
        NUMANG=MAXANG

        DO 10 I = 1, NUMANG
         ANGLE = ( I - 1 )*180. / ( NUMANG - 1 )
         XMU( I ) = COS( PI / 180.*ANGLE )
   10 CONTINUE

        PRNT( 1 ) = .FALSE.
        PRNT( 2 ) = .FALSE.

	call eps_troi(wl,xnr,xni)

	crefin=cmplx(xnr,-xni)
	deltaa=(alog10(amax)-alog10(amin))/(na-1)
	wlmic=wl

         do ia=1,na
         al=alog10(amin)+deltaa*(ia-1)
         amic=10.**al

         aa=amic*1.e-4
         xx=2.*pi*amic/wlmic
	
	if (xx.lt.1.e4) then
c	call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opatroi: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 	call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opatroi: ', d_and_t(12:19)
        else
        xx=1.e4
c			call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opatroi: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 	call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opatroi: ', d_and_t(12:19)
        endif


        x(ia)=alog(amic)

        y4(ia)=amic**4*amic**(-p)
        y(ia)=y4(ia)/aa*qext
        y2(ia)=y4(ia)/aa*qsca
        y3(ia)=y4(ia)/aa*(qext-qsca)
	enddo

        dnorm=trapz(x,y4,na)
        dopa=(
     *trapz(x,y,na)*3./4./rhotroi*abtroi)/dnorm
        dscatt=(
     *trapz(x,y2,na)*3./4./rhotroi*abtroi)/dnorm

	
	alb=dscatt/dopa
        opatroi=dopa
	return
	end
c................................................
	function opaice(wl,alb)
	real MIMCUT
	character*24 d_and_t
        parameter(na=50)
        PARAMETER  ( MAXANG = 7, MOMDIM = 200 )
        COMPLEX  CREFIN, SFORW, SBACK, S1( MAXANG ), S2( MAXANG ),
     $          TFORW( 2 ), TBACK( 2 )

        LOGICAL  PERFCT,ANYANG, PRNT(2)
	common/sizes/amin,amax,p
        dimension XMU(MAXANG)
        dimension x(na),y(na),y2(na),y3(na),y4(na)


c	data abice/0.0056/
c.	prueba
	data abice/0.00056/
	data rhoice/0.92/
	data pi/3.1416/

	PERFCT = .FALSE.
        ANYANG = .TRUE.
        MIMCUT = 1.E-6
        NUMANG=MAXANG

        DO 10 I = 1, NUMANG
         ANGLE = ( I - 1 )*180. / ( NUMANG - 1 )
         XMU( I ) = COS( PI / 180.*ANGLE )
   10 CONTINUE

        PRNT( 1 ) = .FALSE.
        PRNT( 2 ) = .FALSE.

	call eps_ice(wl,xnr,xni)
	crefin=cmplx(xnr,-xni)
	deltaa=(alog10(amax)-alog10(amin))/(na-1)
	wlmic=wl

         do ia=1,na
         al=alog10(amin)+deltaa*(ia-1)
         amic=10.**al

         aa=amic*1.e-4
         xx=2.*pi*amic/wlmic
	
	if (xx.lt.1.e4) then
c		call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opaice: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 	call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opaice: ', d_and_t(12:19)
        else
        xx=1.e4
c			call fdate(d_and_t)
c	write(*,*)'TIME calling mievo in opaice: ', d_and_t(12:19)
        call MIEV0( XX, CREFIN, PERFCT, MIMCUT,
     &  ANYANG, NUMANG, XMU,
     &         PRNT, QEXT, QSCA, GQSC,
     &         SFORW, SBACK, S1, S2, TFORW, TBACK,
     &         SPIKE )
c	 	call fdate(d_and_t)
c	write(*,*)'TIME end mievo in opaice: ', d_and_t(12:19)
        endif


        x(ia)=alog(amic)

        y4(ia)=amic**4*amic**(-p)
        y(ia)=y4(ia)/aa*qext
        y2(ia)=y4(ia)/aa*qsca
        y3(ia)=y4(ia)/aa*(qext-qsca)
	enddo

        dnorm=trapz(x,y4,na)
        dopa=(
     *trapz(x,y,na)*3./4./rhoice*abice)/dnorm
        dscatt=(
     *trapz(x,y2,na)*3./4./rhoice*abice)/dnorm

	
	alb=dscatt/dopa
        opaice=dopa
	return
	end
c...................................................
	subroutine eps_org(xl,xnr,xni)
c.	lee tabla de funciones dielectricas de draine 
c.	e interpola en lambda
c	implicit double precision (a-h,o-z)
	parameter (nwlmax=2000)
c	character*40 nada
c	complex*16 epsilon


	xxl=alog10(xl)

	if (xxl.lt.-1.0) then
	xni=10.**(-0.14)
	elseif (xxl.ge.-1.0.and.xxl.lt.0.) then
	xni=10.**(-1.83*xxl-1.97)
	elseif (xxl.ge.0..and.xxl.lt.1.4) then
	xni=10.**(1.135*xxl-1.97)
	elseif (xxl.ge.1.4.and.xxl.lt.2.57) then
	xni=10.**(-0.87*xxl+0.837)
	elseif (xxl.ge.2.57.and.xxl.lt.3.0) then
	xni=10.**(-2.58*xxl+5.2317)
	else
	xni=10.**(-0.17*xxl-1.9983)
	endif

	
	
	if (xxl.lt.-1.0) then
	xnr=1.75
	elseif (xxl.ge.-1.0.and.xxl.lt.1.18) then
	xnr=-0.08*xxl+1.67
	elseif (xxl.ge.1.18.and.xxl.lt.1.685) then
	xnr=1.42*xxl-0.1
	else
	xnr=2.2927
	endif
	

c        ep1=xnr**2-xni**2
c       ep2=xnr*xni*2.
	return
	end
c.....................................................

	subroutine eps_troi(xlam,xnrr,xnii)
c.	lee tabla de funciones dielectricas de draine 
c.	e interpola en lambda
c	implicit double precision (a-h,o-z)
	parameter (nwlmax=2000)
c	character*40 nada
	character*120 filesil,filetroi,fileice1,fileice2
        common/dustfiles/filesil,filetroi,fileice1,fileice2

	common/retienetroi/wl(nwlmax),xnr(nwlmax),xni(nwlmax),nwl
	

	if (wl(1).eq.0.) then
c	open(51,file='troilita.dat',
	open(51,file=filetroi,
     *	status='old',form='formatted')

        do il=1,nwlmax
        read(51,*,end=101)wl(il),xnr(il),xni(il)
	wl(il)=wl(il)*1.e4
c	eps1(il)=xnr**2-xni**2
c	eps2(il)=2.*xnr*xni
	
        enddo
101     nwl=il-1
	write(*,*)'nwl=',nwl,wl(1),wl(nwl)
	endif


	xl=alog10(xlam)

	if (xl.lt.alog10(wl(1))) then

		if (xl.le.-0.2) then
		xnrr=10.**(0.156)
		else
		xnrr=10.**(0.597796*xl+0.27556)
		endif

		if (xl.le.0.3) then
		xnii=10.**(0.52*xl+0.35)
		elseif (xl.gt.0.3.and.xl.le.0.57) then
		xnii=10.**(-6.27*xl+2.387)
		else
		xnii=10.**(-0.113429*xl-1.12225)
		endif

	return

	elseif (xl.gt.alog10(wl(nwl))) then

	if (xl.le.2.98214) then
	xnrr=10.**(0.917253)
	elseif (xl.gt.2.98214.and.xl.le.4.5) then
	xnrr=10.**(0.86*xl-1.64739)
	else
	xnrr=10.**(0.5*xl-0.02739)
	endif

	if (xl.le.4.5) then
	xnii=10.**(0.86*xl-1.64739)
	else
	xnii=10.**(0.5*xl-0.02739)
        endif

	

c	ep1=xnrr**2-xnii**2
c	ep2=xnrr*xnii*2.

	return

	else
	do il=1,nwl
	if (xlam.ge.wl(il)) then
	il1=il
	il2=il+1
	endif
	enddo
	
	x=xl
	x1=alog10(wl(il1))
	x2=alog10(wl(il2))
c	y1=alog10(eps1(il1))
c	y2=alog10(eps1(il2))

	y1=xnr(il1)
	y2=xnr(il2)
	ep1=y1+(y2-y1)*(x-x1)/(x2-x1)

c	write(*,*)'real:'
c	write(*,*)'xl,x,x1,x2=',xl,x,x1,x2
c	write(*,*)'y1,y2,ep1=',y1,y2,ep1
c	write(*,*)'eps1(il1),eps1(il2)=',eps1(il1),eps1(il2)

	y1=(xni(il1))
	y2=(xni(il2))
	ep2=y1+(y2-y1)*(x-x1)/(x2-x1)
c	write(*,*)'imag:'
c	write(*,*)'y1,y2,ep1=',y1,y2,ep1
c	write(*,*)'eps2(il1),eps2(il2)=',eps2(il1),eps2(il2)

	endif

	
	xnii=ep2
        xnrr=ep1
c	write(*,*)'xnii,xnrr=',xnii,xnrr
	return
	end
c...........................................................
	
	subroutine eps_ice(xl,xnreal,xnimag)
c.	lee tabla de funciones dielectricas de draine 
c.	e interpola en lambda
c	implicit double precision (a-h,o-z)
	parameter (nwlmax=2000,nt=4)
	character*40 nada
	character*120 filesil,filetroi,fileice1,fileice2
        common/dustfiles/filesil,filetroi,fileice1,fileice2

c	dimension xnr(nwlmax,nt),xni(nwlmax,nt)
	common/retieneice/wl(nwlmax),xnr(nwlmax,nt),
     *	xni(nwlmax,nt),nwl
	dimension xn1(nt),xn2(nt)
	

	if (wl(1).eq.0.) then
c	open(50,file='warren1984.dat',
	open(50,file=fileice1,
     *	status='old')

        do il=1,nwlmax
        read(50,*,end=101)wl(il),xnr(il,1),xni(il,1)
        do it=2,4
        xnr(il,it)=xnr(il,1)
        xni(il,it)=xni(il,1)
        enddo
        enddo
101     nwl1=il-1
        close(50)

c       open(50,file='warren1984.temp.dat',
        open(50,file=fileice2,
     *	status='old')

        read(50,'(a)')nada
        read(50,'(a)')nada
        read(50,'(a)')nada

        do il=nwl1+1,nwlmax
        read(50,*,end=102)wl(il),xnr(il,1),xni(il,1),
     *	xnr(il,2),xni(il,2),
     *  xnr(il,3),xni(il,3),xnr(il,4),xni(il,4)
        enddo
102     nwl=il-1

	endif

c	write(*,*)'xl,wl(1),wl(nwl)=',xl,wl(1),wl(nwl)
	do it=1,nt
	
	if (xl.le.wl(1)) then
	x=alog10(xl)
	x1=alog10(wl(1))
	x2=alog10(wl(2))

	y1=alog10(xnr(1,it))
	y2=alog10(xnr(2,it))
	xn1(it)=y1+(y2-y1)*(x-x1)/(x2-x1)
	xn1(it)=10.**xn1(it)

	y1=alog10(xni(1,it))
	y2=alog10(xni(2,it))
	xn2(it)=y1+(y2-y1)*(x-x1)/(x2-x1)
	xn2(it)=10.**xn2(it)

	elseif (xl.gt.wl(nwl)) then
	x=alog10(xl)
	x1=alog10(wl(nwl-1))
	x2=alog10(wl(nwl))

	y1=alog10(xnr(nwl-1,it))
	y2=alog10(xnr(nwl,it))
	xn1(it)=y1+(y2-y1)*(x-x1)/(x2-x1)
	xn1(it)=10.**xn1(it)

	y1=alog10(xni(nwl-1,it))
	y2=alog10(xni(nwl,it))
	xn2(it)=y1+(y2-y1)*(x-x1)/(x2-x1)
	xn2(it)=10.**xn2(it)

	else
	do il=1,nwl
	if (xl.ge.wl(il)) then
	il1=il
	il2=il+1
	endif
	enddo
	
	x=alog10(xl)
	x1=alog10(wl(il1))
	x2=alog10(wl(il2))

	

	y1=xnr(il1,it)
	y2=xnr(il2,it)
	xn1(it)=y1+(y2-y1)*(x-x1)/(x2-x1)

	y1=xni(il1,it)
	y2=xni(il2,it)
	xn2(it)=y1+(y2-y1)*(x-x1)/(x2-x1)

	endif
	enddo
	xnreal=xn1(4)
	xnimag=xn2(4)

	return
	end


c..........................................

	subroutine eps_sil_jager2003(xl,xnr,xni)
c	implicit double precision (a-h,o-z)
	parameter (nwlmax=4000)
	character*50 nada
	common/retsil3/wl(nwlmax),rem(nwlmax),xim(nwlmax),nwl
	dimension wll(nwlmax),reml(nwlmax),ximl(nwlmax)
	

	if (wl(1).eq.0.) then
	open(45,file='mg24.vnk',status='old',form='formatted')
        do il=1,nwlmax
        read(45,*,end=101)wll(il),reml(il),ximl(il)
        wll(il)=1./wll(il)*1.e4
	enddo
101     nwl=il-1
c.	invertimos
	do il=1,nwl
	k=nwl-il+1
	wl(il)=wll(k)
	rem(il)=reml(k)
	xim(il)=ximl(k)
c	write(40,*)il,wl(il),rem(il),xim(il)
	enddo
	close(45)
	write(*,*)'lei mg24.vnk'
	endif

	if (xl.gt.wl(nwl)) then
	xnr=rem(nwl)
	xni=xim(nwl)/(xl/wl(nwl))
	elseif (xl.le.wl(1)) then
	xnr=rem(1)
	xni=xim(1)
	else
	do il=1,nwl
	if (xl.gt.wl(il)) then
	il1=il
	il2=il+1
	endif
	enddo

	x=alog10(xl)
	x1=alog10(wl(il1))
	x2=alog10(wl(il2))

	y1=rem(il1)
	y2=rem(il2)
	xnr=y1+(y2-y1)*(x-x1)/(x2-x1)

	y1=(xim(il1))
	y2=(xim(il2))
	xni=y1+(y2-y1)*(x-x1)/(x2-x1)
c	write(43,*)xl,xnr,xni
	endif
	return
	end
c.................................................
