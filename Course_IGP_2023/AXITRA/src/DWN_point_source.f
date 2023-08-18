C ********** POINT SOURCE IN A LAYERED MEDIUM ***************

c  THE CODE USES THE DISCRETE WAVENUMBER METHOD (Bouchon,
c  BSSA, 71, 959-971, 1981) AND THE METHOD OF REFLECTIVITY
c  AND TRANSMISSIVITY MATRICES OF Kennett (Muller, Journal
c  of Geophysics, 158, 153-174, 1988).

c Parameters to be defined to dimension the arrays:
c - nlmax = maximum value of NL0
c - nrmax = maximum value of NR
c - ntmax = maximum value of ntime
c - nxmax = maximum value of NX
c - nzmax = maximum value of NZ
c ici NX=NZ=1 : POINT SOURCE !
c - mmax = maximum value of M+1 (use mmax=2000)
      parameter (nlmax=10)
!Hugo increased number of receivers nrmax
      parameter (nrmax=525)
      parameter (ntmax=1024)
      parameter (nxmax=1)
      parameter (nzmax=1)
      parameter (mmax=4000)

      parameter (NLMX=NLMAX+2+NZMAX)

      DIMENSION TH(NLMX),AL0(NLMX),BE0(NLMX),DENS(NLMX),
     $  QP(NLMX),QS(NLMX),
     $  R0(NRMAX),AZ(NRMAX),XR0(NRMAX),YR0(NRMAX),
     $  LAY(NZMAX),
     $  R(NXMAX,NZMAX,NRMAX),
     $  S2T(NXMAX,NZMAX,NRMAX),C2T(NXMAX,NZMAX,NRMAX),
     $  ST(NXMAX,NZMAX,NRMAX),CT(NXMAX,NZMAX,NRMAX),
     $  AJ0(NXMAX,NZMAX,NRMAX),AJ1(NXMAX,NZMAX,NRMAX),
     $  yyy(ntmax),y(ntmax+ntmax+3),y2(ntmax+ntmax+3),
     $  sy(ntmax,nrmax,3,2),fctsou(ntmax,2)

      COMPLEX*16 ALPHA(NLMX),BETA(NLMX),CTH(NLMX),EMU(NLMX),
     $  EMU2(NLMX),CMU(NLMX),WA2(NLMX),WB2(NLMX),
     $  WZA(NLMX),WZB(NLMX),AQP(NLMX),AQS(NLMX),
     $  U(NTMAX/2,NRMAX,3),
     $  DU1(NZMAX,5),DU2(NZMAX,5),DU3(NZMAX,5),
     $  AMPSV(NZMAX),AMSH(NZMAX),
     $  EV(NXMAX,NZMAX),TEV(NXMAX,NZMAX),
     $  bs(2),bl(2),su(3,5),sd(3,5),
     $  woa(nlmx),wob(nlmx),
     $  rd(2,2,nlmx),ru(2,2,nlmx),td(2,2,nlmx),tu(2,2,nlmx),
     $  mb(2,2,nlmx),mt(2,2,nlmx),nb(2,2,nlmx),nt(2,2,nlmx),
     $  quu(2,2,nlmx),qud(2,2,nlmx),
     $  tup(2,2,nlmx,nlmx),
     $  g(2,2,nlmx),ee(2,2,nlmx),c(2,2),d(2,2),e(2,2),
     $  rdsh(nlmx),rush(nlmx),tdsh(nlmx),tush(nlmx),
     $  mbsh(nlmx),mtsh(nlmx),nbsh(nlmx),ntsh(nlmx),
     $  quush(nlmx),qudsh(nlmx),
     $  tupsh(nlmx,nlmx),
     $  gsh(nlmx),eesh(nlmx),
     $  source(ntmax/2,2)

      COMPLEX  EJ1(NXMAX,NRMAX,NZMAX),EJ2(NXMAX,NRMAX,NZMAX),
     $  EJ3(NXMAX,NRMAX,NZMAX)

      COMPLEX*16 AI,AIPI,OMEGA,DOM,XLNF,AIK,AIK2,
     $  C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C14,C15,C16,C17,
     $  C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C35,C36,
     $  cu,cu2,cc,d1d,d2d,d1u,d2u,dd,bssh,blsh

      REAL*8 ARG,BEJ0,BEJ1,BEY0,BEY1

      real*8 ak,ak2,pi,pi2,pil,freq,rw,aw

      complex*8 q1,q2
      
      character*3 carstat
      character*30 cnord,cest,cvert
      character*30 dir
      integer ndir

      logical rt

      AI=(0.D0,1.D0)
      PI=3.141592653589793D0
      PI2=PI+PI
      EPS=1.E-3

c DEFINE THE INPUT AND OUTPUT FILES:

c The input data file:
      open(10,file="DWN_point_source.in",form="formatted")
      open(18,file="glissement.out")
      open(19,file="vitesse.out") 
      open(22,file="nomstat.in")



c ***********************************************************************
c *                          INPUT                                      *
c ***********************************************************************

c - NL0 : number of layers.
c - TH : layer thickness
c - AL0, BE0 : P and S wave velocities
c - DENS : density
c - QP, QS : Q of P and S waves
c - ZH  : hypocentral depth
c   The location of the receivers is defined in a cartesian
c   coordinate system (x,y,z), where (x,y) lies on the earth surface and is
c   centered at the epicenter. x is North, y is East, and z is positive downward. 
c - STRIKE : strike of the fault measured clockwise from North. 
c   DIP : dip of the fault. 
c   RAKE : direction of slip of the hanging wall relatively to the foot wall
c   It is measured counterclockwise from the strike direction.
c   If faulting is right-lateral (like the San Andreas or the North Anatolian), rake=0.
c - NR : number of receivers
c   for each receiver distance and azimuth in degree
c - NTIME : number of points of each seismogram
c - TL : length of time window
c - RMOMENT = seismic moment (N/m)
c - TSOURCE : rise time of the slip (s)
c   LOOK at the various possibilities for source time functions  (STF)shapes in the program:
c          -  integral of 'Bouchon like' smooth ramp
c          -  isocele triangular pulse
c          -  gaussian STF : in this case tsource is the standard error of the gaussian
c - TSHIFT : shift of the source time fucntion
c           - for 'bouchon like' smooth ramp, shift of the beginning of the STF
c           - for isocele triangular pulse, shift of the beginning of the STF
c           - for gaussian STF, shift of the center of the gaussian (use for example 2*TSOURCE)

c - RT : determines if the horizontal seismograms calculated are (radial, tangential) with
c   respect to the hypocenter in which case RT=.TRUE. , or (north, east) in which case 
c   RT=.FALSE.

c  All lengths (that is: distances, thicknesses, wave velocities) should 
c  be expressed in km.
c  if moment is expressed in N.m , seismograms are in m

      READ(10,*) dir
      READ(10,*) NL0
      DO 1 L=1,NL0
    1 READ(10,*) TH(L),AL0(L),BE0(L),DENS(L),QP(L),QS(L)
      READ(10,*) ZH
      READ(10,*) STRIKE,DIP,RAKE
      READ(10,*) NR
      do ir=1,NR
      READ(10,*) r0(ir),az(ir)
      enddo
      READ(10,*) NTIME,TL
      READ(10,*) RMOMENT
      READ(10,*) TSOURCE
      READ(10,*) TSHIFT
      RT=.FALSE.
 

c  ***************************** END OF INPUT *********************

      WRITE(6,*) NL0
      DO 11 L=1,NL0
   11 WRITE(6,*) TH(L),AL0(L),BE0(L),DENS(L),QP(L),QS(L)
      WRITE(6,*) ZH
      WRITE(6,*) STRIKE,DIP,RAKE

       WRITE(*,*) STRIKE,DIP,RAKE
                  
      DO 16 IR=1,NR
   16 WRITE(6,*) R0(IR),AZ(IR)
      WRITE(6,*) NTIME,TL
      WRITE(6,*) RMOMENT
      WRITE(6,*) TSOURCE
      WRITE(6,*) TSHIFT

cc NX,NZ,ZFO, ZF1 FAULTL1,FAULTL2,RMOMENT are defined for adaptation to the extended source
cc add 180 degre to rake for compatibility with other programs
      NX=1
      NZ=1
      ZF0=ZH-0.7
      ZF1=ZH+0.7
      FAULTL1=0.7
      FAULTL2=0.7
c     PRUEBA
      RMOMENT=RMOMENT/1.e15
      RAKE=RAKE + 180.
      TSOURCE=TSOURCE/4.0

c - M : maximum truncation index of the series.
c - XL : periodicity length.
      M=MMAX-1
      XL=AL0(NL0)*TL
      DO 1112 IR=1,NR
      XX=AL0(NL0)*TL+R0(IR)
      IF(XX.GT.XL) XL=XX
 1112 CONTINUE
      

      DO 15 IR=1,NR
      XR0(IR)=R0(IR)*COS(AZ(IR)*PI/180.)
   15 YR0(IR)=R0(IR)*SIN(AZ(IR)*PI/180.)
      CS=COS(STRIKE*PI/180.)
      SS=SIN(STRIKE*PI/180.)
      CDI=COS(DIP*PI/180.)
      SDI=SIN(DIP*PI/180.)
      CR=COS(RAKE*PI/180.)
      SR=SIN(RAKE*PI/180.)
      AS1=CR*CS+SR*CDI*SS
      AS2=CR*SS-SR*CDI*CS
      AS3=-SR*SDI
      WRITE(6,*) 'SLIP VECTOR', AS1, AS2, AS3
      AN1=-SDI*SS
      AN2=SDI*CS
      AN3=-CDI
      WRITE(6,*) 'NORMAL VECTOR', AN1, AN2, AN3
      CM11=-2.*AS1*AN1
      CM22=-2.*AS2*AN2
      CM33=-2.*AS3*AN3
      CM12=-(AS1*AN2+AS2*AN1)
      CM13=-(AS1*AN3+AS3*AN1)
      CM23=-(AS2*AN3+AS3*AN2)

      NL=NL0+1
      DO 110 L=NL,2,-1
      L1=L-1
      TH(L)=TH(L1)
      AL0(L)=AL0(L1)
      BE0(L)=BE0(L1)
      DENS(L)=DENS(L1)
      QP(L)=QP(L1)
  110 QS(L)=QS(L1)
      TH(NL)=0.

c We add a upper half-space which has the elastic properties of air:
      TH(1)=0.
      AL0(1)=0.340
      BE0(1)=1e-6
      DENS(1)=1.3E-3
      QP(1)=1.E6
      QS(1)=1.E6
c **

      DZ=(ZF1-ZF0)/DFLOAT(NZ)
      Z0=ZF0+DZ/2.
      Z1=ZF1-DZ/2.

      DL=(FAULTL1+FAULTL2)/DFLOAT(NX)
      DX=DL*CS
      DY=DL*SS
      DL2=DL/2.
      X0=-(FAULTL1-DL2)*CS
      X1=(FAULTL2-DL2)*CS
      Y0=-(FAULTL1-DL2)*SS
      Y1=(FAULTL2-DL2)*SS

      do 112 iz=1,nz
      zz=z0+dz*float(iz-1)
      zt=0.
      do 114 l=2,nl
      zb=zt+th(l)
      zt=zb
      if(zz.le.zb) go to 115
  114 continue
      nl=nl+1
      l1=nl-1
      al0(nl)=al0(l1)
      be0(nl)=be0(l1)
      dens(nl)=dens(l1)
      qp(nl)=qp(l1)
      qs(nl)=qs(l1)
      th(nl)=0.
      th(l1)=zz-zb
      lay(iz)=nl
      go to 112
      
  115 nl=nl+1
      do 117 ll=nl,l+1,-1
      l1=ll-1
      th(ll)=th(l1)
      al0(ll)=al0(l1)
      be0(ll)=be0(l1)
      dens(ll)=dens(l1)
      qp(ll)=qp(l1)
  117 qs(ll)=qs(l1)
      th(l+1)=zb-zz
      th(l)=th(l)-th(l+1)
      lay(iz)=l+1
  112 continue
      nl1=nl-1


      Q=TL
      PIL=PI2/XL                    
      AIPI=AI*PIL                   
      DT=TL/DFLOAT(NTIME)              
      DFREQ=1.D0/TL                 
      NFREQ=NTIME/2                    
      AW=-PI/Q 

c The introduction of anelastic attenuation produces some
c dispersion of the velocities, although very light. freq0
c is the frequency of reference for the velocities.
c For more... see the discussion in Aki and Richards on Q.
      DO 10 L=1,NL
      AQP(L)=(1.+AI/(QP(L)+QP(L)))/(1.+.25/QP(L)**2)*AL0(L)
      AQS(L)=(1.+AI/(QS(L)+QS(L)))/(1.+.25/QS(L)**2)*BE0(L)
   10 CTH(L)=-AI*TH(L)

      DO 14 IR=1,NR    
      DO 14 IX=1,NX
      DO 14 IZ=1,NZ
      ZZ=ZH-(Z0+DZ*FLOAT(IZ-1))
      A1=XR0(IR)-(X0+DX*FLOAT(IX-1))-ZZ*CDI/SDI*SS
      A2=YR0(IR)-(Y0+DY*FLOAT(IX-1))+ZZ*CDI/SDI*CS
      E0=SQRT(A1**2+A2**2)      
      IF(E0.LT.DL2) E0=DL2      
      R(IX,IZ,IR)=E0 
      TETA=ATAN2(A2,A1)*180./PI
      ST(IX,IZ,IR)=SIN(TETA*PI/180.)
      CT(IX,IZ,IR)=COS(TETA*PI/180.)
      S2T(IX,IZ,IR)=2.*ST(IX,IZ,IR)*CT(IX,IZ,IR)
      C2T(IX,IZ,IR)=CT(IX,IZ,IR)**2-ST(IX,IZ,IR)**2     
   14 A0=A0-DX 

      DOM=-AI*DFREQ*PI2             
      DO 72 IZ=1,NZ      
      DO 72 IX=1,NX
      XX=FAULTL1-DL2-DL*FLOAT(IX-1)
      ZZ=Z0+DZ*FLOAT(IZ-1)-ZH
c      E0=SQRT(XX**2+(ZZ/SDI)**2)/VR
      E0=0.
c - pas de propagation de la rupture: point source
      EV(IX,IZ)=EXP(E0*AW)
c - pas de glissement : remplace par le moment
      EV(IX,IZ)=EV(IX,IZ)*RMOMENT
      C0=E0*DOM          
      TEV(IX,IZ)=EXP(C0)
   72 CONTINUE

      FREQ=0.D0

c** Begin loop on frequencies


      DO 100 IF=1,NFREQ             
      RW=PI2*FREQ                   
      OMEGA=CMPLX(RW,AW)           

c  source is the Fourier transform of a smooth ramp function of
c  rise time equal to tsource:
c      tstart=-tsource-TSHIFT
c      c1=exp(omega*pi*tsource/4.)
c      source(if,1)=-ai*pi*tsource/2./(c1-1./c1)
c     &*exp(ai*omega*tstart)
c      source(if,2)=ai*omega*source(if,1)

c     fourier tranform of a ramp
c      source(if,1)=exp(-ai*omega*TSHIFT)*
c     &4*ai*(1-exp(-ai*omega*tsource/2.))**2/
c     &(omega**3*tsource**2)
c      source(if,2)=ai*omega*source(if,1)

c     fourier transform of a gaussian
      source(if,2)=
     &exp(-ai*omega*TSHIFT)*exp(-(omega*tsource/2.)**2)
      source(if,1)=source(if,2)/(ai*omega)

      FREQ0=2.5
      ZOM=SQRT(RW**2+AW**2)/PI2
      IF(IF.EQ.1) PPHI=-PI/2.
      IF(IF.NE.1) PPHI=ATAN(AW/RW)
      XLNF=(AI*PPHI+ALOG(ZOM)-ALOG(FREQ0))/PI
      DO 82 L=1,NL                   
      ALPHA(L)=AQP(L)/(1.-XLNF/QP(L))
      IF(FREQ.EQ.0.) ALPHA(L)=AL0(L)
      BETA(L)=AQS(L)/(1.-XLNF/QS(L))
      IF(FREQ.EQ.0.) BETA(L)=BE0(L)
      EMU(L)=BETA(L)**2*DENS(L)
      EMU2(L)=EMU(L)+EMU(L)
      CMU(L)=AI*EMU(L)
      WA2(L)=(OMEGA/ALPHA(L))**2    
   82 WB2(L)=(OMEGA/BETA(L))**2     

      DO 86 IZ=1,NZ
      L=LAY(IZ)
      A0=DFREQ/((XL+XL)*DENS(L))
      AMPSV(IZ)=A0/(OMEGA**2)
   86 AMSH(IZ)=A0/(BETA(L)**2)

      DO 40 J=1,3
      DO 40 IR=1,NR                 
   40 U(IF,IR,J)=0.D0

      DO 578 IZ=1,NZ
      DO 578 IX=1,NX
      DO 578 IR=1,NR
      EJ1(IX,IR,IZ)=0.
      EJ2(IX,IR,IZ)=0.
  578 EJ3(IX,IR,IZ)=0.


c For each frequency, one calculates the contribution
c of each wavenumber k. The following do loop on wavenumbers
c represents the discrete sum on wavenumbers which replaces
c the wavenumber integral.

      nzk=nz
      nlk=nl
      nlk1=nlk-1
      do 50 k=1,m
      ak=pil*dfloat(k)
      ak2=ak**2
      aik=ai*ak
      aik2=ai*ak2

      DO 556 IR=1,NR
      DO 556 IX=1,NX
      DO 556 IZ=1,NZ
      ARG=AK*R(IX,IZ,IR)
      CALL FF01AD(BEJ0,BEY0,ARG,0)
      CALL FF02AD(BEJ1,BEY1,ARG,0)
      AJ0(IX,IZ,IR)=BEJ0
  556 AJ1(IX,IZ,IR)=BEJ1

      do 30 l=1,nlk
      c1=wa2(l)-ak2
      wza(l)=sqrt(c1)
      q1=wza(l)
      if(aimag(q1).gt.0.) wza(l)=-wza(l)
      c2=wb2(l)-ak2
      wzb(l)=sqrt(c2)
      q1=wzb(l)
      if(aimag(q1).gt.0.) wzb(l)=-wzb(l)
      woa(l)=wza(l)/omega
      wob(l)=wzb(l)/omega
   30 continue

c --------------------------------------------------------------------------
c
c Calculation of the reflection and transmission matrices for each value of k 
c
c --------------------------------------------------------------------------
      cu=ak/omega
      cu2=cu**2
      do 301 l=1,nlk1
      l1=l+1
      cc=emu2(l)-emu2(l1)
      c1=cc*cu2
      c2=c1-dens(l)
      c3=c1+dens(l1)
      c4=c1-dens(l)+dens(l1)
      c5=c2*c2
      c6=c3*c3
      c7=c4*c4*cu2
      a1=dens(l)*dens(l1)
      c8=woa(l)*wob(l)
      c9=woa(l)*wob(l1)
      c10=woa(l1)*wob(l)
      c11=woa(l1)*wob(l1)
      c14=a1*c9
      c15=a1*c10
      c16=cc*c1*c8*c11
      c17=c5*c11
      c18=c6*c8
      d1d=c7+c17+c15
      d2d=c16+c18+c14
      d1u=c7+c18+c14
      d2u=c16+c17+c15
      c19=c3*wob(l)-c2*wob(l1)
      c20=c3*woa(l)-c2*woa(l1)
      dd=d1d+d2d
      rd(1,1,l1)=(d2d-d1d)/dd
      ru(1,1,l1)=(d2u-d1u)/dd
      c21=(cu+cu)*woa(l)
      c22=(cu+cu)*wob(l)
      c23=(cu+cu)*woa(l1)
      c24=(cu+cu)*wob(l1)
      c25=(c4*c3+cc*c2*c11)/dd
      rd(2,1,l1)=-c21*c25
      c35=(c4*c2+cc*c3*c8)/dd
      ru(2,1,l1)=c23*c35
      c26=dens(l)/dd
      td(1,1,l1)=(c26+c26)*woa(l)*c19
      td(2,1,l1)=-c26*c21*(c4+cc*c10)
      c27=(a1+a1)*(c10-c9)
      rd(2,2,l1)=(d2d-d1d+c27)/dd
      rd(1,2,l1)=c22*c25
      td(2,2,l1)=(c26+c26)*wob(l)*c20
      td(1,2,l1)=c26*c22*(c4+cc*c9)
      c36=dens(l1)/dd
      tu(1,1,l1)=(c36+c36)*woa(l1)*c19
      tu(2,1,l1)=-c36*c23*(c4+cc*c9)
      ru(2,2,l1)=(d2u-d1u-c27)/dd
      ru(1,2,l1)=-c24*c35
      tu(2,2,l1)=(c36+c36)*wob(l1)*c20
      tu(1,2,l1)=c36*c24*(c4+cc*c10)
  301 continue

      do 304 i=1,2
      do 304 j=1,2
      mt(i,j,nlk)=0.d0
      mb(i,j,nlk1)=rd(i,j,nlk)
      nb(i,j,1)=0.d0
      nt(i,j,2)=ru(i,j,2)
  304 g(i,j,1)=tu(i,j,2)
      do 303 l=2,nlk1
      c1=cth(l)*wza(l)
      c2=cth(l)*wzb(l)
      ee(1,1,l)=exp(c1)
      ee(2,2,l)=exp(c2)
      ee(1,2,l)=0.d0
  303 ee(2,1,l)=0.d0
      do 306 l=nlk1,2,-1
      l1=l-1
      c1=ee(1,1,l)
      c2=ee(2,2,l)
      c3=c1*c2
      mt(1,1,l)=mb(1,1,l)*c1*c1
      mt(1,2,l)=mb(1,2,l)*c3
      mt(2,1,l)=mb(2,1,l)*c3
      mt(2,2,l)=mb(2,2,l)*c2*c2
      do 308 i=1,2
      do 308 j=1,2
      c(i,j)=0.d0
      do 308 ij=1,2
  308 c(i,j)=c(i,j)+mt(i,ij,l)*ru(ij,j,l)
      e(1,1)=1.d0-c(1,1)
      e(2,2)=1.d0-c(2,2)
      e(1,2)=-c(1,2)
      e(2,1)=-c(2,1)
      call inv2(e)
      do 310 i=1,2
      do 310 j=1,2
      c(i,j)=0.d0
      do 310 ij=1,2
  310 c(i,j)=c(i,j)+tu(i,ij,l)*e(ij,j)
      do 312 i=1,2
      do 312 j=1,2
      e(i,j)=0.d0
      do 312 ij=1,2
  312 e(i,j)=e(i,j)+c(i,ij)*mt(ij,j,l)
      do 314 i=1,2
      do 314 j=1,2
      mb(i,j,l1)=rd(i,j,l)
      do 314 ij=1,2
  314 mb(i,j,l1)=mb(i,j,l1)+e(i,ij)*td(ij,j,l)
  306 continue
      do 316 l=2,nlk1
      l1=l+1
      c1=ee(1,1,l)
      c2=ee(2,2,l)
      c3=c1*c2
      nb(1,1,l)=nt(1,1,l)*c1*c1
      nb(1,2,l)=nt(1,2,l)*c3
      nb(2,1,l)=nt(2,1,l)*c3
      nb(2,2,l)=nt(2,2,l)*c2*c2
      do 318 i=1,2
      do 318 j=1,2
      c(i,j)=0.d0
      do 318 ij=1,2
  318 c(i,j)=c(i,j)+nb(i,ij,l)*rd(ij,j,l1)
      e(1,1)=1.d0-c(1,1)
      e(2,2)=1.d0-c(2,2)
      e(1,2)=-c(1,2)
      e(2,1)=-c(2,1)
      call inv2(e)
      do 320 i=1,2
      do 320 j=1,2
      c(i,j)=0.d0
      do 320 ij=1,2
  320 c(i,j)=c(i,j)+td(i,ij,l1)*e(ij,j)
      do 322 i=1,2
      do 322 j=1,2
      e(i,j)=0.d0
      do 322 ij=1,2
  322 e(i,j)=e(i,j)+c(i,ij)*nb(ij,j,l)
      do 324 i=1,2
      do 324 j=1,2
      nt(i,j,l1)=ru(i,j,l1)
      do 324 ij=1,2
  324 nt(i,j,l1)=nt(i,j,l1)+e(i,ij)*tu(ij,j,l1)
  316 continue
      do 350 l=2,nlk1
      do 451 i=1,2
      do 451 j=1,2
      c(i,j)=0.d0
      do 451 ij=1,2
  451 c(i,j)=c(i,j)+mt(i,ij,l)*nt(ij,j,l)
      e(1,1)=1.d0-c(1,1)
      e(1,2)=-c(1,2)
      e(2,1)=-c(2,1)
      e(2,2)=1.d0-c(2,2)
      call inv2(e)
      c(1,1)=ee(1,1,l)*mb(1,1,l)
      c(1,2)=ee(1,1,l)*mb(1,2,l)
      c(2,1)=ee(2,2,l)*mb(2,1,l)
      c(2,2)=ee(2,2,l)*mb(2,2,l)
      do 452 i=1,2
      do 452 j=1,2
      quu(i,j,l)=e(i,j)
      qud(i,j,l)=0.d0
      do 452 ij=1,2
  452 qud(i,j,l)=qud(i,j,l)+e(i,ij)*c(ij,j)
  350 continue
      do 370 i=1,2
      do 370 j=1,2
      quu(i,j,1)=0.d0
      quu(i,j,nlk)=0.d0
      qud(i,j,1)=mb(i,j,1)
  370 qud(i,j,nlk)=0.d0
      do 371 i=1,2
  371 quu(i,i,nlk)=1.d0
      do 380 l=1,nlk1
      l1=l+1
      do 381 i=1,2
      do 381 j=1,2
      c(i,j)=0.d0
      d(i,j)=0.d0
      do 381 ij=1,2
      c(i,j)=c(i,j)-rd(i,ij,l1)*nb(ij,j,l)
  381 d(i,j)=d(i,j)-ru(i,ij,l1)*mt(ij,j,l1)
      c(1,1)=1.d0+c(1,1)
      c(2,2)=1.d0+c(2,2)
      d(1,1)=1.d0+d(1,1)
      d(2,2)=1.d0+d(2,2)
      call inv2(c)
      call inv2(d)
      do 382 i=1,2
      do 382 j=1,2
      g(i,j,l)=0.d0
      do 382 ij=1,2
  382 g(i,j,l)=g(i,j,l)+c(i,ij)*tu(ij,j,l1)
  380 continue

      lr=1
      do 512 iz=1,nz
      l=lay(iz)
      l1=l-1
      lu=lr
      do 5100 i=1,2
      do 5100 j=1,2
 5100 tup(i,j,l,lu)=0.d0
      tup(1,1,l,lu)=1.0d0
      tup(2,2,l,lu)=1.0d0
      do 512 ll=lu,l1
      if(ll.eq.lu) go to 514
      tup(1,1,l,lu)=tup(1,1,l,lu)*ee(1,1,ll)
      tup(1,2,l,lu)=tup(1,2,l,lu)*ee(2,2,ll)
      tup(2,1,l,lu)=tup(2,1,l,lu)*ee(1,1,ll)
      tup(2,2,l,lu)=tup(2,2,l,lu)*ee(2,2,ll)
  514 do 511 i=1,2
      do 511 j=1,2
      d(i,j)=0.d0
      do 511 ij=1,2
  511 d(i,j)=d(i,j)+tup(i,ij,l,lu)*g(ij,j,ll)
      do 515 i=1,2
      do 515 j=1,2
  515 tup(i,j,l,lu)=d(i,j)
  512 continue
c ------------------------------------------------------------------------
c
c
c ------------------------------------------------------------------------

      do 1024 iz=1,nzk
      l=lay(iz)
      c1=cm12*ampsv(iz)*aik2
      sd(1,1)=-c1*ak/wza(l)
      su(1,1)=sd(1,1)
      sd(2,1)=c1
      su(2,1)=-sd(2,1)

      c2=ampsv(iz)*(ak2+ak2)
      c3=ampsv(iz)*ak*(ak2/wzb(l)-wzb(l))
      sd(1,2)=cm13*c2
      su(1,2)=-sd(1,2)
      sd(2,2)=cm13*c3
      su(2,2)=sd(2,2)

      sd(1,3)=cm23*c2
      su(1,3)=-sd(1,3)
      sd(2,3)=cm23*c3
      su(2,3)=sd(2,3)

      c4=ampsv(iz)*aik2/2.d0*(cm11-cm22)
      sd(1,4)=-c4*ak/wza(l)
      su(1,4)=sd(1,4)
      sd(2,4)=c4
      su(2,4)=-sd(2,4)

      sd(1,5)=ampsv(iz)*aik*(ak2/wza(l)*(cm11+cm22)/2.d0+wza(l)*cm33)
      su(1,5)=sd(1,5)
      sd(2,5)=ampsv(iz)*aik2*(cm33-(cm11+cm22)/2.d0)
      su(2,5)=-sd(2,5)

      do 1024 im=1,5
      do 801 ip=1,2
      bs(ip)=0.d0
      do 801 jp=1,2
  801 bs(ip)=bs(ip)+quu(ip,jp,l)*su(jp,im)+qud(ip,jp,l)*
     $  sd(jp,im)*ee(jp,jp,l)
      do 806 ip=1,2
      bl(ip)=0.d0
      do 806 jp=1,2
  806 bl(ip)=bl(ip)+tup(ip,jp,l,1)*bs(jp)

      du1(iz,im)=bl(1)+wzb(1)/ak*bl(2)
      du2(iz,im)=ai*(wza(1)*bl(1)-ak*bl(2))

 1024 continue

c ** sh case:

      do 1301 l=1,nlk1
      l1=l+1
      c1=emu(l)*wob(l)
      c2=emu(l1)*wob(l1)
      rdsh(l1)=(c1-c2)/(c1+c2)
      tdsh(l1)=(c1+c1)/(c1+c2)
      rush(l1)=-rdsh(l1)
 1301 tush(l1)=(c2+c2)/(c1+c2)
      mtsh(nlk)=0.d0
      mbsh(nlk1)=rdsh(nlk)
      nbsh(1)=0.d0
      ntsh(2)=rush(2)
      gsh(1)=tush(2)
      do 1303 l=2,nlk1
      c2=cth(l)*wzb(l)
 1303 eesh(l)=exp(c2)
      do 1306 l=nlk1,2,-1
      l1=l-1
      mtsh(l)=mbsh(l)*eesh(l)**2
 1306 mbsh(l1)=rdsh(l)+tdsh(l)*tush(l)*mtsh(l)/(1.d0-rush(l)*mtsh(l))
      do 1316 l=2,nlk1
      l1=l+1
      nbsh(l)=ntsh(l)*eesh(l)**2
 1316 ntsh(l1)=rush(l1)+tdsh(l1)*tush(l1)*nbsh(l)/
     $  (1.d0-rdsh(l1)*nbsh(l))
      do 1350 l=2,nlk1
      quush(l)=1.d0/(1.d0-mtsh(l)*ntsh(l))
      qudsh(l)=eesh(l)*mbsh(l)/(1.d0-mtsh(l)*ntsh(l))
 1350 continue
      quush(1)=0.d0
      qudsh(1)=mbsh(1)
      qudsh(nlk)=0.d0
      quush(nlk)=1.d0
      do 1380 l=2,nlk1
      l1=l+1
 1380 gsh(l)=tush(l1)/(1.d0-rdsh(l1)*nbsh(l))

      lr=1
      do 1512 iz=1,nzk
      l=lay(iz)
      l1=l-1
      lu=lr
      tupsh(l,lu)=1.d0
      do 1512 ll=lu,l1
      if(ll.eq.lu) go to 1514
      tupsh(l,lu)=tupsh(l,lu)*eesh(ll)
 1514 tupsh(l,lu)=tupsh(l,lu)*gsh(ll)
 1512 continue

      do 1025 iz=1,nzk
      l=lay(iz)
      sd(3,1)=cm12*amsh(iz)*aik/wzb(l)
      su(3,1)=sd(3,1)
      sd(3,2)=cm13*amsh(iz)
      su(3,2)=-sd(3,2)
      sd(3,3)=-cm23*amsh(iz)
      su(3,3)=-sd(3,3)
      sd(3,4)=amsh(iz)*aik/(2.d0*wzb(l))*(cm22-cm11)
      su(3,4)=sd(3,4)
      do 1025 im=1,4
      bssh=quush(l)*su(3,im)+qudsh(l)*sd(3,im)*eesh(l)
      blsh=tupsh(l,1)*bssh

      du3(iz,im)=blsh

 1025 continue
      
      DO 580 IZ=1,NZK
      DO 580 IX=1,NX
      DO 580 IR=1,NR

      AJ1R=AJ1(IX,IZ,IR)/R(IX,IZ,IR)
      AJKR=AK*AJ0(IX,IZ,IR)-AJ1R
      AJ2=(AJ1R+AJ1R)/AK-AJ0(IX,IZ,IR)
      AJ2R=(AJ2+AJ2)/R(IX,IZ,IR)
      AJ1K=AK*AJ1(IX,IZ,IR)
      DAJ2=AJ1K-AJ2R

      EJ1(IX,IR,IZ)=EJ1(IX,IR,IZ)+
     $  S2T(IX,IZ,IR)*(DAJ2*DU1(IZ,1)-AJ2R*DU3(IZ,1))
     $  +CT(IX,IZ,IR)*(AJKR*DU1(IZ,2)+AJ1R*DU3(IZ,2))
     $  +ST(IX,IZ,IR)*(AJKR*DU1(IZ,3)-AJ1R*DU3(IZ,3))
     $ +C2T(IX,IZ,IR)*(DAJ2*DU1(IZ,4)+AJ2R*DU3(IZ,4))
     $                -AJ1K*DU1(IZ,5)

      EJ2(IX,IR,IZ)=EJ2(IX,IR,IZ)+
     $  C2T(IX,IZ,IR)*(AJ2R*DU1(IZ,1)-DAJ2*DU3(IZ,1))
     $  -ST(IX,IZ,IR)*(AJ1R*DU1(IZ,2)+AJKR*DU3(IZ,2))
     $  +CT(IX,IZ,IR)*(AJ1R*DU1(IZ,3)-AJKR*DU3(IZ,3))
     $ -S2T(IX,IZ,IR)*(AJ2R*DU1(IZ,4)+DAJ2*DU3(IZ,4))

  580 EJ3(IX,IR,IZ)=EJ3(IX,IR,IZ)+
     $  AJ2*(S2T(IX,IZ,IR)*DU2(IZ,1)+C2T(IX,IZ,IR)*DU2(IZ,4))
     $ +AJ1(IX,IZ,IR)*(CT(IX,IZ,IR)*DU2(IZ,2)+ST(IX,IZ,IR)*DU2(IZ,3))
     $ +AJ0(IX,IZ,IR)*DU2(IZ,5)

c ** Test of convergence of the wavenumber series: If the current
c terms of the series (for both the vertical and the horizontal
c displacements) are smaller than the current sum of the series
c times a small epsilon value (eps), convergence is assumed and
c the calculation of the wavenumber series is stopped.

      IZ=NZK
      DO 581 IR=1,NR
      DO 581 IX=1,NX,10
      A10=EPS*ABS(EJ1(IX,IR,IZ))
      A11=EPS*ABS(EJ2(IX,IR,IZ))
      A12=EPS*ABS(EJ3(IX,IR,IZ))

      AJ1R=AJ1(IX,IZ,IR)/R(IX,IZ,IR)
      AJKR=AK*AJ0(IX,IZ,IR)-AJ1R
      AJ2=(AJ1R+AJ1R)/AK-AJ0(IX,IZ,IR)
      AJ2R=(AJ2+AJ2)/R(IX,IZ,IR)
      AJ1K=AK*AJ1(IX,IZ,IR)
      DAJ2=AJ1K-AJ2R

      C1=
     $  S2T(IX,IZ,IR)*(DAJ2*DU1(IZ,1)-AJ2R*DU3(IZ,1))
     $  +CT(IX,IZ,IR)*(AJKR*DU1(IZ,2)+AJ1R*DU3(IZ,2))
     $  +ST(IX,IZ,IR)*(AJKR*DU1(IZ,3)-AJ1R*DU3(IZ,3))
     $ +C2T(IX,IZ,IR)*(DAJ2*DU1(IZ,4)+AJ2R*DU3(IZ,4))
     $                -AJ1K*DU1(IZ,5)

      C2=
     $  C2T(IX,IZ,IR)*(AJ2R*DU1(IZ,1)-DAJ2*DU3(IZ,1))
     $  -ST(IX,IZ,IR)*(AJ1R*DU1(IZ,2)+AJKR*DU3(IZ,2))
     $  +CT(IX,IZ,IR)*(AJ1R*DU1(IZ,3)-AJKR*DU3(IZ,3))
     $ -S2T(IX,IZ,IR)*(AJ2R*DU1(IZ,4)+DAJ2*DU3(IZ,4))

      C3=
     $  AJ2*(S2T(IX,IZ,IR)*DU2(IZ,1)+C2T(IX,IZ,IR)*DU2(IZ,4))
     $ +AJ1(IX,IZ,IR)*(CT(IX,IZ,IR)*DU2(IZ,2)+ST(IX,IZ,IR)*DU2(IZ,3))
     $ +AJ0(IX,IZ,IR)*DU2(IZ,5)

      A14=ABS(C1)
      A15=ABS(C2)
      A16=ABS(C3)
      IF((A14.GT.A10.OR.A15.GT.A11).OR.A16.GT.A12) GO TO 50
  581 CONTINUE
      IF(NZK.EQ.1) GO TO 20
      NLK=LAY(NZK)-1
      NLK1=NLK-1
      NZK=NZK-1
   50 CONTINUE 
   20 CONTINUE

      DO 586 IR=1,NR
      DO 585 J=1,3
  585 U(IF,IR,J)=0.D0
      DO 584 IX=1,NX
      DO 584 IZ=1,NZ
      C1=EJ1(IX,IR,IZ)*EV(IX,IZ)
      C2=EJ2(IX,IR,IZ)*EV(IX,IZ)
      C3=C1*CT(IX,IZ,IR)-C2*ST(IX,IZ,IR)
      C4=C2*CT(IX,IZ,IR)+C1*ST(IX,IZ,IR)
      C7=EJ3(IX,IR,IZ)*EV(IX,IZ)
      IF(RT) C3=C1
      IF(RT) C4=C2
      U(IF,IR,1)=U(IF,IR,1)+C3
      U(IF,IR,2)=U(IF,IR,2)+C4
 584  U(IF,IR,3)=U(IF,IR,3)+C7

 586  continue

****************************************************************
*	               output                                  *
****************************************************************

c writes the current frequency index and the number of wavenumbers
c considered in the wavenumber series:
      WRITE(6,*) IF,K

      DO 102 IZ=1,NZ                
      DO 102 IX=1,NX
  102 EV(IX,IZ)=EV(IX,IZ)*TEV(IX,IZ)

  100 FREQ=FREQ+DFREQ    

c ** End of loop on frequencies


c ** Calculation of the time domain solution:

      n3=ntime+ntime+3
      tex1=-aw*dt
      tex1=exp(tex1)
      ex7=1.0d0
      do 140 i=1,ntime
      yyy(i)=ex7
  140 ex7=ex7*tex1
      do 410 is=1,2
      do 410 ir=1,nr
      do 410 j=1,3
      do 517 i=1,n3
      y2(i)=0.   
  517 y(i)=0.
      freq=0.
      do 518 i=1,nfreq
         RW=PI2*FREQ                   
         OMEGA=CMPLX(RW,AW)
c Convolution of the impulse response with the source time
c function (=spectral multiplication) 
c + calculation of ITF of source function 
      c1=u(i,ir,j)*source(i,is)
      i2=i+i
      i1=i2-1
      q1=c1

      if ((ir.eq.1).and.(j.eq.1))then
         c2=source(i,is)*exp(-ai*omega*0*tsource)*dfreq
         q2=c2
         y2(i1)=real(q2)
         y2(i2)=aimag(q2)
         y2(n3-i1)=-y2(i2)
         y2(n3-i2)=y2(i1)
      endif

      y(i1)=real(q1)
      y(i2)=aimag(q1)
      y(n3-i1)=-y(i2)
      y(n3-i2)=y(i1)
 518  freq=freq+dfreq
c FFT:
      call four1(y,ntime,1)
      if ((ir.eq.1).and.(j.eq.1)) call four1(y2,ntime,1)
      do 411 i=1,ntime
c yy(i) corrects the FFT result for the time attenuation which results
c from the imaginary part of the frequency:
      sy(i,ir,j,is)=y(i+i-1)*yyy(i)-y(1)
      if ((ir.eq.1).and.(j.eq.1))
     & fctsou(i,is)=y2(i+i-1)*yyy(i)-y2(1)
  411 continue
  410 continue

c OUTPUT: **********************************************************
c - ir  :index of receiver considered
c - sy  :seismogram at receiver ir
c - j=1 :x-component of displacement (North)
c - j=2 :y-component of displacement (East)
c - j=3 :z-component of displacement (positive upward)
c - is=1 :displacement
c - is=2 :velocity
c
c
      ndir=LEN_TRIM(dir)
c     
!HUGO WRITES IN TWO CYCLES
      in=25
      do 505 is=2,2 
      do 505 ir=1,nr !must be nr in original code
       read(22,*)carstat
       cnord = dir(1:ndir) // "/" // carstat // "nord"
       cest = dir(1:ndir) // "/" // carstat // "est"
       cvert = dir(1:ndir) // "/" // carstat // "vert"
c      cnord=carstat//"nord"
c      cest=carstat//"est"
c      cvert=carstat//"vert"
       open (in,file=cnord)
       open (in+1,file=cest)
       open (in+2,file=cvert)
       
      do 504 j=1,3
      do 504 i=1,ntime   
      if (j.eq.1) write(in,*) sy(i,ir,1,is)
      if (j.eq.2) write(in+1,*) sy(i,ir,2,is)

cc  signe - pour le remettre vers le haut
      if (j.eq.3) write(in+2,*) -sy(i,ir,3,is)

      

 504  continue
       in=in+6
 505  continue

      
      do i=1,ntime
         write(18,*)fctsou(i,1)
         write(19,*)fctsou(i,2)
      enddo   

   70 format(10e12.5)


     
c ******************************************************************

      STOP     
      END     


c INVERTS A 2X2 MATRIX:
      subroutine inv2(a)
      complex*16 a(2,2),b(2,2),det
      do 1 i=1,2
      do 1 j=1,2
    1 b(i,j)=a(i,j)
      det=b(1,1)*b(2,2)-b(1,2)*b(2,1)
      a(1,1)=b(2,2)/det
      a(2,2)=b(1,1)/det
      a(1,2)=-b(1,2)/det
      a(2,1)=-b(2,1)/det
      return
      end

C/     ADD NAME=FF01AD          HSL     F77     DOUBLE
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS FF01AD
      SUBROUTINE FF01AD(VJ0,VY0,XD,N)
C  STANDARD FORTRAN 66(A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION VJ0,VY0,X,Y,Z,Q1,Q2,Q3,FX,X1,X2,X3,
     1                 X4,XD,XLG,A,B,C,D,E
      DIMENSION A(73),B(18),C(19),D(18),E(18)
      EQUIVALENCE (A(1),B(1)),(A(19),C(1)),(A(38),D(1)),(A(56),E(1))
      DATA XLG /1.0D+70/
      DATA B(1),B(2),B(3),B(4),B(5),B(6),B(7),B(8),B(9),B(10),B(11),
     1     B(12),B(13),B(14),B(15),B(16),B(17),B(18)    /
     1   -.17D-18                  , .1222D-16             ,
     2   -.75885D-15               , .4125321D-13          ,
     3   -.194383469D-11           , .7848696314D-10       ,
     4   -.267925353056D-8         , .7608163592419D-7     ,
     5   -.176194690776215D-5      , .3246032882100508D-4  ,
     6   -.46062616620627505D-3    , .48191800694676045D-2 ,
     7   -.34893769411408885D-1    , .15806710233209726D0  ,
     8   -.37009499387264978D0     , .26517861320333681D0  ,
     9   -.87234423528522213D-2    , .31545594294978024D0  /
      DATA C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),
     1     C(12),C(13),C(14),C(15),C(16),C(17),C(18),C(19)    /
     A   -.1D-19                   , .39D-18               ,
     B   -.2698D-16                , .164349D-14           ,
     C   -.8747341D-13             , .402633082D-11        ,
     D   -.15837552542D-9          , .524879478733D-8      ,
     E   -.14407233274019D-6       , .32065325376548D-5    ,
     F   -.5632079141056987D-4     , .75311359325777423D-3 ,
     G   -.72879624795520792D-2    , .47196689595763387D-1 ,
     H   -.17730201278114358D0     , .26156734625504664D0  ,
     I    .17903431407718266D0     ,-.27447430552974527D0  ,
     J   -.66292226406569883D-1     /
      DATA D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),D(10),D(11),
     1     D(12),D(13),D(14),D(15),D(16),D(17),D(18)    /
     K   -.1D-19                   , .2D-19                ,
     L   -.11D-18                  , .55D-18               ,
     M   -.288D-17                 , .1631D-16             ,
     N   -.10012D-15               , .67481D-15            ,
     O   -.506903D-14              , .4326596D-13          ,
     O   -.43045789D-12            , .516826239D-11        ,
     P   -.7864091377D-10          , .163064646352D-8      ,
     Q   -.5170594537606D-7        , .307518478751947D-5   ,
     R   -.53652204681321174D-3    , .19989206986950373D1 /
      DATA E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10),E(11),
     1     E(12),E(13),E(14),E(15),E(16),E(17),E(18)   /
     S    .1D-19                   ,-.3D-19                ,
     T    .13D-18                  ,-.62D-18               ,
     U    .311D-17                 ,-.1669D-16             ,
     V    .9662D-16                ,-.60999D-15            ,
     W    .425523D-14              ,-.3336328D-13          ,
     X    .30061451D-12            ,-.320674742D-11        ,
     Y    .4220121905D-10          ,-.72719159369D-9       ,
     Z    .1797245724797D-7        ,-.74144984110606D-6    ,
     1    .683851994261165D-4      ,-.31111709210674018D-1 /
      X=XD
      Y=DABS(X)
      Z=Y*.125D0
      IF(Z .LE.1.0D0)GO TO 10
      Z=1.0D0/Z
      X2=4.0D0*Z*Z-2.0D0
      N1=38
      N2=55
      GO TO 70
   10 IF(Z .EQ. 0.0D0)GO TO  78
      X2=4.0D0*Z*Z-2.0D0
      N1=1
      N2=18
   70 DO 80 J=1,2
      Q3=0.0D0
      Q2=0.0D0
      DO 40  I=N1,N2
      Q1=Q2
      Q2=Q3
   40 Q3=X2*Q2-Q1+A(I)
      FX=(Q3-Q1)*.5D0
      IF(N1-19)50,51,52
   50 VJ0=FX
      IF(N .LE. 0)GO TO 75
      N1=19
      N2=37
      GO TO 80
   52 IF(N1.EQ.56)GO TO 53
      X1=FX
      N1=56
      N2=73
   80 CONTINUE
   78 VJ0=1.0D0
      VY0=-XLG
      GO TO 75
   51 VY0=.6366197723675813D0*DLOG(Y)*VJ0+FX
      GO TO 75
   53 X2=DCOS(Y-.7853981633974483D0)
      X3=DSIN(Y-.7853981633974483D0)
      X4=.7978845608028654D0/DSQRT(Y)
      FX=FX*Z
      VJ0=X4*(X1*X2-FX*X3)
      VY0=X4*(FX*X2+X1*X3)
   75 RETURN
      END

C/     ADD NAME=FF02AD          HSL     F77     DOUBLE
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS FF02AD
      SUBROUTINE FF02AD(VJ1,VY1,XD,N)
C  STANDARD FORTRAN 66(A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION VJ1,VY1,X,Y,Z,Q1,Q2,Q3,FX,X1,X2,X3,X4,
     1                 XD,XLG,A,B,C,D,E
      DIMENSION A(72),B(18),C(18),D(18),E(18)
      EQUIVALENCE (A(1),B(1)),(A(19),C(1)),(A(37),D(1)),(A(55),E(1))
      DATA XLG/1.0D+70 /
      DATA B(1),B(2),B(3),B(4),B(5),B(6),B(7),B(8),B(9),B(10),B(11),
     1     B(12),B(13),B(14),B(15),B(16),B(17),B(18)   /
     1   -.4D-19                   , .295D-17              ,
     2   -.19554D-15               , .1138572D-13          ,
     3   -.57774042D-12            , .2528123664D-10       ,
     4   -.94242129816D-9          , .2949707007278D-7     ,
     5   -.76175878054003D-6       , .1588701923993213D-4  ,
     6   -.26044438934858068D-3    , .32402701826838575D-2 ,
     7   -.29175524806154208D-1    , .17770911723972828D0  ,
     8   -.66144393413454325D0     , .12879940988576776D1  ,
     9   -.11918011605412169D1     , .12967175412105298D1  /
      DATA C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),
     1     C(12),C(13),C(14),C(15),C(16),C(17),C(18)   /
     A    .9D-19                   ,-.658D-17              ,
     B    .42773D-15               ,-.2440949D-13          ,
     C    .121143321D-11           ,-.5172121473D-10       ,
     D    .187547032473D-8         ,-.5688440039919D-7     ,
     E    .141662436449235D-5      ,-.283046401495148D-4   ,
     F    .44047862986709951D-3    ,-.51316411610610848D-2 ,
     G    .42319180353336904D-1    ,-.22662499155675492D0  ,
     H    .67561578077218767D0     ,-.76729636288664594D0  ,
     I   -.12869738438135000D0     , .40608211771868508D-1 /
      DATA D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),D(10),D(11),
     1     D(12),D(13),D(14),D(15),D(16),D(17),D(18)   /
     J    .1D-19                   ,-.2D-19                ,
     K    .12D-18                  ,-.58D-18               ,
     L    .305D-17                 ,-.1731D-16             ,
     M    .10668D-15               ,-.72212D-15            ,
     N    .545267D-14              ,-.4684224D-13          ,
     O    .46991955D-12            ,-.570486364D-11        ,
     P    .881689866D-10           ,-.187189074911D-8      ,
     Q    .6177633960644D-7        ,-.398728430048891D-5   ,
     R    .89898983308594085D-3    , .20018060817200274D1  /
      DATA E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10),E(11),
     1     E(12),E(13),E(14),E(15),E(16),E(17),E(18)   /
     S   -.1D-19                   , .3D-19                ,
     T   -.14D-18                  , .65D-18               ,
     U   -.328D-17                 , .1768D-16             ,
     V   -.10269D-15               , .65083D-15            ,
     W   -.456125D-14              , .3596777D-13          ,
     X   -.32643157D-12            , .351521879D-11        ,
     Y   -.4686363688D-10          , .82291933277D-9       ,
     Z   -.2095978138408D-7        , .91386152579555D-6    ,
     1   -.9627723549157079D-4     , .93555574139070650D-1 /
      X=XD
      Y=DABS(X)
      Z=Y*.125D0
      IF(Z.LE.1.0D0)GO TO 10
      Z=1.0D0/Z
      X2=4.0D0*Z*Z-2.0D0
      N1=37
      N2=54
      GO TO 70
   10 IF(Z .LE. 0.0D0)GO TO 78
      X2=4.0D0*Z*Z-2.0D0
      N1=1
      N2=18
   70 DO 80 J=1,2
      Q3=0.0D0
      Q2=0.0D0
      DO 40 I=N1,N2
      Q1=Q2
      Q2=Q3
      Q3=X2*Q2-Q1+A(I)
   40 CONTINUE
      FX=(Q3-Q1)*.5D0
      IF(N1-19)50,51,52
   50 VJ1=FX*Z
      IF(N.LE.0)GO TO 75
      N1=19
      N2=36
      GO TO 80
   52 IF(N1.EQ.55)GO TO 53
      X1=FX
      N1=55
      N2=72
   80 CONTINUE
   78 VJ1=0.0D0
      VY1=-XLG
      GO TO 75
   51 VY1=.6366197723675813D0*(DLOG(Y)*VJ1-1.0D0/Y)+FX*Z
      GO TO 75
   53 X2=DCOS(Y-2.356194490192345D0)
      X3=DSIN(Y-2.356194490192345D0)
      X4=.7978845608028654D0/DSQRT(Y)
      FX=FX*Z
      VJ1=X4*(X1*X2-FX*X3)
      VY1=X4*(FX*X2+X1*X3)
   75 RETURN
      END

c  FFT
      subroutine four1(data,n,isign)
      dimension data(2050)
      ip0=2
      ip3=ip0*n
      i3rev=1
      do 50 i3=1,ip3,ip0
      if(i3-i3rev) 10,20,20
   10 tempr=data(i3)
      tempi=data(i3+1)
      data(i3)=data(i3rev)
      data(i3+1)=data(i3rev+1)
      data(i3rev)=tempr
      data(i3rev+1)=tempi
   20 ip1=ip3/2
   30 if(i3rev-ip1) 50,50,40
   40 i3rev=i3rev-ip1
      ip1=ip1/2
      if(ip1-ip0) 50,30,30
   50 i3rev=i3rev+ip1
      ip1=ip0
   60 if(ip1-ip3) 70,100,100
   70 ip2=ip1*2
      theta=6.283185307
      theta=theta/float(isign*ip2/ip0)
      sinth=sin(theta/2.0)
      wstpr=-2.0*sinth*sinth
      wstpi=sin(theta)
      wr=1.
      wi=0.
      do 90 i1=1,ip1,ip0
      do 80 i3=i1,ip3,ip2
      i2a=i3
      i2b=i2a+ip1
      tempr=wr*data(i2b)-wi*data(i2b+1)
      tempi=wr*data(i2b+1)+wi*data(i2b)
      data(i2b)=data(i2a)-tempr
      data(i2b+1)=data(i2a+1)-tempi
      data(i2a)=data(i2a)+tempr
   80 data(i2a+1)=data(i2a+1)+tempi
      tempr=wr
      wr=wr*wstpr-wi*wstpi+wr
   90 wi=wi*wstpr+tempr*wstpi+wi
      ip1=ip2
      go to 60
  100 return
      end

