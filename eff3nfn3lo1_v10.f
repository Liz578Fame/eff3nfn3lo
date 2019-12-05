      subroutine eff3nfn3lo1
c        effective three-nucleon force up to n3lo
c
c        the n2lo part is the version of 04/22/2012 
c        and based upon
c        Holt, Kaiser, Weise, PRC 81, 024002 (2010).
c
c        the n3lo part is to be added.
c
c
c        authors of code:
c                   r. machleidt et al.
c                   department of physics
c                   university of idaho
c                   machleid@uidaho.edu
c                   
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
c
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine.
c
c        specifications for these two common blocks
c
      logical heform,sing,trip,coup,endep
c
c
c        common block for all eff-subroutines
c
      common /ceff/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
      common /crrreff/ rrr
c
c
c        further specifications
c
      dimension vl(4),adminv(4,4),ldminv(4),mdminv(4)
      character*4 mesong(40)
      logical index
      logical indmg(40)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            'c   ','ss  ','ls  ','sq  ','sk  ',
     3            'sl  ',
     4         24*'    '/
      data index/.false./
      data indmg/40*.false./
      data jj/-1/
      data pi/3.141592653589793d0/
      save
c
c
c
c
      inter=1
c
c
c
c
c        call subroutine effpar once and only once
c
c
      if (index) go to 30
      index=.true.
c
c        if you want the potential to be zero for very large momenta,
c        choose rrr=1000.
c
      rrr=1000.
c
c
      call effpar
c
c
      iftgo=ift(inter)+1
      dwn=1.d0/wnn(inter)
c
c
      iman=imaa(inter)
      imen=imea(inter)
c
c
      ez1=ezz1(inter)
      ez2=ezz2(inter)
c
c
c        prepare constant over-all factor
c
      fac=pi/(2.d0*pi)**3*dwn*dwn
c     ---------------------------
c
c
c
c
   30 if (j.eq.jj) go to 50
      jj=j
      if (j.eq.0) go to 50
      aj=dble(j)
      aj1=dble(j+1)
      a2j1=dble(2*j+1)
      aaj6=dsqrt(aj*aj1)
c
c        coefficient matrix for the translations into lsj formalism
c
      adminv(1,1)=aj1
      adminv(1,2)=aj
      adminv(1,3)=-aaj6
      adminv(1,4)=-aaj6
      adminv(2,1)=aj
      adminv(2,2)=aj1
      adminv(2,3)=aaj6
      adminv(2,4)=aaj6
      adminv(3,1)=aaj6
      adminv(3,2)=-aaj6
      adminv(3,3)=aj1
      adminv(3,4)=-aj
      adminv(4,1)=aaj6
      adminv(4,2)=-aaj6
      adminv(4,3)=-aj
      adminv(4,4)=aj1
c
c       inversion
c
      call dminveff (adminv,4,deter,ldminv,mdminv)
c
c
c
c
c        prepare expressions depending on x and y
c        ----------------------------------------
c        ----------------------------------------
c
c
c
c
   50 xa=xmev*dwn
      ya=ymev*dwn
      indxy=.false.
      x=xa
      xx=x*x
      y=ya
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=dsqrt(1.d0+xx)
      ey=dsqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c
c
      xy=xy2*0.5d0
      ee=ex*ey
      ree=dsqrt(ee)
      eem1=ee-1.d0
      eme=ex-ey
      emeh=eme*0.5d0
      emehq=emeh*emeh
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
c
      xxpyyh=xxpyy*0.5d0
      xy3=xy*3.d0
      xy4=xy*4.d0
c
c
c
c
      do 63 iv=1,6
   63 v(iv)=0.d0
      do 65 il=iman,imen
      do 65 iv=1,32
   65 vj(iv,il)=0.d0
c
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72,71,72,75,76),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 80
c
c        minimal relativity
c
   71 fff=fac/ree
      go to 80
c
c        factor m/e*m/e
c
   72 fff=fac/ee
      go to 80
c
c        sharp cutoff
c
   75 if (xmev.gt.ez1.or.ymev.gt.ez1) then
      return
      else
      fff=fac
      end if
      go to 80
c
c        exponential cutoff
c
   76 expo=(xmev/ez1)**(2.d0*ez2)+(ymev/ez1)**(2.d0*ez2)
      if (expo.gt.rrr) then
      expo=rrr
      end if
      fff=fac*dexp(-expo)
c
c
   80 continue
c
c
c
c
c        contributions
c        -------------
c        -------------
c
c
c
c
      do 5995 img=1,mge
      mg=mggo(img,inter)
      if (mg.gt.16) go to 9000
      if (mg.eq.0) go to 8000
      me=mgg(mg,inter)
      go to (9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,
     1       1100,1200,1300,1400,1500,1600),mg
c
c
c
c
c        c   ,central force
c        ------------------
c
c
c
c
 1100 mc=1
c
      ff=1.d0
      f(1)=2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
c
      call effstr(1,1,me)
      go to 5995
c
c
c
c
c        ss  , spin-spin force
c        ---------------------
c
c
c
c
 1200 mc=1
c
      ff=1.d0
      f(1)=-6.d0
      f(2)=0.d0
      f(3)=2.d0
      f(4)=0.d0
      f(5)=0.d0
      f(6)=f(3)
      f(7)=-f(3)
      f(8)=f(7)
c
      call effstr(1,1,me)
      go to 5995
c
c
c
c
c        ls  , spin-orbit force
c        ----------------------
c
c
c
c
 1300 mc=1
c
      ff=1.d0
      f(1)=0.d0
      f(2)=0.d0
      f(3)=0.d0
      f(4)=-xy2
      f(5)=-xy2
      f(6)=0.d0
      f(7)=0.d0
      f(8)=0.d0
      f(9)=0.d0
      f(10)=+xy2
      f(11)=-xy2
c
      call effstr(2,1,me)
      go to 5995
c
c
c
c
c        sq  , sq tensor force (where q denotes the momentum transfer)
c        ---------------------
c
c
c
c
 1400 mc=1
c
      ff=1.d0
      f(1)=-xxpyy*2.0d0
      f(2)=xy*4.d0
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=(xx-yy)*2.0d0
      f(8)=-f(7)
c
      call effstr(1,1,me)
      go to 5995
c
c
c
c
c        sk  , sk tensor force (where k denotes the average momentum)
c        ---------------------
c
c
c
c
 1500 mc=1
c
      ff=0.25d0
      f(1)=-xxpyy*2.0d0
      f(2)=-xy*4.d0
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=(xx-yy)*2.0d0
      f(8)=-f(7)
c
      call effstr(1,1,me)
      go to 5995
c
c
c
c
c        sl  , "quadratic spin-orbit force"
c               or sigma-l operator
c        ----------------------------------
c
c
c
c
 1600 mc=1
c
      ff=1.d0
      f(1)=-xxyy*2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=-f(1)
      f(7)=f(1)
      f(8)=f(7)
      f(9)=f(6)*2.d0
c
      call effstr(4,1,me)
      go to 5995
c
c
c
c
c
c        this has been the end of the contributions
c        ------------------------------------------
c
c
c
c
c        errors and warnings
c        -------------------
c
c
c
c
 9000 if (indmg(mg)) go to 5995
      write (kwrite,19000) mesong(mg)
19000 format(1h ////' warng. in eff3nf: contribution ',a4,' does not exi
     1st in this program.'/' contribution ignored. execution continued.'
     2////)
      indmg(mg)=.true.
c
c
c
c
 5995 continue
c
c
c
c
c        add up contributions
c        --------------------
c
c
c
c
 8000 continue
c
c
c
c        add up terms
c        -----------
c
      do 8155 il=iman,imen
      do 8155 iv=1,6
 8155 v(iv)=v(iv)+vj(iv,il)
c
c
c
c
      if (j.eq.0.or..not.heform) go to 8900
c
c
c         translation into (combinations of) helicity states
c
c
      do 8505 i=1,4
 8505 vl(i)=v(i+2)
c
      do 8520 ii=1,4
      iii=ii+2
      v(iii)=0.d0
c
      do 8515 i=1,4
 8515 v(iii)=v(iii)+adminv(ii,i)*vl(i)
 8520 v(iii)=v(iii)*a2j1
c
c
c
c
 8900 return
      end
      subroutine effpar
c
c        effpar reads, writes, and stores the parameter for all 
c        eff-subroutines.
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all eff-subroutines
c
      common /ceff/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
      common /cpareff/ cb1,cb2,cb3,cb4,css,ctt,cdd,cee,qf
c
c
c
c        further specifications
c
      dimension cc(5),cca(5)
      integer name(3),nname(15)
      integer imga(3)
      integer cut,cuta,fun,end
      integer mesong(40)
      logical index
      logical indca
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            'c   ','ss  ','ls  ','sq  ','sk  ',
     3            'sl  ',
     4         24*'    '/
      data index/.false./
      data pi/3.141592653589793d0/
      data uf/197.327d0/
      data cut/'cut '/,cuta/'cuta'/
      data fun/'fun '/,end/'end '/
      save
c
c
c
c
c
10000 format (2a4,a2,15a4)
10001 format (1h )
10002 format (1h /' contribution  par_1     par_2     mass     iso-spin
     1    iprop'/9x,'fun-type  f u n c t i o n   p a r a m e t e r s')
10003 format (2a4,a2,5f10.6)
10022 format (2a4,a2,6f10.6)
10023 format (1h, 2a4,a2,6f10.6)
10004 format (1h ,2a4,a2,f12.6,f10.6,1x,f10.5,2(f7.1,3x))
10005 format (1h ,2a4,a2,f4.1,1x,2f10.5,f13.5,f10.5)
10006 format (2a4,a2,3i3)
10007 format (1h ,2a4,a2,3i3)
10008 format (1h ,60(1h-))
10009 format (2a4,a2,i3,2f10.2)
10010 format (1h ,2a4,a2,i3,2f10.2)
10011 format (//' eff3nfn3lo: effective three-nucleon force up to n3lo')
10015 format (' input-parameter-set:'/1h ,20(1h-))
10016 format (1h ,2a4,a2,15a4)
10020 format (1h ,2a4,a2,f12.6,4f9.2)
10021 format (1h ,2a4,a2,5f10.6)
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
c
c
c        maxima of certain indices related to the dimension as follows:
c        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
c                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
c                  ima(mee,mge,3)
c
      mge=40
      mee=30
      mme=20
      mice=20
      mindce=2
      imb=1
      ime=0
      imee=270
c        mme always ge mice, mindce
c
c        set all parameters and indices to zero or .false.
c
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
c
c
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0.d0
      endep=.false.
c
c
      pi2=pi*pi
      pi4=pi2*pi2
c
c
c
c
c        reading and writing of first 4 cards
c        --------------------------------------
c        --------------------------------------
c
c
c
c        write headline and read and write name of parameter set
c
   50 continue
      write (kwrite,10011)
      write (kwrite,10008)
      write (kwrite,10015)
      read  (kread, 10000) name,nname
      write (kwrite,10016) name,nname
      if (inter.eq.1) label=name(1)
      indpar(inter)=.true.
      indca=.false.
c
c        read and write index-parameter concerning the factor for the
c        whole potential and cutoff mass
c
c
      read  (kread, 10009) name,ift(inter),ezz1(inter),ezz2(inter)
      write (kwrite,10010) name,ift(inter),ezz1(inter),ezz2(inter)
      iftyp=ift(inter)
      if (iftyp.lt.0.or.iftyp.gt.6) go to 9003
c
c        read and write parameters for numerical integration
c
      read  (kread, 10006) name,mint(inter),maxt(inter)
      write (kwrite,10007) name,mint(inter),maxt(inter)
c
c        read and write mass of nucleon
c
      read  (kread, 10003) name,wn,qf
      write (kwrite,10004) name,wn,qf
      wnq=wn*wn
      dwn=1.d0/wn
      dwnq=dwn*dwn
      wnn(inter)=wn
      qf=qf*uf*dwn
      rho=2.d0*qf**3/(3.d0*pi2)
c
c        read and write ga and fpi
c
      read  (kread, 10003) name,ga,fpi
      write (kwrite,10004) name,ga,fpi
      ga2=ga*ga
      ga4=ga2*ga2
      ga6=ga4*ga2
      fpi=fpi*dwn
      fpi2=fpi*fpi
      fpi4=fpi2*fpi2
      fpi6=fpi4*fpi2
      gaa(inter)=ga2
      fpia(inter)=fpi2
c
c        read and write LECs of the pi-N Lagrangian
c
c        the c_i LECs
c     read  (kread, 10003) name,cc
c     write (kwrite,10021) name,cc
      read  (kread, 10022) name,(cc(i), i=1,4), css,ctt 
      write (kwrite, 10023) name,(cc(i), i=1,4), css,ctt 
      cb1=cc(1)*wn*1.d-3
      cb2=cc(2)*wn*1.d-3
      cb3=cc(3)*wn*1.d-3
      cb4=cc(4)*wn*1.d-3
      css=css*wnq*1.d-2 
      ctt=ctt*wnq*1.d-2 
c
c        the LECs cd and ce and their lambda
      read  (kread, 10003) name,cc
      write (kwrite,10004) name,cc
      cdd=cc(1)*wn/(fpi2*cc(3))
      cee=cc(2)*wn/(fpi4*cc(3))
c
c
c
c        write headline for parameters
c
      write (kwrite,10002)
      write (kwrite,10008)
c
c
c
c        read, write and store parameters
c        --------------------------------
c        --------------------------------
c
c
c
   61 read  (kread, 10003) name,cc
c
c        check if end of input
c
      if (name(1).eq.end) go to 7000
c
c        check if data-card just read contains cut-off or
c        function parameters
c
      if (name(1).eq.cut.or.name(1).eq.fun) go to 70
c
      if (name(1).eq.cuta) then
      write (kwrite,10005) name,cc
      indca=.true.
      do i=1,5
      cca(i)=cc(i)
      end do
      go to 61
      end if

c
c
c
c
c        write parameters which are no cut-off or function parameters
c        ------------------------------------------------------------
c
c
c
c
      write (kwrite,10004) name,cc
c
c
c        find out number of contribution mg
c
      do 63 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 64
   63 continue
      go to 9000
c
c
c
c
c        store parameters which are no cut-off or function parameters
c        ------------------------------------------------------------
c
c
c
c
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 65
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   65 continue
c
c
      c(1,ime)=cc(1)
c
c
      if (mg.ge.11.and.cc(2).ne.0.d0) then
      c(1,ime)=(cc(1)/(2.d0*cc(2))*wn)**2
      if (cc(1).lt.0.d0) c(1,ime)=-c(1,ime)
      end if
c
c
c        store meson mass square in units of nucleon mass square
      c(4,ime)=cc(3)*cc(3)*dwnq
c
c        test iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
c         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
c        store and test iprsp
      icc=cc(5)
      ic(1,ime)=icc
      if (iabs(ic(1,ime)).gt.1) go to 9005
c
c        index values for further storing
      mi=4
      mm=5
c
c
c        check if there is a `cutall' cutoff
c
      if (indca) then
      name(1)=cut
      do i=1,5
      cc(i)=cca(i)
      end do
      go to 72
      else
      go to 61
      end if
c
c
c
c
c        write cut-off or function parameters
c        ------------------------------------
c
c
c
c
   70 write (kwrite,10005) name,cc
c
c
   72 continue
c
c
c
c
c        store parameters
c        ----------------
c
c
c
      ityp=cc(1)
c
      if (ityp.eq.0) go to 5995
      if (ityp.lt.1.or.ityp.gt.38) go to 9002
c
      im=ime
c
c        store typ of cut-off or function
      ic(mi,im)=ityp
c
      if (ityp.le.10) then
c        store and test typ of propagator of cut-off
      ic(mi+1,im)=cc(2)
      if (ic(mi+1,im).lt.0.or.ic(mi+1,im).gt.1) go to 9006
      end if
c
      go to (100,100,300,9002,500,600,9002,9002,9002,9002,
     1 1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     2 2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     3 3100,3200,3300,3400,3500,3600,3700,3800),ityp
c
c
c
c
c        cut-off of dipole type
c        **********************
c
c
c        store and test exponent of cut-off
  100 ic(mi+2,im)=cc(3)
      if (ic(mi+2,im).lt.0) go to 9009
      if (ic(mi+2,im).gt.0) go to 101
c        exponent is zero, omit cut-off
      ic(mi,im)=0
      ic(mi+1,im)=0
      go to 5995
c        store cut-off mass for denominator
  101 c(mm+1,im)=cc(4)*cc(4)*dwnq
c        store numerator of cut-off
      c(mm,im)=c(mm+1,im)
      if (ityp.eq.2)     c(mm,im)=c(mm,im)-c(4,im)
      mi=mi+3
      mm=mm+2
      go to 5995
c
c
c
c
c        exponential form factor of momentum transfer
c        ********************************************
c
c
c        check exponent
  300 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 301
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 5995
c        store exponent
  301 c(mm+1,im)=cc(3)
c        compute constant factor for argument of exponential function
      c(mm,im)=wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+2
      go to 5995
c
c
c
c
c        sharp cutoff in x and y
c        ***********************
c
c
  500 c(mm,im)=cc(4)*dwn
      mi=mi+2
      mm=mm+1
      go to 5995
c
c
c
c
c        exponential form factor of xx and yy
c        ************************************
c
c
c        check exponent
  600 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 601
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 5995
c        store exponent
  601 c(mm+1,im)=cc(3)
c        compute constant factor for argument of exponential function
      c(mm,im)=wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+2
      go to 5995
c
c
c
c
c        function q^2 (momentum-transfer squared)
c        ************
c
c
 1100 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function k^2 (average-momentum squared)
c        ************



c
c
 1200 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,1
c        *******************
c
c
 1300 c(mm,im)=rho*ga2/(2.d0*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,2
c        *******************
c
c
 1400 c(mm,im)=ga2/(8.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1500 c(mm,im)=ga2/(16.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1600 c(mm,im)=ga2/(16.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1700 c(mm,im)=ga2/(16.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1800 c(mm,im)=ga2/(16.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1900 c(mm,im)=ga2/(16.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,4
c        *******************
c
c
 2000 c(mm,im)=-ga*cdd*rho/(8.d0*fpi2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2100 c(mm,im)=ga*cdd/(16.d0*pi2*fpi2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2200 c(mm,im)=ga*cdd/(16.d0*pi2*fpi2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2300 c(mm,im)=-ga*cdd/(16.d0*pi2*fpi2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2400 c(mm,im)=ga*cdd/(16.d0*pi2*fpi2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function for vmed,6
c        *******************
c
c
 2500 c(mm,im)=-3.d0*cee*rho/2.d0
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c

c        function for vmed,7
c        *******************
c
c
cc 2600 c(mm,im)=ga2*ctt*rho/(8.d0*pi*fpi4)
cc      mi=mi+1
cc      mm=mm+1
cc      go to 5995
c
c
c        function for vmed,7
c        *******************
c
c
 2600 c(mm,im)=ga2*ctt/(12.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c    function for vmed,8
c    *******************
 2700  c(mm,im)= ga2*ctt*3.d0/(720.d0*pi2*pi*fpi4)
       mi=mi+1
       mm=mm+1
       go to 5995
c
c
c    function for vmed,8
c    *******************
 2800  c(mm,im)= ga2*ctt/(720.d0*pi2*pi*fpi4)
       mi=mi+1
       mm=mm+1
       go to 5995
c
c
c
c
c    function for vmed,9
c    *******************
 2900  c(mm,im)= ga4*ctt*3.d0/(24.d0*pi2*pi*fpi4)
       mi=mi+1
       mm=mm+1
       go to 5995
c
c
c
c
c    function for vmed,9
c    *******************
 3000  c(mm,im)= ga4*ctt*(-3.d0)/(24.d0*pi2*pi*fpi4)
       mi=mi+1
       mm=mm+1
       go to 5995
c
c
c
c
c
c    function for vmed,9
c    *******************
 3100  c(mm,im)= ga4*ctt/(12.d0*pi2*pi*fpi4)
       mi=mi+1
       mm=mm+1
       go to 5995
c
c
c
c
c       function for vmed,10
c        *******************
c
c
 3200 c(mm,im)=ga4*ctt/(80.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c       function for vmed,10
c        *******************
c
 3300 c(mm,im)=ga4*ctt*(-3.d0)/(360.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c       function for vmed,10
c        *******************
c
 3400 c(mm,im)=ga4*ctt*(-1.d0)/(360.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c       function for vmed,10
c        *******************
c
 3500 c(mm,im)=ga4*ctt/(64.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c
c
c       function for vmed,10
c        *******************
c
 3600 c(mm,im)=(ga4*ctt*(-2.d0))/(64.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c
c
c       function for vmed,10
c        *******************
c
 3700 c(mm,im)=ga4*ctt*2.d0/(64.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c
c       function for vmed,10
c        *******************
c
 3800 c(mm,im)=ga4*ctt*2.d0/(64.d0*pi2*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c

c        end cut-offs and functions
c        **************************
c
c        test dimensions
 5995 if (mi.gt.mice.or.mm.gt.mme) go to 9010
c
      go to 61


c
c
c        conclusions
c        -----------
c        -----------
c
c
c        write end
 7000 write (kwrite,10004) name
      write (kwrite,10008)
      write (kwrite,10008)
c
      imaa(inter)=imb
      imea(inter)=ime
      imb=ime+1
c
      return
c
c
c
c
c        errors
c        ------
c        ------
c
c
 9000 write (kwrite,19000) name(1)
19000 format (1h ////' error in effpar:  contribution  ',a4,'   does not
     1 exist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9001 write (kwrite,19001)
19001 format (1h ////' error in effpar:too many contributions within a g
     1roup with respect to'/' the given dimensions. execution terminated
     2.'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1)
19002 format (1h ////' error in effpar: cut/fun typ',f10.4,'  does not e
     1xist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9003 write (kwrite,19003) iftyp
19003 format (1h ////' error in effpar: factor typ has the non-permissib
     1le value',i4,' .'/' execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (1h ////' error in effpar: isospin has the non-permissible
     1value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9005 write (kwrite,19005) cc(5)
19005 format (1h ////' error in effpar: iprop/spe has the non-permissibl
     1e value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9006 write (kwrite,19006) cc(2)
19006 format (1h ////' error in effpar: the index for the propagator of
     1the cut-off has the'/' non-permissible value',f10.4,'  . execution
     2 terminated.'////)
      go to 9999
c
c
 9009 write (kwrite,19009)
19009 format (1h ////' error in effpar: the exponent of the cut-off is l
     1ess than zero.'/' execution terminated.'////)
      go to 9999
c
c
 9010 write (kwrite,19010)
19010 format (1h ////' error in effpar: too many cut/fun parameters with
     1 respect to the given'/' dimensions. execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (1h ////' error in effpar:  too many contr. with respect to
     1 the dimensions given'/' to this program. execution terminated.'
     2////)
      go to 9999
c
c
 9999 stop
      end
      subroutine effstr (icase,max,mex)
c
c        effstr computes the structure of contributions to the nuclear force
c
c
      implicit real*8 (a-h,o-z)
c
c
c        common blocks
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all eff-subroutines
c
      common /ceff/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c     further specifications
c
      dimension vv(32)
      dimension tt(2,3)
      logical index
      logical indiso
      data jj/-1/
      data index/.false./
      save
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      tt(1,1)=1.d0
      tt(2,1)=-3.d0
c
      do 1 ii=2,3
      do 1 i=1,2
    1 tt(i,ii)=1.d0
c
c
c
c
c
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
c
c
      if (mc.ne.1) go to 60
c
c
c
c
c        call integrals
c        --------------
c
c
c
c
      call effai
c
c
c
c
   60 if (mc.lt.1) mc=1
c
      if (c(mc,im).eq.0.d0) go to 1095
c
c
c
c
c        nn-nn helicity amplitudes /combinations/
c        ----------------------------------------
c
c
c
c
c        basic structure (a factor of 2 is included in v5 and v6)
c
c
      ive=6
c
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
c
c
      go to (1000,120,130,140),icase
c
c
c        additional terms required for the tensor coupling
c        of the rho-meson or for certain operators,
c        like, the spin-orbit operator (`ls  ')
c
c
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
c
c
c        additional terms in case of 2+ mesons
c        not needed here
c
c
  130 Continue
      go to 1000
c
c
c        additional terms needed for the sigma-l operator (`sl  ')
c
c
  140 vv(1)=vv(1)+f(6)*ai(5,m)
      vv(2)=vv(2)+f(1)*ai(5,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(1)*ai(11,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(1)*ai(12,m)
      vv(5)=vv(5)+f(6)*ai(13,m)
      vv(6)=vv(6)+f(6)*ai(13,m)
c
c
c
c
 1000 continue
c
c
c
c
c        set certain cases to zero in case of inter=1
c
      if (j.ne.0) go to 1021
      vv(2)=0.d0
      vv(4)=0.d0
      vv(5)=0.d0
      vv(6)=0.d0
c
 1021 if (.not.sing) vv(1)=0.d0
      if (.not.trip) vv(2)=0.d0
      if (coup) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0.d0
c
 1030 continue
c
c
c        transformation into lsj-formalism
c 
      if (j.eq.jj) go to 1035
      jj=j
      aj=dfloat(j)
      aj1=dfloat(j+1)
      d2j1=1.d0/dfloat(2*j+1)
      arjj1=dsqrt(aj*aj1)
c
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=-arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34-aj1*v5+aj*v6)
      vv(6)=d2j1*(v34+aj*v5-aj1*v6)
c
c
c        possible different sign depending on the convention used
      vv(5)=-vv(5)
      vv(6)=-vv(6)
c
c
c
c
c        multiply with factors
c        ---------------------
c
c
c
c
 1040 is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
c
c        multiply with coupling-constant and factors fff and ff
c
      vv(iv)=vv(iv)*fc
c
c        multiply with isospin factor
c
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
c
c     add up in case of several couplings for one meson-exchange
c     and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
c
c
 1095 continue
c
c
      return
      end
      subroutine effai
c
c        effai integrates over theta
c
c
      implicit real*8 (a-h,o-z)
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all eff-subroutines
c
      common /ceff/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
      dimension gi(7)
c
      dimension pj(7,96)
      real*4 axy2,aomq,am
      logical indj
      data ige/7/
      data nnt/-1/,iinter/-1/,jj/-1/
      save
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
      min=mint(inter)
      max=maxt(inter)
c
      igeint=7
c
      wn=wnn(inter)
      dwn=1.d0/wn
      wnq=wn*wn
c
c
c
c
   60 if (j.eq.jj) go to 70
      jj=j
      indj=.false.
c
c
      aj=dfloat(j)
      aj1=dfloat(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=dsqrt(ajdj1)
c
c
      aj2=dfloat(j+2)
      ajm1=dfloat(j-1)
c
c
      ajj1=aj*aj1
      ajj2=ajm1*aj2
      ajjb=aj*ajm1
c
      aajj=0.d0
      if (j.gt.1)
     1aajj=aj/dsqrt(ajj1*ajj2)
c
      aaj1=aajj*ajm1
      aaj2=aajj*aj1
      aaj3=aajj*2.d0
c
      if (j.gt.1) go to 62
      aajj=0.d0
      go to 63
   62 aajj=1.d0/(aj1*dsqrt(ajj2))
c
   63 aaj4=aajj*ajjb
      aaj5=aajj*aj1*2.d0
      aaj6=aajj*(ajj1+2.d0)
      aaj7=aajj*ajj2
c
c
c
c
c        find out appropriate number of gauss-points
c        -------------------------------------------
c
c
   70 c4=c(4,im)
      if (c4.eq.0.d0) then
      c4=(138.*dwn)**2
      end if
      iprsp=ic(1,im)
c
c
c        compute am
c
      axy2=xy2
      if (iprsp.ne.1) go to 91
      aomq=eem12+c4
      go to 92
   91 aomq=xxpyy+c4
c
   92 am=axy2/aomq
c
c
c        compute number of gausspoints (nt)
c
c
      if (am.gt.0.999) go to 94
c
c
      if (am.gt.0.85) am=am**(-alog(1.-am)-0.9)
c
c
      nt=float(min)/(1.-am)+0.9
c
c
      if (nt.gt.max) nt=max
      go to 95
c
c
   94 nt=max
c
c
   95 nt=nt+j
c
c        compute nt, which is suitable for gset
c
      if (nt.le.16) go to 98
      if (nt.gt.24) go to 96
      nt=4*(nt/4)
      go to 98
   96 if (nt.gt.48) go to 97
      nt=8*(nt/8)
      go to 98
   97 nt=16*(nt/16)
      if (nt.gt.96) nt=96
c
   98 if (nt.eq.nnt.and.indj) go to 100
c
c
c
c
c        call gauss-points
c        -----------------
c
c
c
c
      call gseteff (-1.d0,1.d0,nt,ct,wt)
      nnt=nt
c
c
c
c
c        call legendre-polynoms if necessary
c        -----------------------------------
c
c
c
c
      indj=.true.
      do 99 i=1,nt
      t=ct(i)
      call legpeff (pj(1,i),pj(3,i),t,j)
      pj(2,i)=pj(1,i)*t
      pj(4,i)=pj(2,i)*t
      pj(6,i)=pj(4,i)*t
      pj(5,i)=pj(3,i)*t
   99 pj(7,i)=pj(5,i)*t
c
c
c
c
c        call integrand
c        --------------
c
c
c
c
  100 call effaa
c
c
c
c
c        prepare for integration
c
c
c
c
      do 2001 ig=1,igeint
 2001 gi(ig)=0.d0
c
c
c
c
c        integration-loop of theta
c        -------------------------
c
c
c
c
      do 2005 i=1,nt
      do 2005 ig=1,igeint
 2005 gi(ig)=gi(ig)+pj(ig,i)*aa(i)
c
c
c
      if (j.ne.0) go to 2010
      gi(3)=0.d0
      gi(5)=0.d0
      gi(7)=0.d0
c
c
c
c
c        combinations of integrals
c        -------------------------
c
c
c
c
 2010 ai(1,m)=gi(1)
c
      ai(2,m)=gi(2)
      ai(3,m)= ajdj1*gi(2)+dj1*gi(3)
      gi23m  =gi(2)-gi(3)
      ai(4,m)=aaj*gi23m
c
c
      ai(5,m)=gi(4)
      ai(6,m)= ajdj1*gi(4)+dj1*gi(5)
      gi45m  =gi(4)-gi(5)
      ai(7,m)=aaj*gi45m
c
c
      ai( 8,m)= aaj1*gi(4)-aaj2*gi(1)+aaj3*gi(5)
      aai1    = aaj4*gi(4)+aaj5*gi(1)-aaj6*gi(5)
      aai2    = aaj7*gi23m
      ai( 9,m)= aai2+aai1
      ai(10,m)= aai2-aai1
c
c
      ai(11,m)=gi(6)
      ai(12,m)=ajdj1*gi(6)+dj1*gi(7)
      ai(13,m)=aaj*(gi(6)-gi(7))
c
c
      return
      end
      subroutine effaa
c
c        effaa computes propagators, cutoffs, and functions
c
c
      implicit real*8 (a-h,o-z)
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all eff-subroutines
c
      common /ceff/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
      common /cpareff/ cb1,cb2,cb3,cb4,css,ctt,cdd,cee,qf
c
c
      common /crrreff/ rrr
c
c
c
c        further specifications
      dimension deltaq(96,7)
      dimension ak(96),wk(96)
      dimension gint0(96),gints(96),gintss(96)
      dimension gint1(96),gint1s(96),gint2(96),gint3(96)
      data iinter/-1/
      data nnt/-1/
      data cc4/-1.d0/
      data pi/3.141592653589793d0/
      save
c
c
c
c
      if (inter.eq.iinter) go to 10
      iinter=inter
      ga2=gaa(inter)
      ga4=ga2*ga2
      fpi2=fpia(inter)
      pi2=pi*pi
   10 continue
c
c
c        prepare functions of x, y, and c4, but no angle
c        -----------------------------------------------
c
c
      c4=c(4,im)
      if (indxy.and.cc4.eq.c4) go to 30
      indxy=.true.
      cc4=c4
      rootc4=dsqrt(c4)
      yyoff=0.5d0*(xx+yy)
      yoff=dsqrt(yyoff)
      qfq=qf*qf
      qf3=qfq*qf
      
c
c
      gamma0=qf-rootc4*(datan((qf+yoff)/rootc4)
     1      +datan((qf-yoff)/rootc4))
     2      +(c4+qfq-yyoff)/(4.d0*yoff)*
     3       dlog((c4+(qf+yoff)**2)/(c4+(qf-yoff)**2))
      gamma1=qf/(4.d0*yyoff)*(c4+qfq+yyoff)-gamma0
     1      -1.d0/(16.d0*yyoff*yoff)*(c4+(qf+yoff)**2)*
     2       (c4+(qf-yoff)**2)*
     3       dlog((c4+(qf+yoff)**2)/(c4+(qf-yoff)**2))
      gamma2=qf3/9.d0+1.d0/6.d0*(qfq-c4-yyoff)*gamma0
     1      +1.d0/6.d0*(c4+qfq-yyoff)*gamma1
      gamma3=qf3/(3.d0*yyoff)
     2      -(c4+qfq+yyoff)/(2.d0*yyoff)*gamma0
     3      -(c4+qfq+3.d0*yyoff)/(2.d0*yyoff)*gamma1
c
      gam013=gamma0+2.d0*gamma1+gamma3
     


c       Adding H-functions, H_(0-3)
c      ***********************************
     

       rc4=dsqrt(c4)
       c4r3=c4*rc4
       c4r4=c4*c4
       c4r6=c4*c4r4
       pkp=yoff+qf
       pkp2=pkp**2
       pkp4=pkp**4
       pkm=yoff-qf
       pkm2=pkm**2
       pkm4=pkm**4
       qf4=qfq**2
       qf5=qf*qf4
       qf7=qf5*qfq
       y3off=yyoff*yoff
       y4off=yyoff**2
       y5off=y4off*yoff
       q2p10c4=qfq+(10.d0*c4)
       q2p2c4=qfq+(2.d0*c4)
       q2mc4=(13.d0*qfq)-(3.d0*c4)
       q4mp4=(15.d0*qf4)-y4off
       p4mq4=y4off-(15.d0*qf4)
       c4mp2=(6.d0*c4)-(19.d0*yyoff)
       p2mq2=(5.d0*yyoff)-qfq
       qp2mpq4=(6.d0*qfq*yyoff)-(5.d0*y4off)-qf4



       haitch0=(qf*c4*(yyoff+(21.d0*qfq)-(19.d0*c4)))
     1   +((c4r3/4.d0*yoff)*((35.d0*yyoff)-(35.d0*qfq)-(44.d0*c4)))*
     2    (dlog(((4.d0*c4)+pkp2)/((4.d0*c4)+pkm2)))
     3   +(((qf3/yoff)*q2p10c4)-(5.d0*yyoff*c4)+(5.d0*qf3*yoff)
     4   +((5.d0/2.d0)*yyoff*qfq)+(15.d0*c4*q2p2c4)
     5   +((1.d0/4.d0)*q4mp4))* (datan(pkp/(2.d0*c4)))
     6   +(((qf3/yoff)*q2p10c4)+(5.d0*yyoff*c4)+(5.d0*qf3*yoff)
     7   -((5.d0/2.d0)*yyoff*qfq)-(15.d0*c4*q2p2c4)
     8   +((1.d0/4.d0)*p4mq4))* (datan(pkm/(2.d0*c4)))

       haitch1=(qf*c4*(yyoff-(9.d0*qfq)+(36.d0*c4)))
     1   +(((2.d0*c4r3)/yoff)*((5.d0*qfq)+(12.d0*c4)-(5.d0*yyoff)))*
     2    (dlog(((4.d0*c4)+pkp2)/((4.d0*c4)+pkm2)))
     3   +((qf5/yoff)-(60.d0*c4r4)+(5.d0*qf3*yoff)
     4   +((5.d0/2.d0)*yyoff*qfq)+((1.d0/4.d0)*q4mp4))*
     5    (datan(pkp/(2.d0*c4)))
     6   +((qf5/yoff)+(60.d0*c4r4)+(5.d0*qf3*yoff)
     7   -((5.d0/2.d0)*yyoff*qfq)+((1.d0/4.d0)*p4mq4))*
     8    (datan(pkm/(2.d0*c4)))


       haitch2=(((qf*c4)/(7.d0*yyoff))*(((19.d0*qf4)/10.d0)
     1  +(((16.d0*qfq)/15.d0)*c4mp2)+(8.d0*c4r4)+(68.d0*yyoff*c4)
     2  +(y4off/2.d0))) +((rc4/(8.d0*y3off))*(((-64.d0/7.d0)*c4r6)
     3  +((48.d0/5.d0)*c4r4*p2mq2)+(4.d0*c4*qp2mpq4)+((yyoff-qfq)**3)))*
     4   (dlog(((4.d0*c4)+pkp2)/((4.d0*c4)+pkm2)))
     5  +(((2.d0*qf7)/(35.d0*y3off))+(2.d0*qf3*yoff)
     6  +((6.d0/5.d0)*yyoff*qfq)+ qf4-(y4off/7.d0)-(16.d0*c4r4))
     7  *(datan(pkp/(2.d0*c4)))+(((2.d0*qf7)/(35.d0*y3off))
     8  +(2.d0*qf3*yoff)-((6.d0/5.d0)*yyoff*qfq)-qf4
     9  +(y4off/7.d0)+(16.d0*c4r4))*(datan(pkm/(2.d0*c4)))


       haitch3=(((qf*c4)/(140.d0*y4off))*((41.d0*y4off)
     1 +(8.d0*yyoff*q2mc4)-(3.d0*((19.d0*qf4)+(64.d0*qfq*c4)
     2 +(80.d0*c4r4)))))+((rc4/y5off)*(((c4/4.d0)*((3.d0*qf4)
     3 -(2.d0*qfq*yyoff)-(y4off)))+((3.d0/5.d0)*c4r4*((3.d0*qfq)+yyoff))
     4 +((12.d0/7.d0)*c4r6)+((3.d0/16.d0)*((qfq-yyoff)**3)))*
     5 (dlog(((4.d0*c4)+pkp2)/((4.d0*c4)+pkm2))))
     6 +((pkp4/(35.d0*y5off))*((4.d0*y3off)-(16.d0*qf*yyoff)
     7 +(12.d0*qfq*yoff)-(3.d0*qf3))*(datan(pkp/(2.d0*c4))))
     8 -((pkm4/(35.d0*y5off))*((4.d0*y3off)+(16.d0*qf*yyoff)
     9 +(12.d0*qfq*yoff)+(3.d0*qf3))*(datan(pkm/(2.d0*c4))))


      go to 31
c
c
c        prepare functions of x, y, and angle
c        -----------------------------------
c
c
   30 if (nt.eq.nnt) go to 50
   31 continue
c
      nqf=32 
      call gseteff (0.d0,qf,nqf,ak,wk)
c
      nnt=nt
      do 45 i=1,nt
      xy2t=xy2*ct(i)
c
c        function  -q^2 (- momentum-transfer-squared)
c        --------------
c        retardation ignored
      deltaq(i,1)=xy2t-xxpyy
c
c        retardation incorporated
      deltaq(i,2)=xy2t-eem12
c
c        function  +k^2 (average-momentum squared)
c        --------------
      deltaq(i,3)=(xy2t+xxpyy)*0.25d0
c
c        function  q^4 (momentum-transfer to the power of 4)
c        -------------
      deltaq(i,4)=deltaq(i,1)*deltaq(i,1)
c
c        function  k^4 (average-momentum to the power of 4)
c        -------------
      deltaq(i,5)=deltaq(i,3)*deltaq(i,3)
c
c        function  +q^2*k^2
c        -----------------
      deltaq(i,6)=-deltaq(i,1)*deltaq(i,3)
c
c        function  (\vec q x \vec k)^2
c        -----------------------------
      deltaq(i,7)=xx*yy*(1.d0-ct(i)*ct(i))
c
c
c        calculate the G-integrals
c        -------------------------
c        -------------------------
c
      gint0(i)=0.d0
      gints(i)=0.d0
      gintss(i)=0.d0
c
      aqq=-deltaq(i,1)
      aq=dsqrt(aqq)
c
      do 35 ik=1,nqf
      ak1=ak(ik)
      ak2=ak1*ak1
      ak3=ak2*ak1
      ak5=ak3*ak2
      ap=(c4+(ak1+yoff)**2)*(c4+(ak1-yoff)**2)
      ap1=dsqrt(ap)
      ap2=dsqrt(ap+aqq*ak2)
      aint=1.d0/ap2*dlog((aq*ak1+ap2)/ap1)
c
      gint0(i)=gint0(i)+ak1*aint*wk(ik)
      gints(i)=gints(i)+ak3*aint*wk(ik)
      gintss(i)=gintss(i)+ak5*aint*wk(ik)
c
   35 continue
c
      gint0(i)=gint0(i)*2.d0/aq
      gints(i)=gints(i)*2.d0/aq
      gintss(i)=gintss(i)*2.d0/aq
c
      c4yy=c4+yyoff
      denom=4.d0*yyoff-aqq
c
      gint1(i)=(gamma0-c4yy*gint0(i)-gints(i))/denom
      gint1s(i)=(3.d0*gamma2+yyoff*gamma3-c4yy*gints(i)-gintss(i))
     1          /denom
      gint2(i)=c4yy*gint1(i)+gints(i)+gint1s(i)
      gint3(i)=(0.5d0*gamma1-2.d0*c4yy*gint1(i)-2.d0*gint1s(i)-gints(i))
     1         /denom
   45 continue
c
c
c
c
c        propagator
c        ----------
c        ----------
c
c
   50 continue
      iprsp=ic(1,im)
      if (iprsp.lt.0) go to 60
      iret=iprsp+1
c
c         propagator for the nn case
      do 55 i=1,nt
   55 aa(i)=wt(i)/(c4-deltaq(i,iret))
      go to 80
c
c
c        "no propagator"
c
   60 do 65 i=1,nt
   65 aa(i)=wt(i)
c
c
   80 continue
c
c
c
c
c
c        cut-offs and functions
c        ----------------------
c        ----------------------
c
c
      mi=4
      mm=5
c
c
 5999 ityp=ic(mi,im)
      if (ityp.eq.0) go to 8000
      if (ityp.le.10) then
      iprspc=ic(mi+1,im)
      iret=iprspc+1
      end if
      go to (100,100,300,9002,500,600,9002,9002,9002,9002,
     1 1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     2 2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     3 3100,3200,3300,3400,3500,3600,3700,3800),ityp
c
c
c
c
c        cut-off of dipole type
c        **********************
c
c
  100 c5=c(mm,im)
      c6=c(mm+1,im)
      nexp=ic(mi+2,im)
c
      do 105 i=1,nt
c
      aaa=c5/(c6-deltaq(i,iret))
c     -------------------------
c
      do 105 ii=1,nexp
  105 aa(i)=aa(i)*aaa
c
c
      mi=mi+3
      mm=mm+2
      go to 5999
c
c
c
c
c        exponential form factor of momentum transfer
c        ********************************************
c
c
  300 c5=c(mm,im)
      c6=c(mm+1,im)
      do 305 i=1,nt
c
      expo=(c5*dabs(deltaq(i,iret)))**c6
c     ----------------------------
c
      if (expo.gt.rrr) expo=rrr
c
      aa(i)=aa(i)*dexp(-expo)
c     ----------------------
c
  305 continue
      mi=mi+2
      mm=mm+2
      go to 5999
c
c
c
c
c        sharp cutoff of x and y
c        ***********************
c
c
  500 c5=c(mm,im)
c
      if (x.gt.c5.or.y.gt.c5) then
c     ----------------------------
      do 505 i=1,nt
  505 aa(i)=0.d0
      end if
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c
c
c        exponential form factor of xx and yy
c        ************************************
c
c
  600 c5=c(mm,im)
      c6=c(mm+1,im)
c
      expo=(c5*xx)**c6+(c5*yy)**c6
c     ----------------------------
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
c     ------------------
c
      do 605 i=1,nt
  605 aa(i)=aa(i)*expexp
      mi=mi+2
      mm=mm+2
      go to 5999
c
c
c
c
c        function +q^2 (momentum-transfer squared)
c        *************
c
c
 1100 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 1105 i=1,nt
 1105 aa(i)=-aa(i)*deltaq(i,1)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function k^2 (average-momentum squared)
c        ************
c
c
 1200 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 1205 i=1,nt
 1205 aa(i)=aa(i)*deltaq(i,3)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,1
c        *******************
c
c
 1300 continue
      c5=c(mm,im)
      do 1305 i=1,nt
      aqq=-deltaq(i,1)
 1305 aa(i)=aa(i)*c5*(2.d0*cb1*c4+cb3*aqq)/(c4+aqq)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,2
c        *******************
c
c
 1400 continue
      c5=c(mm,im)
      do 1405 i=1,nt
      aqq=-deltaq(i,1)
      confac=-4.d0*cb1*c4*(gamma0+gamma1)-(cb3+cb4)*
     1       (aqq*gam013+4.d0*gamma2)
     2       +4.d0*cb4*(2.d0*qf3/3.d0-c4*gamma0)
 1405 aa(i)=aa(i)*c5*confac
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1500 continue
      c5=c(mm,im)
      do 1505 i=1,nt
      aqq=-deltaq(i,1)
      c42qq=2.d0*c4+aqq
      facct=-12.d0*cb1*c4*(2.d0*gamma0-c42qq*gint0(i))
     1       -cb3*(8.d0*qf3-12.d0*c42qq*gamma0-6.d0*aqq*gamma1
     2       +3.d0*c42qq**2*gint0(i))
 1505 aa(i)=aa(i)*c5*facct
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1600 continue
      c5=c(mm,im)
      do 1605 i=1,nt
 1605 aa(i)=aa(i)*c5*cb4*gint2(i)*4.d0
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function for vmed,3
c        *******************
c
c
 1700 continue
      c5=c(mm,im)
      do 1705 i=1,nt
      aqq=-deltaq(i,1)
      c42qq=2.d0*c4+aqq
      facct=2.d0*cb4*(2.d0*gamma0+2.d0*gamma1
     1      -c42qq*(gint0(i)+2.d0*gint1(i)))
 1705 aa(i)=aa(i)*c5*facct
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1800 continue
      c5=c(mm,im)
      do 1805 i=1,nt
      aqq=-deltaq(i,1)
      c42qq=2.d0*c4+aqq
      term1=6.d0*cb3*(2.d0*gamma0+2.d0*gamma1
     1      -c42qq*(gint0(i)+2.d0*gint1(i)))
      term2=24.d0*cb1*c4*(gint0(i)+2.d0*gint1(i))
 1805 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,3
c        *******************
c
c
 1900 continue
      c5=c(mm,im)
      do 1905 i=1,nt
      facct=4.d0*cb4*(gint0(i)+4.d0*gint1(i)+4.d0*gint3(i))
 1905 aa(i)=aa(i)*c5*facct
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,4
c        *******************
c
c
 2000 continue
      c5=c(mm,im)
      do 2005 i=1,nt
 2005 aa(i)=aa(i)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2100 continue
      c5=c(mm,im)
      do 2105 i=1,nt
      aqq=-deltaq(i,1)
 2105 aa(i)=aa(i)*c5*(2.d0*gamma2+(2.d0*yyoff-0.5d0*aqq)*gam013)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2200 continue
      c5=c(mm,im)
      do 2205 i=1,nt
      aqq=-deltaq(i,1)
 2205 aa(i)=aa(i)*c5*(1.d0-2.d0*yyoff/aqq)*gam013
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2300 continue
      c5=c(mm,im)
      do 2305 i=1,nt
      aqq=-deltaq(i,1)
 2305 aa(i)=aa(i)*c5*2.d0/aqq*gam013
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,5
c        *******************
c
c
 2400 continue
      c5=c(mm,im)
      confac=4.d0*qf3-6.d0*c4*gamma0
      do 2405 i=1,nt
 2405 aa(i)=aa(i)*c5*confac
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,6
c        *******************
c
c
 2500 continue
      c5=c(mm,im)
      do 2505 i=1,nt
 2505 aa(i)=aa(i)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c        function for vmed,7
c        *******************
c
c
c 2600 continue
c      c5=c(mm,im)
c      do 2605 i=1,nt
c       aqq=-deltaq(i,1)
c      c42qq=2.d0*c4 + aqq
c      facct= rootc4+c42qq*(1.d0/(2.d0*aq)*
c     1       datan(aq/(2.d0*rootc4)))
c 2605 aa(i)=aa(i)*c5*facct*qf3
c      mi=mi+1
c      mm=mm+1
c      go to 5999

c
c        function for vmed,7
c        *******************
c
c
 2600 continue
      c5=c(mm,im)
      do 2605 i=1,nt
      aqq=-deltaq(i,1)
      aq1 = dsqrt(aqq) 
      rc4= dsqrt(c4)
      c42qq=2.d0*c4 + aqq
      capA=(1.d0/(2.d0*aq1))*datan(aq1/(2.d0*rc4))
      cff= rc4+(c42qq*capA)
 2605 aa(i)=aa(i)*c5*cff*qf3
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c    function for vmed,8
c    *******************
 2700  continue
       c5=c(mm,im)
       do 2705 i=1,nt
       rc4=dsqrt(c4)
       c4r3=c4*rc4
       pkp= yoff+qf
       pkp2=pkp**2
       pkm=yoff-qf
       pkm2=pkm**2
       t1=qf*rc4*(yyoff+11.d0*qfq-4.d0*c4)
       t2a=((qf/yoff)*(qfq+10.d0*c4))-5.d0*c4
       t2b=(7.d0*qfq/4.d0)+(qf*yoff/2.d0)-(yyoff/4.d0)
       t2ab=t2a+t2b 
       t2=pkp2*t2ab*datan(pkp/(2.d0*rc4))
       t3a=((qf/yoff)*(qfq+10.d0*c4))+5.d0*c4
       t3b=-(7.d0*qfq/4.d0)+(qf*yoff/2.d0)+(yyoff/4.d0)
       t3ab=t3a+t3b
       t3=pkm2*t3ab*datan(pkm/(2.d0*rc4))
       t4a=(5.d0*yyoff)-(5.d0*qfq)+(4.d0*c4) 
       t4b=dlog((4.d0*c4+pkp2)/(4.d0*c4+pkm2))
       t4=(c4r3/yoff)*t4a*t4b      
       t1234=t1+t2+t3+t4
 2705  aa(i)=aa(i)*c5*t1234
       mi=mi+1
       mm=mm+1
       go to 5999
c
c
c
c
c
c    function for vmed,8
c    *******************
 2800  continue
       c5=c(mm,im)
       do 2805 i=1,nt
       rc4=dsqrt(c4)
       c4r3=c4*rc4
       pkp= yoff+qf
       pkp2=pkp**2
       pkm=yoff-qf
       pkm2=pkm**2
       t1=qf*rc4*(yyoff+11.d0*qfq-4.d0*c4)
       t2a=((qf/yoff)*(qfq+10.d0*c4))-5.d0*c4
       t2b=(7.d0*qfq/4.d0)+(qf*yoff/2.d0)-(yyoff/4.d0)
       t2ab=t2a+t2b
       t2=pkp2*t2ab*datan(pkp/(2.d0*rc4))
       t3a=((qf/yoff)*(qfq+10.d0*c4))+5.d0*c4
       t3b=-(7.d0*qfq/4.d0)+(qf*yoff/2.d0)+(yyoff/4.d0)
       t3ab=t3a+t3b
       t3=pkm2*t3ab*datan(pkm/(2.d0*rc4))
       t4a=(5.d0*yyoff)-(5.d0*qfq)+(4.d0*c4)
       t4b=dlog((4.d0*c4+pkp2)/(4.d0*c4+pkm2))
       t4=(c4r3/yoff)*t4a*t4b
       t1234=t1+t2+t3+t4
 2805  aa(i)=aa(i)*c5*t1234
       mi=mi+1
       mm=mm+1
       go to 5999
c
c
c
c
c
c        function for vmed,9
c        *******************
c
c
 2900 continue
      c5=c(mm,im)
      do 2905 i=1,nt
      aqq=-deltaq(i,1)
      aq1 = dsqrt(aqq)
      rc4= dsqrt(c4)
      capA=(1.d0/(2.d0*aq1))*datan(aq1/(2.d0*rc4))
 2905 aa(i)=aa(i)*c5*qf3*capA
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,9
c        *******************
c
c
 3000 continue
      c5=c(mm,im)
      do 3005 i=1,nt
      aqq=-deltaq(i,1)
      aq1 = dsqrt(aqq)
      rc4= dsqrt(c4)
      capA=(1.d0/(2.d0*aq1))*datan(aq1/(2.d0*rc4))
 3005 aa(i)=aa(i)*c5*qf3*aqq*capA
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c        function for vmed,9
c        *******************
c
c
 3100 continue
      c5=c(mm,im)
      do 3105 i=1,nt
      aqq=-deltaq(i,1)
      aq1 = dsqrt(aqq)
      rc4= dsqrt(c4)
      r3c4= c4*rc4
      c42qq=2.d0*c4 + aqq
      c44qq=(4.d0*c4)+aqq
      capA=(1.d0/(2.d0*aq1))*datan(aq1/(2.d0*rc4))
      t1=(r3c4/c44qq)-(3.d0*rc4)-(2.d0*c42qq*capA)
 3105 aa(i)=aa(i)*c5*qf3*t1
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c
c
c        function for vmed,10
c        *******************
c
 3200 continue
      c5=c(mm,im)
      do 3205 i=1,nt
       aqq=-deltaq(i,1) 
      
 3205 aa(i)=aa(i)*c5*haitch1
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function for vmed,10
c        *******************
c
 3300 continue
      c5=c(mm,im)
      do 3305 i=1,nt
       aqq=-deltaq(i,1)

 3305 aa(i)=aa(i)*c5*haitch0
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c
c        function for vmed,10
c        *******************
c
 3400 continue
      c5=c(mm,im)
      do 3405 i=1,nt
       aqq=-deltaq(i,1)

 3405 aa(i)=aa(i)*c5*haitch0
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c
c        function for vmed,10
c        *******************
c
 3500 continue
      c5=c(mm,im)
      do 3505 i=1,nt
       aqq=-deltaq(i,1)

 3505 aa(i)=aa(i)*c5*haitch2
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c
c        function for vmed,10
c        *******************
c
 3600 continue
      c5=c(mm,im)
      do 3605 i=1,nt
       aqq=-deltaq(i,1)

 3605 aa(i)=aa(i)*c5*haitch3*(1.d0/aqq)
      mi=mi+1
      mm=mm+1
      go to 5999

c
c
c
c
c        function for vmed,10
c        *******************
c
 3700 continue
      c5=c(mm,im)
      do 3705 i=1,nt
       aqq=-deltaq(i,1)

 3705 aa(i)=aa(i)*c5*haitch3*(1.d0/aqq)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c
c        function for vmed,10
c        *******************
c
 3800 continue
      c5=c(mm,im)
      do 3805 i=1,nt
       aqq=-deltaq(i,1)

 3805 aa(i)=aa(i)*c5*haitch3*(1.d0/aqq)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c
c
c
c
 9002 write (kwrite,19002) ityp


19002 format (1h ////' error in effaa:  cut/fun typ',i10  ,'  does not e
     1xist in this program.'/' execution terminated.'////)
      stop
c
c
c
c
c
 8000 return
      end
      subroutine legpeff (pj,pjm1,x,j)
c
c
c        subroutine legp   computes the legendre polynominals
c
      real*8 pj,pjm1,x,a,b
c
c
c
c        compute legendre polynom for j equals zero
c
c
      if (j.gt.0) go to 1
      pj=1.d0
      pjm1=0.d0
      if (j.lt.0) pj=0.d0
      return
c
c
c
c        compute legendre polynoms for j equals one
c
c
c
    1 pj=x
      pjm1=1.d0
      if (j.eq.1) return
c
c
c
c        compute legendre polynom for j greater or equal two
c
c
c
      do 2 i=2,j
      a=x*pj
      b=a-pjm1
      pjm1=pj
    2 pj=-b/dfloat(i)+b+a
c
c
      return  
      end
      subroutine gseteff (ax,bx,n,z,w)
c
c
c        this code has been obtained from the CERN computer library
c        in the year of the lord 1972.
c
c
      implicit real*8 (a-h,o-z)
c
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c
      dimension     a(273),x(273),ktab(96)
      dimension z(1),w(1)
c
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c
c**** n=2
      data x(1)/0.577350269189626  d0/, a(1)/1.000000000000000  d0/
c**** n=3
      data x(2)/0.774596669241483  d0/, a(2)/0.555555555555556  d0/
      data x(3)/0.000000000000000  d0/, a(3)/0.888888888888889  d0/
c**** n=4
      data x(4)/0.861136311594053  d0/, a(4)/0.347854845137454  d0/
      data x(5)/0.339981043584856  d0/, a(5)/0.652145154862546  d0/
c**** n=5
      data x(6)/0.906179845938664  d0/, a(6)/0.236926885056189  d0/
      data x(7)/0.538469310105683  d0/, a(7)/0.478628670499366  d0/
      data x(8)/0.000000000000000  d0/, a(8)/0.568888888888889  d0/
c**** n=6
      data x(9)/0.932469514203152  d0/, a(9)/0.171324492379170  d0/
      data x(10)/0.661209386466265 d0/, a(10)/0.360761573048139 d0/
      data x(11)/0.238619186083197 d0/, a(11)/0.467913934572691 d0/
c**** n=7
      data x(12)/0.949107912342759 d0/, a(12)/0.129484966168870 d0/
      data x(13)/0.741531185599394 d0/, a(13)/0.279705391489277 d0/
      data x(14)/0.405845151377397 d0/, a(14)/0.381830050505119 d0/
      data x(15)/0.000000000000000 d0/, a(15)/0.417959183673469 d0/
c**** n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c**** n=9
      data x(20)/0.968160239507626 d0/, a(20)/0.081274388361574 d0/
      data x(21)/0.836031107326636 d0/, a(21)/0.180648160694857 d0/
      data x(22)/0.613371432700590 d0/, a(22)/0.260610696402935 d0/
      data x(23)/0.324253423403809 d0/, a(23)/0.312347077040003 d0/
      data x(24)/0.000000000000000 d0/, a(24)/0.330239355001260 d0/
c**** n=10
      data x(25)/0.973906528517172 d0/, a(25)/0.066671344308688 d0/
      data x(26)/0.865063366688985 d0/, a(26)/0.149451349150581 d0/
      data x(27)/0.679409568299024 d0/, a(27)/0.219086362515982 d0/
      data x(28)/0.433395394129247 d0/, a(28)/0.269266719309996 d0/
      data x(29)/0.148874338981631 d0/, a(29)/0.295524224714753 d0/
c**** n=11
      data x(30)/0.978228658146057 d0/, a(30)/0.055668567116174 d0/
      data x(31)/0.887062599768095 d0/, a(31)/0.125580369464905 d0/
      data x(32)/0.730152005574049 d0/, a(32)/0.186290210927734 d0/
      data x(33)/0.519096129206812 d0/, a(33)/0.233193764591990 d0/
      data x(34)/0.269543155952345 d0/, a(34)/0.262804544510247 d0/
      data x(35)/0.000000000000000 d0/, a(35)/0.272925086777901 d0/
c**** n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c**** n=13
      data x(42)/0.984183054718588 d0/, a(42)/0.040484004765316 d0/
      data x(43)/0.917598399222978 d0/, a(43)/0.092121499837728 d0/
      data x(44)/0.801578090733310 d0/, a(44)/0.138873510219787 d0/
      data x(45)/0.642349339440340 d0/, a(45)/0.178145980761946 d0/
      data x(46)/0.448492751036447 d0/, a(46)/0.207816047536889 d0/
      data x(47)/0.230458315955135 d0/, a(47)/0.226283180262897 d0/
      data x(48)/0.000000000000000 d0/, a(48)/0.232551553230874 d0/
c**** n=14
      data x(49)/0.986283808696812 d0/, a(49)/0.035119460331752 d0/
      data x(50)/0.928434883663574 d0/, a(50)/0.080158087159760 d0/
      data x(51)/0.827201315069765 d0/, a(51)/0.121518570687903 d0/
      data x(52)/0.687292904811685 d0/, a(52)/0.157203167158194 d0/
      data x(53)/0.515248636358154 d0/, a(53)/0.185538397477938 d0/
      data x(54)/0.319112368927890 d0/, a(54)/0.205198463721296 d0/
      data x(55)/0.108054948707344 d0/, a(55)/0.215263853463158 d0/
c**** n=15
      data x(56)/0.987992518020485 d0/, a(56)/0.030753241996117 d0/
      data x(57)/0.937273392400706 d0/, a(57)/0.070366047488108 d0/
      data x(58)/0.848206583410427 d0/, a(58)/0.107159220467172 d0/
      data x(59)/0.724417731360170 d0/, a(59)/0.139570677926154 d0/
      data x(60)/0.570972172608539 d0/, a(60)/0.166269205816994 d0/
      data x(61)/0.394151347077563 d0/, a(61)/0.186161000015562 d0/
      data x(62)/0.201194093997435 d0/, a(62)/0.198431485327111 d0/
      data x(63)/0.000000000000000 d0/, a(63)/0.202578241925561 d0/
c**** n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/
c**** n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/
c**** n=24
      data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
      data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
      data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
      data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
      data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
      data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
      data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
      data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
      data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
      data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
      data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
      data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=32
      data x(94)/0.997263861849481 d0/, a(94)/0.007018610009470 d0/
      data x(95)/0.985611511545268 d0/, a(95)/0.016274394730905 d0/
      data x(96)/0.964762255587506 d0/, a(96)/0.025392065309262 d0/
      data x(97)/0.934906075937739 d0/, a(97)/0.034273862913021 d0/
      data x(98)/0.896321155766052 d0/, a(98)/0.042835898022226 d0/
      data x(99)/0.849367613732569 d0/, a(99)/0.050998059262376 d0/
      data x(100)/0.794483795967942d0/, a(100)/0.058684093478535d0/
      data x(101)/0.732182118740289d0/, a(101)/0.065822222776361d0/
      data x(102)/0.663044266930215d0/, a(102)/0.072345794108848d0/
      data x(103)/0.587715757240762d0/, a(103)/0.078193895787070d0/
      data x(104)/0.506899908932229d0/, a(104)/0.083311924226946d0/
      data x(105)/0.421351276130635d0/, a(105)/0.087652093004403d0/
      data x(106)/0.331868602282127d0/, a(106)/0.091173878695763d0/
      data x(107)/0.239287362252137d0/, a(107)/0.093844399080804d0/
      data x(108)/0.144471961582796d0/, a(108)/0.095638720079274d0/
      data x(109)/0.048307665687738d0/, a(109)/0.096540088514727d0/
c**** n=40
      data x(110)/0.998237709710559d0/, a(110)/0.004521277098533d0/
      data x(111)/0.990726238699457d0/, a(111)/0.010498284531152d0/
      data x(112)/0.977259949983774d0/, a(112)/0.016421058381907d0/
      data x(113)/0.957916819213791d0/, a(113)/0.022245849194166d0/
      data x(114)/0.932812808278676d0/, a(114)/0.027937006980023d0/
      data x(115)/0.902098806968874d0/, a(115)/0.033460195282547d0/
      data x(116)/0.865959503212259d0/, a(116)/0.038782167974472d0/
      data x(117)/0.824612230833311d0/, a(117)/0.043870908185673d0/
      data x(118)/0.778305651426519d0/, a(118)/0.048695807635072d0/
      data x(119)/0.727318255189927d0/, a(119)/0.053227846983936d0/
      data x(120)/0.671956684614179d0/, a(120)/0.057439769099391d0/
      data x(121)/0.612553889667980d0/, a(121)/0.061306242492928d0/
      data x(122)/0.549467125095128d0/, a(122)/0.064804013456601d0/
      data x(123)/0.483075801686178d0/, a(123)/0.067912045815233d0/
      data x(124)/0.413779204371605d0/, a(124)/0.070611647391286d0/
      data x(125)/0.341994090825758d0/, a(125)/0.072886582395804d0/
      data x(126)/0.268152185007253d0/, a(126)/0.074723169057968d0/
      data x(127)/0.192697580701371d0/, a(127)/0.076110361900626d0/
      data x(128)/0.116084070675255d0/, a(128)/0.077039818164247d0/
      data x(129)/0.038772417506050d0/, a(129)/0.077505947978424d0/
c**** n=48
      data x(130)/0.998771007252426d0/, a(130)/0.003153346052305d0/
      data x(131)/0.993530172266350d0/, a(131)/0.007327553901276d0/
      data x(132)/0.984124583722826d0/, a(132)/0.011477234579234d0/
      data x(133)/0.970591592546247d0/, a(133)/0.015579315722943d0/
      data x(134)/0.952987703160430d0/, a(134)/0.019616160457355d0/
      data x(135)/0.931386690706554d0/, a(135)/0.023570760839324d0/
      data x(136)/0.905879136715569d0/, a(136)/0.027426509708356d0/
      data x(137)/0.876572020274247d0/, a(137)/0.031167227832798d0/
      data x(138)/0.843588261624393d0/, a(138)/0.034777222564770d0/
      data x(139)/0.807066204029442d0/, a(139)/0.038241351065830d0/
      data x(140)/0.767159032515740d0/, a(140)/0.041545082943464d0/
      data x(141)/0.724034130923814d0/, a(141)/0.044674560856694d0/
      data x(142)/0.677872379632663d0/, a(142)/0.047616658492490d0/
      data x(143)/0.628867396776513d0/, a(143)/0.050359035553854d0/
      data x(144)/0.577224726083972d0/, a(144)/0.052890189485193d0/
      data x(145)/0.523160974722233d0/, a(145)/0.055199503699984d0/
      data x(146)/0.466902904750958d0/, a(146)/0.057277292100403d0/
      data x(147)/0.408686481990716d0/, a(147)/0.059114839698395d0/
      data x(148)/0.348755886292160d0/, a(148)/0.060704439165893d0/
      data x(149)/0.287362487355455d0/, a(149)/0.062039423159892d0/
      data x(150)/0.224763790394689d0/, a(150)/0.063114192286254d0/
      data x(151)/0.161222356068891d0/, a(151)/0.063924238584648d0/
      data x(152)/0.097004699209462d0/, a(152)/0.064466164435950d0/
      data x(153)/0.032380170962869d0/, a(153)/0.064737696812683d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/
c**** n=80
      data x(186)/0.999553822651630d0/, a(186)/0.001144950003186d0/
      data x(187)/0.997649864398237d0/, a(187)/0.002663533589512d0/
      data x(188)/0.994227540965688d0/, a(188)/0.004180313124694d0/
      data x(189)/0.989291302499755d0/, a(189)/0.005690922451403d0/
      data x(190)/0.982848572738629d0/, a(190)/0.007192904768117d0/
      data x(191)/0.974909140585727d0/, a(191)/0.008683945269260d0/
      data x(192)/0.965485089043799d0/, a(192)/0.010161766041103d0/
      data x(193)/0.954590766343634d0/, a(193)/0.011624114120797d0/
      data x(194)/0.942242761309872d0/, a(194)/0.013068761592401d0/
      data x(195)/0.928459877172445d0/, a(195)/0.014493508040509d0/
      data x(196)/0.913263102571757d0/, a(196)/0.015896183583725d0/
      data x(197)/0.896675579438770d0/, a(197)/0.017274652056269d0/
      data x(198)/0.878722567678213d0/, a(198)/0.018626814208299d0/
      data x(199)/0.859431406663111d0/, a(199)/0.019950610878141d0/
      data x(200)/0.838831473580255d0/, a(200)/0.021244026115782d0/
      data x(201)/0.816954138681463d0/, a(201)/0.022505090246332d0/
      data x(202)/0.793832717504605d0/, a(202)/0.023731882865930d0/
      data x(203)/0.769502420135041d0/, a(203)/0.024922535764115d0/
      data x(204)/0.744000297583597d0/, a(204)/0.026075235767565d0/
      data x(205)/0.717365185362099d0/, a(205)/0.027188227500486d0/
      data x(206)/0.689637644342027d0/, a(206)/0.028259816057276d0/
      data x(207)/0.660859898986119d0/, a(207)/0.029288369583267d0/
      data x(208)/0.631075773046871d0/, a(208)/0.030272321759557d0/
      data x(209)/0.600330622829751d0/, a(209)/0.031210174188114d0/
      data x(210)/0.568671268122709d0/, a(210)/0.032100498673487d0/
      data x(211)/0.536145920897131d0/, a(211)/0.032941939397645d0/
      data x(212)/0.502804111888784d0/, a(212)/0.033733214984611d0/
      data x(213)/0.468696615170544d0/, a(213)/0.034473120451753d0/
      data x(214)/0.433875370831756d0/, a(214)/0.035160529044747d0/
      data x(215)/0.398393405881969d0/, a(215)/0.035794393953416d0/
      data x(216)/0.362304753499487d0/, a(216)/0.036373749905835d0/
      data x(217)/0.325664370747701d0/, a(217)/0.036897714638276d0/
      data x(218)/0.288528054884511d0/, a(218)/0.037365490238730d0/
      data x(219)/0.250952358392272d0/, a(219)/0.037776364362001d0/
      data x(220)/0.212994502857666d0/, a(220)/0.038129711314477d0/
      data x(221)/0.174712291832646d0/, a(221)/0.038424993006959d0/
      data x(222)/0.136164022809143d0/, a(222)/0.038661759774076d0/
      data x(223)/0.097408398441584d0/, a(223)/0.038839651059051d0/
      data x(224)/0.058504437152420d0/, a(224)/0.038958395962769d0/
      data x(225)/0.019511383256793d0/, a(225)/0.039017813656306d0/
c**** n=96
      data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
      data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
      data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
      data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
      data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
      data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
      data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
      data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
      data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
      data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
      data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
      data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
      data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
      data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
      data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
      data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
      data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
      data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
      data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
      data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
      data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
      data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
      data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
      data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
      data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
      data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
      data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
      data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
      data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
      data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
      data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
      data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
      data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
      data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
      data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/
      data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/
      data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/
      data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/
      data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/
      data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/
      data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/
      data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/
      data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/
      data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
      data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
      data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
      data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
      data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
c
c
c-----test n
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c
  100 zn=n
      write(6,200) zn
  200 format(/////' error in gset. n has the non-permissible value',
     1e11.3/' execution terminated.')
      stop
      end
c********************************************************************
c name:    dminv
c        programmbibliothek rhrz bonn        28/11/78       dminv
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c invert a matrix
c
c usage:   call dminv (a,n,d,l,m)
c
c parameters:
c
c a:       input matrix, destroyed in computation and replaced by
c          resultant inverse.
c          double precision required.
c
c n:       order of matrix a
c
c d:       resultant determinant
c          double precision required.
c
c l:       work vector of length n
c
c m:       work vector of length n
c
c remarks: matrix a must be a general matrix
c
c method:
c
c the standard gauss-jordan method is used. the determinant
c is also calculated. a determinant of zero indicates that
c the matrix is singular.
c
c programs required:
c          none
c
c author:  ibm, ssp iii
c
c**********************************************************************
      subroutine dminveff (a,n,d,l,m)
      implicit real*8 (a-h,o-z)
      dimension a(1),l(1),m(1)

c
c
c        search for largest element
c
      d=1.d0
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
   10 if (dabs(biga)-dabs(a(ij)))  15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue
c
c        interchange rows
c
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
c
c        interchange columns
c
   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
c
c        divide column by minus pivot (value of pivot element is
c        contained in biga)
c
   45 if(biga) 48,46,48
   46 d=0.d0
      return
   48 do 55 i=1,n
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
c
c        reduce matrix
c
      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
c
c        divide row by pivot
c
      kj=k-n
      do 75 j=1,n
      kj=kj+n
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
c
c        product of pivots
c
      d=d*biga
c
c        replace pivot by reciprocal
c
      a(kk)=1.d0/biga
   80 continue
c
c        final row and column interchange
c
      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      do 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      if(j-k) 100,100,125
  125 ki=k-n
      do 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 return
      end
c**************** this is the end of the program eff3nfn3lo **********************
