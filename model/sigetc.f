      real function sigocn(t,s)
c
c --- sigma, dsigma/dt, and dsigma/ds as functions of temp and salinity
c
      implicit none
      real t,s,p,bulk
c
      include 'state_eqn.h'
c
ccc   sig=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))*1.e-3                ! cgs
c7    sigocn=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))                !  SI
c9    sigocn=(c1+s*(c3+s*(c8+c9*t))+t*(c2+c5*s+t*(c4+c7*s+c6*t)))       !  SI
c
c  New scheme using Jackett and McDougall's revised equation which
c  computes rho correctly when the input is potential temperature.
c  input is whatever in rho0 (example: 1000.), s is salinity in
c  PSU, t is potential temperature in deg C, and p is pressure in
c  bars.
      p=200.                              ! bars
      bulk = 1.965933e4 + t*(1.444304e2
     & + t*(-1.706103e0 + t*(9.648704e-3
     & + t*(-4.190253e-5))))
     & + s*(5.284855e1 + t*(-3.101089e-1
     & + t*(6.283263e-3 + t*(-5.084188e-5))))
     & + sqrt(s**3)*(3.886640e-1
     & + t*(9.085835e-3 + t*(-4.619924e-4)))
     & + p*(3.186519e0 + t*(2.212276e-2
     & + t*(-2.984642e-4 + t*(1.956415e-6))))
     & + p*s*(6.704388e-3
     & + t*(-1.847318e-4 + t*(2.059331e-7)))
     & + p*sqrt(s**3)
     & *1.480266e-4 + p*p*(2.102898e-4
     & + t*(-1.202016e-5 + t*(1.394680e-7)))
     & + p*p*s
     & *(-2.040237e-6
     & + t*(6.128773e-8 + t*(6.207323e-10)))
      sigocn = ( 999.842594e0
     &           + t*( 6.793952e-2
     &           + t*(-9.095290e-3
     &           + t*( 1.001685e-4
     &           + t*(-1.120083e-6
     &           + t*( 6.536332e-9)))))
     &           + s*( 0.824493e0
     &           + t *(-4.08990e-3
     &           + t *( 7.64380e-5
     &           + t *(-8.24670e-7
     &           + t *( 5.38750e-9)))))
     &           + sqrt(s**3)*(-5.72466e-3
     &           + t*( 1.02270e-4
     &           + t*(-1.65460e-6)))
     &           + 4.8314e-4*s**2)
     &           / (1.0 - p/bulk) - 1000.
      return
      end
c
c
      real function dsigdt(t,s)
      implicit none
      real t,s
c
      include 'state_eqn.h'
c
      dsigdt=(c2+s*(c5+c9*s)+2.*t*(c4+c7*s+1.5*c6*t))           !  SI c9
      return
      end
c
c
      real function dsigds(t,s)
      implicit none
      real t,s
c
      include 'state_eqn.h'
c
c     dsigds=(c3+t*(c5+t*c7))                                        ! SI c7
      dsigds=c3+2*s*(c8+c9*t)+t*(c5+t*c7)                       ! SI c9
      return
      end
c
c
      real function tofsig(sigm,salin)
c
c --- temp (deg c) as a function of sigma and salinity (this routine mimics
c --- the statement function of the same name in state_eqn.h.)
c
      implicit none
      real sigm,salin,s,sqq,a0,a1,a2,cubr,cubq,cuban,cubrl,cubim,athird
      parameter (athird=1./3.)
c
      include 'state_eqn.h'
c
      a0=(c1+salin*(c3+c8*salin))/c6
      a1=(c2+salin*(c5+c9*salin))/c6
      a2=(c4+c7*salin)/c6
c
      cubq=athird*a1-(athird*a2)**2
ccc   cubr=athird*(.5*a1*a2-1.5*(a0-1.e3*sigm/c6))                ! cgs
      cubr=athird*(.5*a1*a2-1.5*(a0-     sigm/c6))                !  SI
     .   -(athird*a2)**3
c
c --- if q**3+r**2>0, water is too dense to yield real root at given
c --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
c --- lowering sigma until a double real root is obtained.
c
      cuban=athird*atan2(sqrt(max(0.,-(cubq**3+cubr**2))),cubr)
      sqq=sqrt(-cubq)
      cubrl=sqq*cos(cuban)
      cubim=sqq*sin(cuban)
      tofsig=-cubrl+sqrt(3.)*cubim-athird*a2
ccc      if (abs(sig(tofsig,salin)-sigm).gt.1.e-6) write (lp,100)
ccc     .   tofsig,salin,sigm,sig(tofsig,salin)
 100  format ('tofsig,sal,old/new sig =',2f9.3,3p,2f9.3)
      return
      end
c
c
      real function sofsig(sigma,tem)
      implicit none
      real sigma,tem,bb,aa,cc
c
      include 'state_eqn.h'
c
c7    sofsig=(sigma                                          ! SI
c7   .  -c1-tem*(c2+tem*(c4+c6*tem)))/(c3+tem*(c5+c7*tem))
      aa=c8+c9*tem
      bb=c3+c5*tem+c7*tem*tem
      cc=c1+c2*tem+c4*tem*tem+c6*tem*tem*tem-sigma
      sofsig=(-bb+sqrt(bb*bb-4.*aa*cc))/(2.*aa)
      return
      end
c
c
      real function qsatur(t)
c
c --- saturation specific humidity (lowe, j.appl.met., 16, 100-103, 1976)
c
      implicit none
      real, intent(in) :: t
      qsatur=.622e-3*(6.107799961e+00+t*(4.436518521e-01
     .            +t*(1.428945805e-02+t*(2.650648471e-04
     .            +t*(3.031240396e-06+t*(2.034080948e-08
     .            +t* 6.136820929e-11))))))
      return
      end
c
      real function kappaf(t,s,prs,th,profil)
c
c --- compressibility coefficient kappa^(theta) from sun et al. (1999)
c
      USE HYCOM_SCALARS, only : thref
      USE HYCOM_ARRAYS_GLOB
      implicit none
      real t,s,prs,tdif,sdif,th,kappa1,kappa2,profil
      integer loc1,loc2
c
      include 'state_eqn.h'
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      kappaf=0.                              ! no thermobaricity
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      loc1=profil
      if (loc1.lt.1 .or. loc1.gt.4) then
        write (*,'(a,f7.3)') 'illegal profile value in kappaf:',profil
        stop '(kappaf)'
      end if
      sdif=max(30.,min(38.,s))-refsal(loc1)
      tdif=max(-2.,min(32.,t))-reftem(loc1)
      kappa1=sclkap * (tdif*(qt(loc1)+tdif*(qtt(loc1)
     .  +tdif*qttt(loc1))+.5*(prs+pref)*(qpt(loc1)
     .  +sdif*qpst(loc1)+tdif*qptt(loc1)))
     .  +sdif*(qs(loc1)+tdif*qst(loc1)))*(prs-pref)
      kappa1=(1./thref+th)*(exp(kappa1)-1.)
c
      loc2=min(4,loc1+1)
      sdif=max(30.,min(38.,s))-refsal(loc2)
      tdif=max(-2.,min(32.,t))-reftem(loc2)
      kappa2=sclkap * (tdif*(qt(loc2)+tdif*(qtt(loc2)
     .  +tdif*qttt(loc2))+.5*(prs+pref)*(qpt(loc2)
     .  +sdif*qpst(loc2)+tdif*qptt(loc2)))
     .  +sdif*(qs(loc2)+tdif*qst(loc2)))*(prs-pref)
      kappa2=(1./thref+th)*(exp(kappa2)-1.)
c
      kappaf=kappa1+(kappa2-kappa1)*(profil-float(loc1))
c
      return
      end
c
      real function sigloc(t,s,prs)
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
c
      implicit none
      include 'state_eqn.h'
c
      real t,s,prs,c1p,c2p,c3p,c4p,c5p,c6p,c7p
      c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
      c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
      c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
      c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
      c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
c
      sigloc=c1p(prs)+c3p(prs)*s+
     &    t*(c2p(prs)+c5p(prs)*s+t*(c4p(prs)+c7p(prs)*s+c6p(prs)*t))
      return
      end
c
c
      real function dsiglocdt(t,s,prs)
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
c
      implicit none
      include 'state_eqn.h'
c
      real t,s,prs,c2p,c4p,c5p,c6p,c7p
ccc   c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
      c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
ccc   c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
      c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
      c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
c
      dsiglocdt=c2p(prs)+c5p(prs)*s+
     &    2.*t*(c4p(prs)+c7p(prs)*s+1.5*c6p(prs)*t)
      return
      end
c
c
      real function dsiglocds(t,s,prs)
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
c
      implicit none
      include 'state_eqn.h'
c
      real t,s,prs,c3p,c5p,c7p
ccc   c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
ccc   c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
      c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
ccc   c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
ccc   c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
c
      dsiglocds=c3p(prs)+t*(c5p(prs)+t*c7p(prs))
      return
      end
c
      subroutine cpy_p(field)
c
c --- exchange information across bering strait seam
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : ipacn,ipacs,jpac,iatln,iatls,jatl,beropn

      implicit none
c
      real field(idm,jdm),sign
c
c --- exchange p-point values (half grid size away from seam)
      if (beropn) then
        field(iatls,jatl)=field(ipacs,jpac)
        field(ipacn,jpac)=field(iatln,jatl)
      endif
      return
      end

      subroutine cpy_p_par(field)
c
c --- exchange information across bering strait seam
c
      USE DOMAIN_DECOMP_1D, only : send_to_j, recv_from_j
      USE HYCOM_DIM, only : ogrid,J_0,J_1,J_0H,J_1H,I_0H,I_1H
      USE HYCOM_SCALARS, only : ipacn,ipacs,jpac,iatln,iatls,jatl,beropn

      implicit none
c
      real field(I_0H:I_1H,J_0H:J_1H)
      real a(1),b(1)
c

c --- exchange p-point values (half grid size away from seam)
      if (beropn) then

        if ( jpac>=J_0 .and. jpac<=J_1 .and.
     &       jatl>=J_0 .and. jatl<=J_1 ) then ! all local
          field(iatls,jatl)=field(ipacs,jpac)
          field(ipacn,jpac)=field(iatln,jatl)
        else                    ! jpac, jatl are on different threads

          if ( jpac>=J_0 .and. jpac<=J_1 ) then
            a(1) = field(ipacs,jpac)
            call send_to_j(ogrid,a,jatl,1)
          endif
          if ( jatl>=J_0 .and. jatl<=J_1 ) then
            b(1) = field(iatln,jatl)
            call send_to_j(ogrid,b,jpac,2)
          endif

          if ( jatl>=J_0 .and. jatl<=J_1 ) then
            call recv_from_j(ogrid,a,jpac,1)
            field(iatls,jatl)=a(1)
          endif
          if ( jpac>=J_0 .and. jpac<=J_1 ) then
            call recv_from_j(ogrid,b,jatl,2)
            field(ipacn,jpac)=b(1)
          endif

        endif
          !field(iatls,jatl)=field(ipacs,jpac)
          !field(ipacn,jpac)=field(iatln,jatl)
      endif

      return
      end
c
      subroutine cpy_mJpacJatL(field)
c
c --- copy information across bering strait seam
c --- from jpac to jatl
c ---     field(iatls,jatl) = - field(ipacn,jpac)
c
      USE DOMAIN_DECOMP_1D, only : send_to_j, recv_from_j
      USE HYCOM_DIM, only : ogrid,J_0,J_1,J_0H,J_1H,I_0H,I_1H
      USE HYCOM_SCALARS, only : ipacn,ipacs,jpac,iatln,iatls,jatl,beropn

      implicit none
c
      real field(I_0H:I_1H,J_0H:J_1H)
      real a(1),b(1)
c
c --- exchange p-point values (half grid size away from seam)
      if ( jpac>=J_0 .and. jpac<=J_1 .and.
     &     jatl>=J_0 .and. jatl<=J_1 ) then ! all local
        field(iatls,jatl) = field(ipacn,jpac)
      else                      ! jpac, jatl are on different threads

        if ( jpac>=J_0 .and. jpac<=J_1 ) then
          a(1) = field(ipacn,jpac)
          call send_to_j(ogrid,a,jatl,1)
        endif
        if ( jatl>=J_0 .and. jatl<=J_1 ) then
          call recv_from_j(ogrid,a,jpac,1)
          field(iatls,jatl)=a(1)
        endif

      endif

      return
      end subroutine cpy_mJpacJatL
c
c
c> Revision history:
c>
c> Jan. 2001 - modified kappaf to allow nonzero temp.ref.value
c> Nov. 2002 - corrected error in kappaf
c> July 2005 - added sig,dsigdt,dsigds functions for in-situ density ('sigloc')
c
c
      real function hyc_pechg1(delp,sig,nunit)
c
c --- calculate change in ocean's available potential energy due to some
c --- process 'X'. call pechg1  b e f o r e  , pechg2  a f t e r  process X.
c --- input: 3-d arrays of  h y b r i d  layer thickness and density anomaly.
c --- results from pechg1 are stored in 'nunit' for later use by pechg2.
c --- use different values of 'nunit' for nested APE process diagnostics.
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : theta,onem,flnmovt,lp,g
      USE HYCOM_ARRAYS_GLOB
      implicit none
      integer i,j,k,l,n
c
      real slithk
      integer nscli,nunit,lgth
      parameter (slithk=1.,nscli=7000./slithk)
      real delp(idm,jdm,kdm),sig(idm,jdm,kdm),
     .     pbfore(idm,jdm,kdm+1),pafter(idm,jdm,kdm+1),
     .     praw(idm,jdm,kdm+1),unused(idm,jdm,kdm)
      real weight(kdm),sliwgt(nscli),slitop,
     .     slibot,slisum,prefbe(kdm),prefaf(kdm),varian(0:kdm),
     .     weightj(jdm),sliwgtj(jdm),varianj(jdm)
      character flnm*60
      data varian(0),varian(kdm)/0.,0./
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
 11   praw(i,j,1)=0.
      do 10 k=1,kk
      do 10 i=ifp(j,l),ilp(j,l)
 10   praw(i,j,k+1)=praw(i,j,k)+delp(i,j,k)
c$OMP END PARALLEL DO
c
c --- transform pressure to isopycnic interface pressure
c
      call reflux(uflx,vflx,sig,praw,
     .            unused,unused,unused,pbfore,theta,kdm,kdm)
c
c --- in preparation for determining the flat reference state (i.e., the
c --- unavailable pot. energy), find total mass above each isopycnic interface
c
      do 1 k=1,kk
      weight(k)=0.
c
c$OMP PARALLEL DO
      do 2 j=1,jj
      weightj(j)=0.
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
 2    weightj(j)=weightj(j)+pbfore(i,j,k+1)*scp2(i,j)
c$OMP END PARALLEL DO
c
      do 1 j=1,jj
 1    weight(k)=weight(k)+weightj(j)
c
c --- divide basin into shallow horizontal slices and find mass of each slice
c
      do 3 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      sliwgt(n)=0.
c$OMP PARALLEL DO
      do 4 j=1,jj
      sliwgtj(j)=0.
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    sliwgtj(j)=sliwgtj(j)+scp2(i,j)*(min(pbfore(i,j,kk+1),slibot)-
     .                                 min(pbfore(i,j,kk+1),slitop))
c$OMP END PARALLEL DO
      do 3 j=1,jj
 3    sliwgt(n)=sliwgt(n)+sliwgtj(j)
c
      do 5 k=1,kk-1
c --- add slices vertically until sum exceeds combined mass of layers 1...k.
c --- this tells us where bottom of layer k will be when flattened
c
      slisum=0.
      do 6 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      slisum=slisum+sliwgt(n)
      if (slisum.ge.weight(k)) go to 7
 6    continue
      write (*,*) 'k =',k,'  error: slisum < weight',slisum,weight(k)
      go to 5
c
c --- interpolate among slices to get precise depth of flattened interface
c
 7    prefbe(k)=(slibot*(weight(k)-slisum+sliwgt(n))-
     .           slitop*(weight(k)-slisum          ))/sliwgt(n)
ccc      write (*,*) 'k =',k,'  flattened interface:',prefbe(k)/onem
 5    continue
c
c --- save results for later use by pechg2
c
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 13
 14   continue
      write (lp,*) 'ape --  cannot find slash in',flnmovt
      stop
 13   write (flnm,'(a,i2.2)') flnmovt(1:lgth)//'ape.',nunit
      open (unit=nunit,file=flnm,form='unformatted',status='unknown')
      write (nunit) pbfore,prefbe
      close (nunit)
c
c --- determine variance of interface pressure relative to flat reference state
c
      do 12 k=1,kk-1
      varian(k)=0.
c
c$OMP PARALLEL DO
      do 9 j=1,jj
      varianj(j)=0.
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    varianj(j)=varianj(j)+scp2(i,j)*
     .  (pbfore(i,j,k+1)-min(pbfore(i,j,kk+1),prefbe(k)))**2
c$OMP END PARALLEL DO
c
      do 12 j=1,jj
 12   varian(k)=varian(k)+varianj(j)
c
      hyc_pechg1=0.
      do 8 k=1,kk
 8    hyc_pechg1=hyc_pechg1+
     .  (varian(k)-varian(k-1))/(1000.+theta(k))
c
c --- report result in units of joules (kg m^2/sec^2)
      hyc_pechg1=.5*hyc_pechg1/g
      return
      end
c
c
      real function hyc_pechg2(delp,sig,nunit)
c
c --- calculate change in ocean's available potential energy due to some
c --- process 'X'. call pechg1  b e f o r e  , pechg2  a f t e r  process X.
c --- input: 3-d arrays of  h y b r i d  layer thickness and density anomaly.
c --- results from pechg1 representing 'before' state are read from 'nunit'.
c --- use different values of 'nunit' for nested APE process diagnostics.
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : theta,onem,flnmovt,lp,g,delt1
      USE HYCOM_ARRAYS_GLOB
      implicit none
      integer i,j,k,l,n
c
      real slithk
      integer nscli,nunit,lgth
      parameter (slithk=1.,nscli=7000./slithk)
      real delp(idm,jdm,kdm),sig(idm,jdm,kdm),
     .     pbfore(idm,jdm,kdm+1),pafter(idm,jdm,kdm+1),
     .     praw(idm,jdm,kdm+1),unused(idm,jdm,kdm)
      real weight(kdm),sliwgt(nscli),slitop,
     .     slibot,slisum,prefbe(kdm),prefaf(kdm),varian(0:kdm),
     .     weightj(jdm),sliwgtj(jdm),varianj(jdm)
      character flnm*60
      data varian(0),varian(kdm)/0.,0./
      common /pechng/ pbfore,pafter,prefbe,prefaf
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
 11   praw(i,j,1)=0.
      do 10 k=1,kk
      do 10 i=ifp(j,l),ilp(j,l)
 10   praw(i,j,k+1)=praw(i,j,k)+delp(i,j,k)
c$OMP END PARALLEL DO
c
c --- transform pressure to isopycnic interface pressure
c
      call reflux(uflx,vflx,sig,praw,
     .            unused,unused,unused,pafter,theta,kdm,kdm)
c
c --- in preparation for determining the flat reference state (i.e., the
c --- unavailable pot. energy), find total mass above each isopycnic interface
c
      do 1 k=1,kk
      weight(k)=0.
c
c$OMP PARALLEL DO
      do 2 j=1,jj
      weightj(j)=0.
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
 2    weightj(j)=weightj(j)+pafter(i,j,k+1)*scp2(i,j)
c$OMP END PARALLEL DO
c
      do 1 j=1,jj
 1    weight(k)=weight(k)+weightj(j)
c
c --- divide basin into shallow horizontal slices and find mass of each slice
c
      do 3 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      sliwgt(n)=0.
c$OMP PARALLEL DO
      do 4 j=1,jj
      sliwgtj(j)=0.
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    sliwgtj(j)=sliwgtj(j)+scp2(i,j)*(min(pafter(i,j,kk+1),slibot)-
     .                                 min(pafter(i,j,kk+1),slitop))
c$OMP END PARALLEL DO
      do 3 j=1,jj
 3    sliwgt(n)=sliwgt(n)+sliwgtj(j)
c
      do 5 k=1,kk-1
c --- add slices vertically until sum exceeds combined mass of layers 1...k.
c --- this tells us where bottom of layer k will be when flattened
c
      slisum=0.
      do 6 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      slisum=slisum+sliwgt(n)
      if (slisum.ge.weight(k)) go to 7
 6    continue
      write (*,*) 'k =',k,'  error: slisum < weight',slisum,weight(k)
      go to 5
c
c --- interpolate among slices to get precise depth of flattened interface
c
 7    prefaf(k)=(slibot*(weight(k)-slisum+sliwgt(n))-
     .           slitop*(weight(k)-slisum          ))/sliwgt(n)
ccc      write (*,*) 'k =',k,'  flattened interface:',prefaf(k)/onem
 5    continue
c
c --- read previously created file representing 'before' state
c
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 13
 14   continue
      write (lp,*) 'ape --  cannot find slash in',flnmovt
      stop
 13   write (flnm,'(a,i2.2)') flnmovt(1:lgth)//'ape.',nunit
      open (unit=nunit,file=flnm,form='unformatted',status='old')
      read (nunit) pbfore,prefbe
      rewind (nunit)
c
c --- write out 'after' state (for potential future use as new 'before' state)
c
      write (nunit) pafter,prefaf
      close (nunit)
c
c --- determine 'before/after' difference of interface pressure variances
c
      do 12 k=1,kk-1
      varian(k)=0.
c
c$OMP PARALLEL DO
      do 9 j=1,jj
      varianj(j)=0.
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    varianj(j)=varianj(j)+scp2(i,j)*
     . ((pafter(i,j,k+1)-min(pafter(i,j,kk+1),prefaf(k)))**2
     . -(pbfore(i,j,k+1)-min(pbfore(i,j,kk+1),prefbe(k)))**2)
c$OMP END PARALLEL DO
c
      do 12 j=1,jj
 12   varian(k)=varian(k)+varianj(j)
c
      hyc_pechg2=0.
      do 8 k=1,kk
 8    hyc_pechg2=hyc_pechg2
     .  +(varian(k)-varian(k-1))/(1000.+theta(k))
c
c --- report result in units of watts (kg m^2/sec^3)
      hyc_pechg2=.5*hyc_pechg2/(g*delt1)
      return
      end
c
c
c> Revision history:
c>
c> Dec. 2004 - fixed bug in loop 9 (excluded interfaces on shallow bottom)
c
      subroutine totals(dp1,field1,dp2,field2,text)
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : lp
      USE HYCOM_ARRAYS_GLOB
      implicit none
c
c --- compute volume integral of 2 fields (field1,field2), each associated
c --- with its own layer thickness field (dp1,dp2)
c
      integer i,j,k,l
c
      real dp1(idm,jdm,kdm),dp2(idm,jdm,kdm),field1(idm,jdm,kdm),
     .     field2(idm,jdm,kdm),sum1j(jdm),sum2j(jdm),sum1,sum2
      character text*(*)
c
c$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
      sum1j(j)=0.
      sum2j(j)=0.
      do 1 k=1,kk
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
      sum1j(j)=sum1j(j)+dp1(i,j,k)*field1(i,j,k)*scp2(i,j)
 1    sum2j(j)=sum2j(j)+dp2(i,j,k)*field2(i,j,k)*scp2(i,j)
c$OMP END PARALLEL DO
c
      sum1=0.
      sum2=0.
c$OMP PARALLEL DO REDUCTION(+:sum1,sum2) SCHEDULE(STATIC,jchunk)
      do 2 j=1,jj
      sum1=sum1+sum1j(j)
 2    sum2=sum2+sum2j(j)
c$OMP END PARALLEL DO
c
      write (lp,'(a,1p,2e19.9)') text,sum1,sum2
      return
      end
c
      subroutine refinp(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- add grid rows near equator to enhance merid. resolution
c --- eqcrs  -- i index of equator in unexpanded grid
c --- eqfin  -- i index of equator in expanded grid
c --- fieldc -- input field containing ii-2*(eqfin-eqcrs) rows
c --- fieldf -- output field containing ii rows
c
c --- use this entry to operate on -p- points
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : lp
      implicit none
      integer i,j
c
      integer numcrs,numfin,newrows,icrs,ifin,ieqcrs,ieqfin,
     .        ipcrs(idm,jdm)
      real fieldc(idm,jdm),fieldf(idm,jdm),
     .     eqcrs,eqfin,wgt,valu0,valu1,x,sumc,sumf
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=1)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0.,1./
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=7)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.25, 0.55, 0.95, 1.50, 2.20, 3.05, 4.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=9)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.400, 0.825, 1.300, 1.850, 2.500, 3.250, 4.100,
ccc     .               5.025, 6.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=11)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.60, 1.20, 1.80, 2.40, 3.02, 3.68, 4.40, 5.20,
ccc     .               6.08, 7.02, 8.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=12)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.375, 0.750, 1.125, 1.500, 1.900, 2.350, 2.875,
ccc     .               3.500, 4.250, 5.100, 6.025, 7.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=13)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.125, 0.250, 0.375, 0.525, 0.725, 1.000, 1.375,
ccc     .               1.875, 2.500, 3.250, 4.100, 5.025, 6.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=13)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.20, 0.40, 0.62, 0.88, 1.20, 1.60, 2.10, 2.70,
ccc     .               3.40, 4.20, 5.08, 6.02, 7.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=13)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.333333, 0.666667, 1.016667, 1.400000, 1.833333,
ccc     .               2.333333, 2.916667, 3.583333, 4.333333, 5.166667,
ccc     .               6.066667, 7.016667, 8.000000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      parameter (numfin=15)
      real coord(0:numfin),xnorm(0:numfin-1)
      data coord/0., 0.30, 0.60, 0.90, 1.20, 1.50, 1.82, 2.18, 2.60,
     .               3.10, 3.70, 4.40, 5.20, 6.08, 7.02, 8.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=16)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.28, 0.56, 0.84, 1.12, 1.40, 1.68, 1.96, 2.24,
ccc     .               2.52, 2.84, 3.24, 3.76, 4.40, 5.16, 6.04, 7.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=25)
ccc      real coord(0:numfin)
ccc      data coord/
ccc  . 0.000, 0.300, 0.600, 0.900, 1.200, 1.500, 1.800, 2.100, 2.400
ccc  .,2.700, 3.000, 3.301, 3.608, 3.928, 4.266, 4.630, 5.024, 5.456
ccc  .,5.931, 6.456, 7.037, 7.680, 8.392, 9.178,10.046,11.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(refinp)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(refinp)'
      end if
c
c$OMP PARALLEL DO PRIVATE(wgt,icrs,ifin,valu0,valu1,sumc,sumf)
c$OMP+ SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
c
c --- make room for 2*(numfin - numcrs) new grid rows
c
c --- rows north of equatorial strip:
      do 2 i=1,ieqcrs-numcrs
 2    fieldf(i,j)=fieldc(i,j)
c
c --- rows south of equatorial strip:
      do 3 i=ieqcrs+numcrs,ii-2*newrows
 3    fieldf(i+2*newrows,j)=fieldc(i,j)
c
c --- equatorial row:
      fieldf(ieqfin,j)=fieldc(ieqcrs,j)
c
c --- fill in remaining rows by interpolation
c
      do 4 i=1,numfin-1
      wgt=coord(i)-int(coord(i))
c
c --- south of equator:
      ifin=ieqfin+i
      icrs=ieqcrs+int(coord(i))
      valu0=fieldc(icrs  ,j)
      valu1=fieldc(icrs+1,j)
c --- make sure no land data are used in interpolation
      if (ipcrs(icrs,j).eq.0 .and. ipcrs(icrs+1,j).eq.1) valu0=valu1
      if (ipcrs(icrs,j).eq.1 .and. ipcrs(icrs+1,j).eq.0) valu1=valu0
      fieldf(ifin,j)=wgt*valu1+(1.-wgt)*valu0
cdiag if (j.eq.6) write (*,'(a,3i4,3f8.2,2i3)')
cdiag. 'i,ifin,icrs,wgt,valu1,valu0:',i,ifin,icrs,wgt,valu1,valu0
cdiag.  ,ipcrs(icrs+1,j),ipcrs(icrs,j)
c
c --- north of equator:
      ifin=ieqfin-i
      icrs=ieqcrs-int(coord(i))
      valu0=fieldc(icrs  ,j)
      valu1=fieldc(icrs-1,j)
c --- make sure no land data are used in interpolation
      if (ipcrs(icrs,j).eq.0 .and. ipcrs(icrs-1,j).eq.1) valu0=valu1
      if (ipcrs(icrs,j).eq.1 .and. ipcrs(icrs-1,j).eq.0) valu1=valu0
      fieldf(ifin,j)=wgt*valu1+(1.-wgt)*valu0
cdiag if (j.eq.6) write (*,'(a,3i4,3f8.2,2i3)')
cdiag. 'i,ifin,icrs,wgt,valu1,valu0:',i,ifin,icrs,wgt,valu1,valu0
cdiag.  ,ipcrs(icrs-1,j),ipcrs(icrs,j)
c
 4    continue
c
 1    continue
c$OMP END PARALLEL DO
      return
c
c
      entry refnap(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- use this entry if 'field' contains values proportional to cell size,
c --- i.e., values that need to be apportioned rather than interpolated
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(refnap)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(refnap)'
      end if
c
c$OMP PARALLEL DO PRIVATE(x)
      do 11 j=1,jj
c
c --- make room for 2*(numfin - numcrs) new grid rows
c
c --- rows north of equatorial strip:
      do 12 i=1,ieqcrs-numcrs
 12   fieldf(i,j)=fieldc(i,j)
c
c --- rows south of equatorial strip:
      do 13 i=ieqcrs+numcrs,ii-2*newrows
 13   fieldf(i+2*newrows,j)=fieldc(i,j)
c
c --- apportion remaining unrefined-grid rows among refined-grid rows
c
      do 15 i=ieqfin-numfin+1,ieqfin+numfin-1
 15   fieldf(i,j)=0.
c
      do 14 icrs=0,numcrs-1
      do 14 ifin=0,numfin-1
      if (icrs+ifin.eq.0) then
        x=.5*coord(1)
      else
        x=min(real(icrs)+.5,max(real(icrs)-.5,
     .     .5*(coord(ifin)+coord(ifin+1))))
     .  - min(real(icrs)+.5,max(real(icrs)-.5,
     .     .5*(coord(ifin)+coord(ifin-1))))
      end if
c
c --- south of equator:
      fieldf(ieqfin+ifin,j)=fieldf(ieqfin+ifin,j)
     .   +fieldc(ieqcrs+icrs,j)*ipcrs(ieqcrs+icrs,j)*x
c
c --- north of equator:
      fieldf(ieqfin-ifin,j)=fieldf(ieqfin-ifin,j)
     .   +fieldc(ieqcrs-icrs,j)*ipcrs(ieqcrs-icrs,j)*x
 14   continue
c
 11   continue
c$OMP END PARALLEL DO
      return
c
c
      entry refinu(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- use this entry to operate on -u- points
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(refinu)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(refinu)'
      end if
c
c$OMP PARALLEL DO PRIVATE(x,wgt,icrs,ifin)
      do 5 j=1,jj
c
c --- make room for 2*(numfin - numcrs) new grid rows
c
c --- rows north of equatorial strip:
      do 6 i=1,ieqcrs-numcrs
 6    fieldf(i,j)=fieldc(i,j)
c
c --- rows south of equatorial strip:
      do 7 i=ieqcrs+numcrs,ii-2*newrows
 7    fieldf(i+2*newrows,j)=fieldc(i,j)
c
c --- fill in remaining rows by interpolation
c
      do 8 i=1,numfin
      x=.5*(coord(i)+coord(i-1)+1.)
      wgt=x-int(x)
c
c --- south of equator:
      ifin=ieqfin+i
      icrs=ieqcrs+int(x)
      fieldf(ifin,j)=wgt*fieldc(icrs+1,j)+(1.-wgt)*fieldc(icrs,j)
cdiag if (j.eq.20) write (*,'(a,3i5,3f6.2)')
cdiag. 'i,ifin,icrs,wgt,fieldc(icrs+1),fieldc(icrs):',
cdiag.  i,ifin,icrs,wgt,fieldc(icrs+1,j),fieldc(icrs,j)
c
c --- north of equator:
      ifin=ieqfin-i+1
      icrs=ieqcrs-int(x)+1
      fieldf(ifin,j)=wgt*fieldc(icrs-1,j)+(1.-wgt)*fieldc(icrs,j)
cdiag if (j.eq.20) write (*,'(a,3i5,3f6.2)')
cdiag. 'i,ifin,icrs,wgt,fieldc(icrs-1),fieldc(icrs):',
cdiag.  i,ifin,icrs,wgt,fieldc(icrs-1,j),fieldc(icrs,j)
 8    continue
c
 5    continue
c$OMP END PARALLEL DO
      return
c
c
      entry unrefp(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- remove extra grid rows introduced by previous calls to 'refinp'
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(unrfin)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(a,2i4,a,i4,f6.2)') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(unrfin)'
      end if
c
c$OMP PARALLEL DO PRIVATE(wgt,icrs,ifin)
      do 21 j=1,jj
c
c --- copy extratropical rows into reduced-size array
c
c --- rows north of equatorial strip:
      do 22 i=1,ieqcrs-numcrs
 22   fieldc(i,j)=fieldf(i,j)
c
c --- rows south of equatorial strip:
      do 23 i=ieqcrs+numcrs,ii-2*newrows
 23   fieldc(i,j)=fieldf(i+2*newrows,j)
c
c --- equatorial row:
      fieldc(ieqcrs,j)=fieldf(ieqfin,j)
c
c --- fill in remaining rows by interpolation
c
      icrs=0
      do 24 i=2,numfin-1
c
      if (int(coord(i)).le.icrs) go to 24
      icrs=coord(i)
      wgt=(coord(i)-icrs)/(coord(i)-coord(i-1))
c
c --- south of equator:
      ifin=ieqfin+i
      fieldc(ieqcrs+icrs,j)=wgt*fieldf(ifin-1,j)+(1.-wgt)*fieldf(ifin,j)
c
c --- north of equator:
      ifin=ieqfin-i
      fieldc(ieqcrs-icrs,j)=wgt*fieldf(ifin+1,j)+(1.-wgt)*fieldf(ifin,j)
 24   continue
      if (icrs.ne.int(coord(numfin-1))) then
        write (lp,'(2(a,i3))') 'icrs=',icrs,'  not',int(coord(numfin-1))
        stop '(unrefp)'
      end if
c
 21   continue
c$OMP END PARALLEL DO
      return
c
      entry unrefu(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- use this entry to remove extra rows of -u- points
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(unrefu)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(unrefu)'
      end if
c
c$OMP PARALLEL DO PRIVATE(x,wgt,icrs,ifin)
      do 25 j=1,jj
c
c --- copy extratropical rows into reduced-size array
c
c --- rows north of equatorial strip:
      do 26 i=1,ieqcrs-numcrs
 26   fieldc(i,j)=fieldf(i,j)
c
c --- rows south of equatorial strip:
      do 27 i=ieqcrs+numcrs,ii-2*newrows
 27   fieldc(i,j)=fieldf(i+2*newrows,j)
c
c --- fill in remaining rows by interpolation
c
      icrs=0
      do 28 i=2,numfin
      x=.5*(coord(i)+coord(i-1)+1.)
      if (int(x).le.icrs) go to 28
      icrs=x
      wgt=(x-icrs)/(.5*(coord(i)-coord(i-2)))
c
c --- south of equator:
      ifin=ieqfin+i
      fieldc(ieqcrs+icrs  ,j)=
     .   wgt*fieldf(ifin-1,j)+(1.-wgt)*fieldf(ifin,j)
c
c --- north of equator:
      ifin=ieqfin-i+1
      fieldc(ieqcrs-icrs+1,j)=
     .   wgt*fieldf(ifin+1,j)+(1.-wgt)*fieldf(ifin,j)
 28   continue
      if (icrs.ne.int(coord(numfin))) then
        write (lp,'(2(a,i3))') 'icrs=',icrs,'  not',int(coord(numfin))
       stop '(unrefu)'
      end if
c
 25   continue
c$OMP END PARALLEL DO
      return
      end
