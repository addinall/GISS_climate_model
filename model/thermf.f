      subroutine thermf(m,n,mm,nn,k1m,k1n)
c
      USE SEAICE, only : fsss
c
c --- hycom version 0.9
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : baclin,thref,watcum,empcum,nstep,nstep0
     &     ,diagno,lp,area,spcifh,avgbot,g,onem,slfcum
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, GLOBALSUM
      implicit none
c
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n
c
      real thknss,radfl,radflw,radfli,vpmx,prcp,prcpw,prcpi,
     .     evap,evapw,evapi,exchng,target,old,
     .     rmean,tmean,smean,vmean,boxvol,emnp(idm,J_0H:J_1H),
     &     slfcol(J_0H:J_1H),watcol(J_0H:J_1H),empcol(J_0H:J_1H),
     &     rhocol(J_0H:J_1H),temcol(J_0H:J_1H),salcol(J_0H:J_1H)
      integer iprime,ktop
      real qsatur,totl,eptt,salrlx
      external qsatur
      data ktop/3/
ccc      data ktop/2/                        !  normally set to 3
css   data salrlx/0.3215e-7/          !  1/(1 yr)
      data salrlx/0.6430e-8/          !  1/(5 yr)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional: weak salinity restoring via corrective surface flux
c
ccc      if (mod(time+.0001,30.).lt..0002) then            !  once a month
ccc        totl=0.
ccc        eptt=0.
cccc$OMP PARALLEL DO REDUCTION(+:totl,eptt) SCHEDULE(STATIC,jchunk)
ccc        do 8 j=1,jj
ccc        do 8 k=1,kk
ccc        do 8 l=1,isp(j)
ccc        do 8 i=ifp(j,l),ilp(j,l)
ccc        eptt=eptt+oemnp(i,j)*scp2(i,j)
ccc 8      totl=totl+saln(i,j,k)*dp(i,j,k)*scp2(i,j)
cccc$OMP END PARALLEL DO
ccc        totl=totl/g                                     !  10^-3 kg
cccc
cccc --- initialize total salt content
ccc        if (saltot.le.0.) saltot=totl
cccc --- determine corrective flux using multi-year relaxation time scale
ccc        pcpcor=(totl-saltot)*salrlx                     !  10^-3 kg/sec
ccc     .   *thref/35.                                     !  => m^3/sec
ccc     .   /eptt                                          !  => rel. units
ccc        write (lp,'(i9,a,f9.5,a,f7.2)') nstep,
ccc     .  '  overall salt gain (psu):',(totl-saltot)*thref/(area*avgbot),
ccc     .  '  => precip bias(%):',100.*pcpcor
ccc      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c$OMP PARALLEL DO PRIVATE(kn) SCHEDULE(STATIC,jchunk)
      do 66 j=J_0,J_1
      do 66 k=1,kk
      kn=k+nn
      do 66 l=1,isp(j)
      do 66 i=ifp(j,l),ilp(j,l)
      if (glue(i,j).gt.1. .and. saln(i,j,kn).gt.40.)                !  Med fudge
     .  saln(i,j,kn)=saln(i,j,kn)+(40.-saln(i,j,kn))*baclin*3.e-8
 66   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c$OMP END PARALLEL DO
c
c --- --------------------------------
c --- thermal forcing of ocean surface
c --- --------------------------------
c
      rmean=0.
      tmean=0.
      smean=0.
c
c$OMP PARALLEL DO PRIVATE(thknss,vpmx,prcp,exchng,
c$OMP. radfl,radflw,radfli,evap,evapw,evapi) SCHEDULE(STATIC,jchunk)
      do 85 j=J_0,J_1
c
      watcol(j)=0.
      empcol(j)=0.
      slfcol(j)=0.
      rhocol(j)=0.
      temcol(j)=0.
      salcol(j)=0.
c
      do 85 l=1,isp(j)
c
      do 85 i=ifp(j,l),ilp(j,l)
c     ustar(i,j)=sqrt(sqrt(taux(i,j)**2+tauy(i,j)**2)*thref)
c    .                             *(1.-oice(i,j))
c
      surflx(i,j)=oflxa2o(i,j)
c
c --- emnp = evaporation minus precipitation over open water (m/sec)
css   emnp(i,j)=(evapw*thref/evaplh + prcp)*(1.-covice(i,j))
css   emnp(i,j)=oemnp(i,j)*(1.+pcpcor)
      emnp(i,j)=oemnp(i,j)
c --- salflx = salt flux (10^-3 kg/m^2/sec) in +p direction, salt from SI added
c     salflx(i,j)=saln(i,j,k1n)*(osalt(i,j)*fsss-emnp(i,j)/thref)
css   salflx(i,j)=-35.0*emnp(i,j)/thref+osalt(i,j)*1.e3
      salflx(i,j)=-saln(i,j,k1n)*emnp(i,j)/thref+osalt(i,j)*1.e3
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,100) nstep,i,j,
cdiag.  '   taux      tauy      surflx      emnp       oice    ustar',
cdiag. taux(i,j),tauy(i,j),surflx(i,j),emnp(i,j),oice(i,j),ustar(i,j)
 100    format(i9,2i5,a/18x,1p,7e10.3/19x,a/18x,5e10.3)
c
      watcol(j)=watcol(j)+surflx(i,j)*scp2(i,j)
      empcol(j)=empcol(j)+emnp(i,j)*scp2(i,j)
      slfcol(j)=slfcol(j)+salflx(i,j)*scp2(i,j)
      rhocol(j)=rhocol(j)+th3d(i,j,k1n)*scp2(i,j)
      temcol(j)=temcol(j)+temp(i,j,k1n)*scp2(i,j)
      salcol(j)=salcol(j)+saln(i,j,k1n)*scp2(i,j)
 85   continue
c$OMP END PARALLEL DO
c
      call GLOBALSUM(ogrid,watcol,watcum, all=.true.)
      call GLOBALSUM(ogrid,empcol,empcum, all=.true.)
      call GLOBALSUM(ogrid,slfcol,slfcum, all=.true.)
c
      if (nstep.eq.nstep0+1 .or. diagno) then

      call GLOBALSUM(ogrid,rhocol,rmean, all=.true.)
      call GLOBALSUM(ogrid,temcol,tmean, all=.true.)
      call GLOBALSUM(ogrid,salcol,smean, all=.true.)
c
      if( AM_I_ROOT() ) then
      write (lp,'(i9,''energy residual (w/m^2)'',f9.2)') nstep,
     .  watcum/(area*(nstep-nstep0))
      write (lp,'(9x,''resulting temp drift (deg/century):'',f9.3)')
     .  watcum*thref/(spcifh*avgbot*area*(nstep-nstep0)) *36500.*86400.
      write (lp,'(i9,''e - p residual (mm/year)'',f9.2)') nstep,
     .  empcum/(area*(nstep-nstep0))*86400.*365000.
      write (lp,'(9x,''saln drift resulting from e-p (psu/century):''
     .                                                     ,f9.3)')
     .  -empcum*35./(avgbot*area*(nstep-nstep0)) *36500.*86400.
css   write (lp,'(i9,''salt residual (T/year)'',3f9.2)') nstep,
css  .  slfcum*365.*86400.*g/onem
      write (lp,'(7x,''saln drift resulting from salfl (psu/century):''
     .                                                     ,f9.3)')
     .  slfcum*36500.*86400.*g/(avgbot*area*onem)
      write (lp,'(i9,a,3f9.3)') nstep,' mean surf. sig,temp,saln:',
     .    rmean/area,tmean/area,smean/area
      end if ! AM_I_ROOT
c
      rmean=0.
      smean=0.
      tmean=0.
      vmean=0.
      ! reusing these arrays for next sums 
      rhocol=0.; temcol=0.; salcol=0.; watcol=0.;
c$OMP PARALLEL DO PRIVATE(kn,boxvol) SCHEDULE(STATIC,jchunk)
      do 84 j=J_0,J_1
      do 84 k=kk,1,-1
      kn=k+nn
      if (nstep.eq.nstep0+1) kn=k+mm
      do 84 l=1,isp(j)
      do 84 i=ifp(j,l),ilp(j,l)
      boxvol=dp(i,j,kn)*scp2(i,j)
      !--changing following to facilitate GLOBALSUM with bitwise reprocibility
      ! rmean=rmean+th3d(i,j,kn)*boxvol
      ! smean=smean+saln(i,j,kn)*boxvol
      ! tmean=tmean+temp(i,j,kn)*boxvol
      ! vmean=vmean+boxvol
      !-----------------------------------------------
      rhocol(j)=rhocol(j)+th3d(i,j,kn)*boxvol
      salcol(j)=salcol(j)+saln(i,j,kn)*boxvol
      temcol(j)=temcol(j)+temp(i,j,kn)*boxvol
      watcol(j)=watcol(j)+boxvol
 84   continue
c$OMP END PARALLEL DO
      call GLOBALSUM(ogrid,rhocol,rmean, all=.true.) 
      call GLOBALSUM(ogrid,salcol,smean, all=.true.) 
      call GLOBALSUM(ogrid,temcol,tmean, all=.true.) 
      call GLOBALSUM(ogrid,watcol,vmean, all=.true.) 
c
      if( AM_I_ROOT() )
     .    write (lp,'(i9,a,3f9.3)') nstep,' mean basin sig,temp,saln:',
     .    rmean/vmean,tmean/vmean,smean/vmean
      end if                                !  diagno = .true.
c
      return
      end
c
c
c> Revision history
c>
c> Mar. 2000 - conversion to SI units
c> Apr. 2000 - added diagnostics of t/s trends implied by srf.flux residuals
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
c> Nov. 2000 - added global density diagnostis to global t/s diagnostics
c> Apr. 2001 - eliminated stmt_funcs.h
c> Apr. 2001 - added diapycnal flux forcing at lateral boundaries
c> June 2001 - eliminated -buoyfl- calculation (now done in mxlayr)
c> Aug. 2001 - improved handling of flux calculations over water and ice
c> Dec. 2001 - fixed bug in -evapw- calculation (replaced max by min)
