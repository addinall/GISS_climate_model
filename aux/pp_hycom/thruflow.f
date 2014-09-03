      real function thrufl(iaa,jaa,ibb,jbb,text)
c
c --- compute thruflow through section (iaa,jaa) - (ibb,jbb)
c --- (it is recommended to put both end points on land)
c
      use hycom_arrays, only: ubavav,vbavav,depths,scuy,scvx
      use hycom_dimen
      use const_proc
c
      implicit none
c
      integer :: ia,ib,ja,jb
      real flo(2)
      integer iaa,jaa,ibb,jbb
      character text*(*)
c
ccc      write (lp,'(a,6i5)')
ccc     .   'thrufl called with iaa,jaa,ibb,jbb =',
ccc     .                       iaa,jaa,ibb,jbb
c --- reverse points if iaa > ibb
      if (iaa.gt.ibb) then
        ia=ibb-iorign+1
        ib=iaa-iorign+1
        ja=jbb-jorign+1
        jb=jaa-jorign+1
      else
        ia=iaa-iorign+1
        ib=ibb-iorign+1
        ja=jaa-jorign+1
        jb=jbb-jorign+1
      end if
c
      flo(1)=0.
      flo(2)=0.                  !                                 x
c                                !                                  x
      if (jb.ge.ja) then         !  cross section orientation:       x
c                                !                                    x
      l=0                        !                                     x
      do 2 j=ja,jb+1,jb-ja+1
      l=l+1
      do 2 i=ia,ib
      if (iv(i,j).eq.1) then
        flo(l)=flo(l)
     .    +vbavav(i,j)*min(depths(i,j),depths(i,j-1))*scvx(i,j)
c     write(*,'(a,2i4,es10.3,2f6.0)')
c    .  'thr1 i,j,v,topo=',i,j,vbavav(i,j),depths(i,j),depths(i,j-1)
ccc        do 6 k=1,kk
ccc 6      flo(l)=flo(l)+vflx(k,i,j)
      end if
 2    continue
c
      l=0
      do 3 i=ib+1,ia,ia-ib-1
      l=l+1
      do 3 j=ja,jb
      if (iu(i,j).eq.1) then
        flo(l)=flo(l)
     .   -ubavav(i,j)*min(depths(i,j),depths(i-1,j))*scuy(i,j)
c     write(*,'(a,2i4,es10.3,2f6.0)')
c    .  'thr1 i,j,u,topo=',i,j,ubavav(i,j),depths(i,j),depths(i-1,j)
ccc        do 7 k=1,kk
ccc 7      flo(l)=flo(l)-uflx(i,j,k)
      end if
 3    continue                   !                                     x
c                                !                                    x
      else       !  (jb < ja)    !  cross section orientation:       x
c                                !                                  x
      l=0                        !                                 x
      do 4 j=jb,ja+1,ja-jb+1
      l=l+1
      do 4 i=ia,ib
      if (iv(i,j).eq.1) then
        flo(l)=flo(l)
     .    +vbavav(i,j)*min(depths(i,j),depths(i,j-1))*scvx(i,j)
c       write(*,'(a,2i4,es10.3,2f6.0)')
c    .  'thr2 i,j,v,topo=',i,j,vbavav(i,j),depths(i,j),depths(i,j-1)
ccc        do 8 k=1,kk
ccc 8      flo(l)=flo(l)+vflx(i,j,k)
      end if
 4    continue
c
      l=0
      do 5 i=ia,ib+1,ib-ia+1
      l=l+1
      do 5 j=jb,ja
      if (iu(i,j).eq.1) then
        flo(l)=flo(l)
     .   +ubavav(i,j)*min(depths(i,j),depths(i-1,j))*scuy(i,j)
c     write(*,'(a,2i4,es10.3,2f6.0)')
c    .  'thr2 i,j,u,topo=',i,j,ubavav(i,j),depths(i,j),depths(i-1,j)
ccc        do 9 k=1,kk
ccc 9      flo(l)=flo(l)+uflx(k,i,j)
      end if
 5    continue
c
      end if
c --- convert to sverdrups (ubavav,vbavav in cm/s):
      flo(1)=flo(1)*sign(1,ibb-iaa)*1.e-8
      flo(2)=flo(2)*sign(1,ibb-iaa)*1.e-8
c
c     write (lp,'(2f6.1,2(a,2i5),a,2x,a)') flo(1),flo(2),
c    .   '  transport between (',iaa,jaa,') and (',ibb,jbb,')',text
c
      thrufl=.5*(flo(1)+flo(2))
      return
      end
