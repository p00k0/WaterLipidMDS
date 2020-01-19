	   program Code_Lipowsky
       implicit none
C **********************************************************************

c This Code is from Lipowsky's article 1998, J.Chem.Phys. 108, No. 17
c      - Gives the dynamics for four particles lipids (first choice)
C***********************************************************************
C      Local parameters:
       integer i,ntherm,iloop
       real*8 ekt,edum1,edum2
C      1) Potential parameters:
       real*8 rcut,sig,eps,drr,flip(200001),vlip(200001),flj(200001),
     $ vlj(200001),fscr(200001),vscr(200001),sigsc
       integer nint
C      2) Box / Subbox parameters:
       real*8 sboxxy,sboxx,sboxy,sboxz,dt,xn(16384),yn(16384),zn(16384),
     $ vx(16384),vy(16384),vz(16384),temp
       integer nwat,nlip,ntot
C      3) Neighbour/link list:
       real*8 rboxxy,rboxz
       integer nbx,nby,nbz,ncell,map(17303),ihead(17303),list(17303),
     $ lip(16384)
C      4) Force:
       real*8 fx(16384),fy(16384),fz(16384),rcut2,sboxx2,sboxy2,sboxz2,
     $ epot,ekx,eky,ekz,ek
C      5) Monte Carlo:
       real*8 eold,random(100000)
C      6) Pressure varibles:
       real*8 Jxx(50),Jyy(50),Jzz(50),Jxxs(50),Jyys(50),Jzzs(50),
     $ vcorxx,vcoryy,vcorzz,volume,cond(50),delc,pxx(50),pyy(50),pzz(50)
     $ ,gamma(50),gamAve(50)
       integer PressC,nsc
C      7) Energy variables:
       real*8 Elipitt(200),Emax
       integer Elevel(10000),nEmax,nEi
C      Input/Output system files:
       open(15,file='Pos.dat',status='unknown')
       open(14,file='Spe.dat',status='unknown')
       open(16,file='Dat.dat',status='unknown')
       open(23,file='StreesP.dat',status='unknown')
C *********************************************
C      1) Parameters:
C *********************************************
       read(16,*) nwat
       read(16,*) nlip
       read(16,*)sboxxy
       write(*,*) sboxxy
       read(16,*)sboxz
       write(*,*) sboxz
       read(16,*) temp
       temp=1.35
       ntot=nwat+5*nlip
       do i=1,ntot
	 read(15,*) xn(i),yn(i),zn(i)
         read(14,*) vx(i),vy(i),vz(i)	
       enddo
       dt=.0005d0
       sboxx=sboxxy
       sboxy=sboxxy

C *********************************************
C      2) Potential Matrices:
C *********************************************
       sig=1.0
       sigsc=1.05
       eps=1.0
       rcut=2.5
       nint=200000
       rcut2=rcut**2
       drr=rcut2/(nint*1.0)

       call F_Lip(nint,drr,flip,vlip,sig,eps)
       call F_SoftCore(drr,nint,fscr,vscr,rcut,sigsc)
       call F_LenardJones(drr,nint,rcut,vlj,flj,sig)
C *********************************************
C      3) Neighbour List:
C *********************************************
       nbx=sboxx/rcut
       nby=sboxy/rcut
       nbz=sboxz/rcut
       ncell=nbx*nby*nbz
       call neighbour(map,ncell,nbx,nby,nbz)
       call liplist(nlip,lip,nwat)


C *********************************************
C      4) Link List and Force:
C *********************************************
       rboxxy=sboxxy/(nbx*1.0)
       rboxz=sboxz/(nbz*1.0)
       call links(xn,yn,zn,ihead,list,ntot,ncell,nbx,nby,nbz,rboxxy,
     $ rboxz)
       sboxx2=sboxx/2.d0
       sboxy2=sboxy/2.d0
       sboxz2=sboxz/2.d0

C *********************************************
C       5) Thermalization/MetroPolis algorithm:
C *********************************************
        call seeds(0)
        edum1=0.0
        edum2=0.0
        call Energy(xn,yn,zn,vlip,vscr,vlj,lip,ntot,nlip,nwat,drr,
     $      rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2,epot,map,
     $      ihead,list,ncell,sig,eps)
        write(*,*)'We"re starting with:',epot
        do i=1,10000
         call randv(random,100000)
          call Monte_Carlo(eps,sig,xn,yn,zn,ntot,epot,vlip,vscr,vlj,
     $     lip,nlip,nwat,drr,rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,
     $     sboxy2,sboxz2,map,ihead,list,ncell,temp,rboxxy,rboxz,nbx,nby,
     $     random)
          call links(xn,yn,zn,ihead,list,ntot,ncell,nbx,nby,nbz,rboxxy,
     $         rboxz)
          call Energy(xn,yn,zn,vlip,vscr,vlj,lip,ntot,nlip,nwat,drr,
     $      rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2,epot,map,
     $      ihead,list,ncell,sig,eps)
          if(mod(i,10).eq.0) then
           edum2=epot
           write(*,*)'The inititial energy:',edum2,edum1-edum2,i
           edum1=edum2
         endif
       enddo

C       do iloop=1,100000
C         call reajusteZ(zn,nlip,ntot,nwat,sboxz)
C	write(*,*) 'boom 1'
C         call moveac(xn,yn,zn,vx,vy,vz,fx,fy,fz,dt,sboxx,sboxy
C     $               ,sboxz,ntot)
C	write(*,*) 'boom 2'
C        call links(xn,yn,zn,ihead,list,ntot,ncell,nbx,nby,nbz,rboxxy,
C     $ rboxz)
C         call force(xn,yn,zn,fx,fy,fz,flip,fscr,vlip,vscr,vlj,
C     $    flj,lip,ntot,nlip,nwat,drr,rcut,rcut2,sboxx,sboxy,sboxz,
C     $   sboxx2,sboxy2,sboxz2,epot,map,ihead,list,ncell,pxx,pyy,
C     $    pzz,sig,eps,delc,nsc)
C         call movebc(ekx,eky,ekz,vx,vy,vz,fx,fy,fz,ntot,dt,vcorxx,
C     $  vcoryy,vcorzz)
C        if(mod(iloop,50).eq.0) then
C          ekt=ekx+eky+ekz
C          write(*,*) iloop,ekt,epot
C        endif
C       call tempscale(vx,vy,vz,temp,ntot)
C       enddo
       rewind(14)
       rewind(15)
       do i=1,ntot
         write(15,*) xn(i),yn(i),zn(i)
         write(14,*) vx(i),vy(i),vz(i)
       enddo
       flush(14)
       flush(15)
C***********************************************
C       6) Obtain Stress profile:
C***********************************************
       call reajusteZ(zn,nlip,ntot,nwat,sboxz)
       call ELevel1(Elipitt,Emax,Elevel,nEi)
       call Pres_1(Jxx,Jyy,Jzz,Jxxs,Jyys,Jzzs,PressC,sboxx,sboxy,
     $ sboxz,volume,nsc,cond,delc,gamma,gamAve)
       call force(xn,yn,zn,fx,fy,fz,flip,fscr,vlip,vscr,vlj,
     $    flj,lip,ntot,nlip,nwat,drr,rcut,rcut2,sboxx,sboxy,sboxz,
     $   sboxx2,sboxy2,sboxz2,epot,map,ihead,list,ncell,pxx,pyy,
     $    pzz,sig,eps,delc,nsc)


       write(*,*)'thepotential from force:',epot,epot*(ntot*1.0)
      ntherm=100000
       do iloop=1,ntherm
         call reajusteZ(zn,nlip,ntot,nwat,sboxz)
         call moveac(xn,yn,zn,vx,vy,vz,fx,fy,fz,dt,sboxx,sboxy
     $               ,sboxz,ntot)

         call links(xn,yn,zn,ihead,list,ntot,ncell,nbx,nby,nbz,rboxxy,
     $ rboxz)
         call force(xn,yn,zn,fx,fy,fz,flip,fscr,vlip,vscr,vlj,
     $    flj,lip,ntot,nlip,nwat,drr,rcut,rcut2,sboxx,sboxy,sboxz,
     $   sboxx2,sboxy2,sboxz2,epot,map,ihead,list,ncell,pxx,pyy,
     $    pzz,sig,eps,delc,nsc)
         call movebc(ekx,eky,ekz,vx,vy,vz,fx,fy,fz,ntot,dt,vcorxx,
     $  vcoryy,vcorzz)
         call Pres_2(Jxx,Jyy,Jzz,pxx,vcorxx,pyy,vcoryy,pzz,vcorzz,
     $    cond,delc,nsc,gamma)
         if(mod(iloop,50).eq.0) then
          ekt=ekx+eky+ekz
          write(*,*) iloop,ekt,epot
          if(mod(iloop,5000).eq.0) then
            call Pres_3(Jxx,Jyy,Jzz,Jxxs,Jyys,Jzzs,PressC,nsc,gamma,
     $           gamAve)
            if(mod(iloop,100000).eq.0) then
              rewind(14)
              rewind(15)
              do i=1,ntot
                write(15,*) xn(i),yn(i),zn(i)
                write(14,*) vx(i),vy(i),vz(i)
              enddo
	      call flush(14)
	      call flush(15)
            endif
          endif
        endif
       call tempscale(vx,vy,vz,temp,ntot)
       enddo
       call Pres_4(Jxxs,Jyys,Jzzs,PressC,volume,iloop,nsc,cond,
     $             gamAve)
C *********************************************
C       8) Write Output positions/velocities:
C *********************************************
       rewind(14)
       rewind(15)
       do i=1,ntot
        write(15,*) xn(i),yn(i),zn(i)
        write(14,*) vx(i),vy(i),vz(i)
       enddo
       call flush(14)
       call flush(15)
       stop
       end
       
       
       
C      Potential matrices:
       subroutine F_Lip(nint,drr,flip,vlip,sig,eps)
       implicit none
       real*8 flip(200001),vlip(200001),rsq,drr,ri,sig,eps,k2
       integer nint,i
       k2=5000.0
       do i=1,nint
        ri=(i-1)*drr+drr/2.d0
        rsq=dsqrt(ri)
        flip(i)=2.0*k2*(rsq-sig)/rsq
        vlip(i)=k2*(rsq-sig)**2
       enddo
       return
       end
       
       subroutine F_SoftCore(drr,nint,fscr,vscr,rcut,sigsc)
       implicit none
       real*8 ri,drr,rsq,fscr(200001),vscr(200001),rcut,sigsc,dum,dum1,
     $ dum2,dum3,dum4
       integer i,nint
       dum=sigsc**9
       do i=1,nint
        ri=(i-1)*drr+drr/2.d0
        rsq=dsqrt(ri)
        dum1=4.0*(dum*(1.0/rsq)*(1.0/ri)**4)
        dum2=4.0*(-9.0)*(dum*(1.0/ri)**5)*(1/rsq)
        dum3=4.0*(dum*(1.0/rcut)**9)
        dum4=4.0*(-9.0)*(dum*(1.0/rcut)**10)
        fscr(i)=dum2-dum4*(1/rsq)
        vscr(i)=dum1-dum3-(rsq-rcut)*dum4
       enddo
       return
       end

       subroutine F_LenardJones(drr,nint,rcut,vlj,flj,sig)
       implicit none
       real*8 sig,ri,rsq,drr,vlj(200001),flj(200001),rcut,rcut2,sig2,
     $ dum1,dum2,dum3,dum4
       integer nint,i
       rcut2=rcut*rcut
       sig2=sig**2
       do i=1,nint
         ri=(i-1)*drr+drr/2.d0
         rsq=dsqrt(ri)
         dum1=4.0*((sig2/ri)**6-(sig2/ri)**3)
         dum2=4.0*((-12.0)*(sig2**6)*(1/rsq)**13-(-6.0)*(sig2**3)*
     $    (1/rsq)**7)*(1/rsq)
         dum3=4.0*((sig2/rcut2)**6-(sig2/rcut2)**3)
         dum4=4.0*((-12.0)*(sig2**6)*(1/rcut)**13-(-6.0)*(sig2**3)*
     $    (1/rcut)**7)
         vlj(i)=dum1-dum3-(rsq-rcut)*dum4
         flj(i)=dum2-dum4*(1/rsq)
       enddo

       return
       end

       subroutine liplist(nlip,lip,nwat)
       implicit none
C      Global
       integer nlip,lip(16384),nwat
C      Local:
       integer count,idumax,ncount,ncount2,i
        count=1
        idumax=5*nlip
        ncount=0
        ncount2=0
       do 11 i=1,idumax
       if (mod(i,5).eq.1) lip(count)=nwat+1+ncount
       if (mod(i,5).eq.2) lip(count+1)=nwat+nlip+1+ncount2
       if (mod(i,5).eq.3) lip(count+2)=nwat+nlip+2+ncount2
       if (mod(i,5).eq.4) lip(count+3)=nwat+nlip+3+ncount2
       if (mod(i,5).eq.0.and.i.gt.1) then
         lip(count+4)=nwat+nlip+4+ncount2
         count=count+5
         ncount=ncount+1
         ncount2=ncount2+4
       endif
11     continue
       return
       end
       
       subroutine neighbour(map,ncell,nbx,nby,nbz)
       implicit none
       integer nbx,nby,nbz,ix,iy,iz,icell,imap,mapsiz,ixp,iyp,izp,ixm
       integer iym,izm
       integer map(17303),ncell
       mapsiz=13*ncell
       do ix=1,nbx
        do iy=1,nby
          do iz=1,nbz
          ixp=ix+1
          if(ixp.gt.nbx) ixp=1
          iyp=iy+1
          if(iyp.gt.nby) iyp=1
          iym=iy-1
          if(iym.lt.1) iym=nby
          izp=iz+1
          if(izp.gt.nbz) izp=1
          izm=iz-1
          if(izm.lt.1) izm=nbz
          ixm=ix-1
          if(ixm.lt.1) ixm=nbx
          imap=13*(icell(ix,iy,iz,nbx,nby)-1)
          map(imap+1)=icell(ixp,iy,iz,nbx,nby)
          map(imap+2)=icell(ixp,iyp,iz,nbx,nby)
          map(imap+3)=icell(ix,iyp,iz,nbx,nby)
          map(imap+4)=icell(ixm,iyp,iz,nbx,nby)
          map(imap+5)=icell(ixp,iy,izm,nbx,nby)
          map(imap+6)=icell(ixp,iyp,izm,nbx,nby)
          map(imap+7)=icell(ix,iyp,izm,nbx,nby)
          map(imap+8)=icell(ixm,iyp,izm,nbx,nby)
          map(imap+9)=icell(ixp,iy,izp,nbx,nby)
          map(imap+10)=icell(ixp,iyp,izp,nbx,nby)
          map(imap+11)=icell(ix,iyp,izp,nbx,nby)
          map(imap+12)=icell(ixm,iyp,izp,nbx,nby)
          map(imap+13)=icell(ix,iy,izp,nbx,nby)
        enddo
        enddo
        enddo
       return
      end

      integer function icell(ix,iy,iz,nbx,nby)
       icell=ix+(iy-1)*nbx+(iz-1)*nbx*nby
      return
      end


C      The frequently called Force/Mouvement Subroutines:

       subroutine links(xn,yn,zn,ihead,list,ntot,ncell,nbx,nby,nbz,
     $ rboxxy,rboxz)
       implicit none
       real*8 xn(16384),yn(16384),zn(16384),rboxxy,rboxz
       integer ihead(17303),list(17303),ncell,ntot,ix,iy,iz,jcell
       integer nbx,nby,nbz,i
       do 10 jcell=1,ncell
        ihead(jcell)=0
10     continue
       do 20 i=1,ntot
        iy=yn(i)/rboxxy
        iz=zn(i)/rboxz
        ix=xn(i)/rboxxy
        jcell=1+ix+iy*nbx+iz*nby*nbx
        list(i)=ihead(jcell)
        ihead(jcell)=i
20     continue
       return
       end

        subroutine force(xn,yn,zn,fx,fy,fz,flip,fscr,vlip,vscr,vlj,
     $    flj,lip,ntot,nlip,nwat,drr,rcut,rcut2,sboxx,sboxy,sboxz,
     $   sboxx2,sboxy2,sboxz2,epot,map,ihead,list,ncell,pxx,pyy,
     $    pzz,sig,eps,delc,nsc)

       implicit none
C Variable global:
       real*8 xn(16384),yn(16384),zn(16384),fx(16384),fy(16384)
     $,fz(16384),flip(200001),fscr(200001),vlip(200001),vscr(200001),
     $ vlj(200001),flj(200001),pxx(50),pyy(50),pzz(50)
       integer lip(16384),ntot,nlip,nwat
       real*8 sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2
       real*8 drr,rcut,eps,sig,rcut2,epot,delc
       integer map(17303),ihead(17303),list(17303),ncell,nsc
C Variable local:
       real*8 rxij,ryij,rzij,rsq,fxi,fyi,fzi,xi,yi,zi
       integer m,n,i,j,ifx,imap,nbh,loc,jc,k1,k2,nbbox,kind,kmn,kmx,idum
     $,kdum1,kdum2
C ----------------------------------------------------------------------
C On initialise les force:
C ----------------------------------------------------------------------

       epot=0.0
       do 1 i=1,ntot
        fx(i)=0.d0
        fy(i)=0.d0
        fz(i)=0.d0
1      continue
       do i=1,50
        pxx(i)=0.d0
        pyy(i)=0.d0
        pzz(i)=0.d0
       enddo
c***********************************************************************
c**********              FORCE DANS LES CHAINES             ************
c***********************************************************************
       do 11 i=1,nlip
         do 110 j=1,4
         m=lip((i-1)*5+j)
         n=lip((i-1)*5+j+1)
         rxij=xn(n)-xn(m)
         ryij=yn(n)-yn(m)
         rzij=zn(n)-zn(m)
C  Condition frontiere:
         if(rxij.gt.sboxx2) rxij=rxij-sboxx
         if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
         if(ryij.gt.sboxy2) ryij=ryij-sboxy
         if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
         if(rzij.gt.sboxz2) rzij=rzij-sboxz
         if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
         rsq=rxij*rxij+ryij*ryij+rzij*rzij
         if (rsq.gt.rcut2) then
            write(*,*) 'bad',m,n,rxij,ryij,rzij
            write(*,*)'Positions:',xn(n),yn(n),zn(i),xn(m),yn(m),zn(m)
            call flush(6)
            stop 'Prog4'
          endif
        ifx=rsq/drr
        epot=epot+vlip(ifx)
        fxi=rxij*flip(ifx)
        fyi=ryij*flip(ifx)
        fzi=rzij*flip(ifx)
C i) Find which box the particles are in:
        k1=zn(m)/delc+1
        k2=zn(n)/delc+1
        nbbox=abs(k1-k2)+1
C ii) Set boundary conditions:
        if(nbbox.ge.(nsc/2))then
          kdum1=k1
          kdum2=k2
          if(kdum1.gt.kdum2)then
            k2=k2+nsc
            nbbox=abs(k1-k2)+1
          endif
          if(kdum1.lt.kdum2)then
            k1=k1+nsc
            nbbox=abs(k1-k2)+1
          endif
        endif
C iii) Find which box is on top/bottom:
        if(k1.ge.k2)then
          kmx=k1
          kmn=k2
        endif
        if(k1.lt.k2)then
          kmx=k2
          kmn=k1
        endif
C iv) Set the values pxx,pyy,pzz:
        do kind=kmn,kmx
          idum=kind
          if(kind.gt.nsc)idum=idum-nsc
          pxx(idum)=pxx(idum)-rxij*fxi/(nbbox*1.0)
          pyy(idum)=pyy(idum)-ryij*fyi/(nbbox*1.0)
          pzz(idum)=pzz(idum)-rzij*fzi/(nbbox*1.0)
        enddo
        fx(m)=fx(m)+fxi
        fy(m)=fy(m)+fyi
        fz(m)=fz(m)+fzi
        fz(n)=fz(n)-fzi
        fx(n)=fx(n)-fxi
        fy(n)=fy(n)-fyi
110     continue
11      continue


c***********************************************************************
c**********          INTERACTIONS DEHORS DES CHAINES           *********
c***********************************************************************

       do 4 m=1,ncell
c Interactions dans la cellule (iptc avec jptc):

       i=ihead(m)
5      if (i.gt.0) then
C----------------------------------------------
C      1) Interactions avec eau: on veut que i=> eau
C----------------------------------------------
       if (i.le.nwat) then
         xi=xn(i)
         yi=yn(i)
         zi=zn(i)
         j=list(i)

C L'indice de la prochaine particule de cellulue:
6        if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then

c 2.1) Eau-eau (L-J):
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vlj(ifx)
               fxi=-eps*rxij*flj(ifx)
               fyi=-eps*ryij*flj(ifx)
               fzi=-eps*rzij*flj(ifx)
               k1=zn(i)/delc+1
               k2=zn(j)/delc+1
               nbbox=abs(k1-k2)+1
               if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
               endif
               if(k1.ge.k2)then
                 kmx=k1
                 kmn=k2
               endif
               if(k1.lt.k2)then
                 kmx=k2
                 kmn=k1
               endif
               do kind=kmn,kmx
                 idum=kind
                 if(kind.gt.nsc)idum=idum-nsc
                 pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                 pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                 pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
               enddo
               fx(i)=fx(i)+fxi
               fy(i)=fy(i)+fyi
               fz(i)=fz(i)+fzi
               fx(j)=fx(j)-fxi
               fy(j)=fy(j)-fyi
               fz(j)=fz(j)-fzi
             endif

c 2.2) Eau-tete (L-J):

             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
                 fxi=-eps*rxij*flj(ifx)
                 fyi=-eps*ryij*flj(ifx)
                 fzi=-eps*rzij*flj(ifx)
                 k1=zn(i)/delc+1
                 k2=zn(j)/delc+1
                 nbbox=abs(k1-k2)+1
                 if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 endif
                 if(k1.ge.k2)then
                   kmx=k1
                   kmn=k2
                 endif
                 if(k1.lt.k2)then
                   kmx=k2
                   kmn=k1
                 endif
                 do kind=kmn,kmx
                   idum=kind
                   if(kind.gt.nsc)idum=idum-nsc
                   pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                   pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                   pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                 enddo
                 fx(i)=fx(i)+fxi
                 fy(i)=fy(i)+fyi
                 fz(i)=fz(i)+fzi
                 fx(j)=fx(j)-fxi
                 fy(j)=fy(j)-fyi
                 fz(j)=fz(j)-fzi
               endif
             endif

c  2.3 Eau-queue (S-C repulsif):

             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
               fxi=-eps*rxij*fscr(ifx)
               fyi=-eps*ryij*fscr(ifx)
               fzi=-eps*rzij*fscr(ifx)
               k1=zn(i)/delc+1
               k2=zn(j)/delc+1
               nbbox=abs(k1-k2)+1
               if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
               endif
               if(k1.ge.k2)then
                 kmx=k1
                 kmn=k2
               endif
               if(k1.lt.k2)then
                 kmx=k2
                 kmn=k1
               endif
               do kind=kmn,kmx
                 idum=kind
                 if(kind.gt.nsc)idum=idum-nsc
                 pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                 pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                 pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
               enddo
               fx(i)=fx(i)+fxi
               fy(i)=fy(i)+fyi
               fz(i)=fz(i)+fzi
               fx(j)=fx(j)-fxi
               fy(j)=fy(j)-fyi
               fz(j)=fz(j)-fzi
             endif
c end r cutoff:
           endif
c On a le prochain indice,
           j=list(j)
           goto 6
c end j>0...
         endif
C ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************
         imap=13*(m-1)
         do 8 nbh=1,13
           loc=imap+nbh
           jc=map(loc)
           j=ihead(jc)
10         if (j.gt.0) then
             rxij=xi-xn(j)
             ryij=yi-yn(j)
             rzij=zi-zn(j)
             if(rxij.gt.sboxx2) rxij=rxij-sboxx
             if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
             if(ryij.gt.sboxy2) ryij=ryij-sboxy
             if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
             if(rzij.gt.sboxz2) rzij=rzij-sboxz
             if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
             rsq=rxij*rxij+ryij*ryij+rzij*rzij
             if(rsq.lt.rcut2) then
c 2.1) Eau-eau:
               if (j.le.nwat) then

                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
                 fxi=-eps*rxij*flj(ifx)
                 fyi=-eps*ryij*flj(ifx)
                 fzi=-eps*rzij*flj(ifx)
                 k1=zn(i)/delc+1
                 k2=zn(j)/delc+1
                 nbbox=abs(k1-k2)+1
                 if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 endif
                 if(k1.ge.k2)then
                   kmx=k1
                   kmn=k2
                 endif
                 if(k1.lt.k2)then
                   kmx=k2
                   kmn=k1
                 endif
                 do kind=kmn,kmx
                   idum=kind
                   if(kind.gt.nsc)idum=idum-nsc
                   pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                   pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                   pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                 enddo
                 fx(i)=fx(i)+fxi
                 fy(i)=fy(i)+fyi
                 fz(i)=fz(i)+fzi
                 fx(j)=fx(j)-fxi
                 fy(j)=fy(j)-fyi
                 fz(j)=fz(j)-fzi
               endif

c 2.2) Eau-tete (L-J):

               if (j.gt.nwat) then
                 if (j.le.(nwat+nlip)) then
                   ifx=rsq/drr
                   epot=epot+eps*vlj(ifx)
                   fxi=-eps*rxij*flj(ifx)
                   fyi=-eps*ryij*flj(ifx)
                   fzi=-eps*rzij*flj(ifx)
                   k1=zn(i)/delc+1
                   k2=zn(j)/delc+1
                   nbbox=abs(k1-k2)+1
                   if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                   endif
                   if(k1.ge.k2)then
                     kmx=k1
                     kmn=k2
                   endif
                   if(k1.lt.k2)then
                     kmx=k2
                     kmn=k1
                   endif
                   do kind=kmn,kmx
                     idum=kind
                     if(kind.gt.nsc)idum=idum-nsc
                     pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                     pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                     pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                   enddo
                   fx(i)=fx(i)+fxi
                   fy(i)=fy(i)+fyi
                   fz(i)=fz(i)+fzi
                   fx(j)=fx(j)-fxi
                   fy(j)=fy(j)-fyi
                   fz(j)=fz(j)-fzi
                 endif
               endif

c  2.3 Eau-queue (S-C repulsif):
               if (j.gt.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
                 fxi=-eps*rxij*fscr(ifx)
                 fyi=-eps*ryij*fscr(ifx)
                 fzi=-eps*rzij*fscr(ifx)
                 k1=zn(i)/delc+1
                 k2=zn(j)/delc+1
                 nbbox=abs(k1-k2)+1
                 if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                  endif
                  if(k1.ge.k2)then
                    kmx=k1
                    kmn=k2
                  endif
                  if(k1.lt.k2)then
                    kmx=k2
                    kmn=k1
                  endif
                  do kind=kmn,kmx
                    idum=kind
                    if(kind.gt.nsc)idum=idum-nsc
                    pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                    pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                    pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                 enddo
                 fx(i)=fx(i)+fxi
                 fy(i)=fy(i)+fyi
                 fz(i)=fz(i)+fzi
                 fx(j)=fx(j)-fxi
                 fy(j)=fy(j)-fyi
                 fz(j)=fz(j)-fzi
               endif
c r cutoff end if:
             endif
14           j=list(j)
             goto 10
c end j>0 (2ieme fois)
           endif
8        continue
c  fermeture de la condition sur i:
       endif

c----------------------------------------------
c 3) Interactions avec tete: on veut que i=> tete
C----------------------------------------------
       if (i.gt.nwat) then
         if (i.le.(nwat+nlip)) then
         xi=xn(i)
         yi=yn(i)
         zi=zn(i)
         j=list(i)

26       if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then
c 2.2) tete-eau (L-J):
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vlj(ifx)
               fxi=-eps*rxij*flj(ifx)
               fyi=-eps*ryij*flj(ifx)
               fzi=-eps*rzij*flj(ifx)
               k1=zn(i)/delc+1
               k2=zn(j)/delc+1
               nbbox=abs(k1-k2)+1
               if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                endif
                if(k1.ge.k2)then
                  kmx=k1
                  kmn=k2
                endif
                if(k1.lt.k2)then
                  kmx=k2
                  kmn=k1
                endif
                do kind=kmn,kmx
                  idum=kind
                  if(kind.gt.nsc)idum=idum-nsc
                  pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                  pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                  pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
               enddo
               fx(i)=fx(i)+fxi
               fy(i)=fy(i)+fyi
               fz(i)=fz(i)+fzi
               fx(j)=fx(j)-fxi
               fy(j)=fy(j)-fyi
               fz(j)=fz(j)-fzi
             endif
c 3.2) tete-tete (L-J):
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
                 fxi=-rxij*eps*flj(ifx)
                 fyi=-ryij*eps*flj(ifx)
                 fzi=-rzij*eps*flj(ifx)
                 k1=zn(i)/delc+1
                 k2=zn(j)/delc+1
                 nbbox=abs(k1-k2)+1
                 if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                  endif
                  if(k1.ge.k2)then
                    kmx=k1
                    kmn=k2
                  endif
                  if(k1.lt.k2)then
                    kmx=k2
                    kmn=k1
                  endif
                  do kind=kmn,kmx
                    idum=kind
                    if(kind.gt.nsc)idum=idum-nsc
                    pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                    pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                    pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                 enddo
                 fx(i)=fx(i)+fxi
                 fy(i)=fy(i)+fyi
                 fz(i)=fz(i)+fzi
                 fx(j)=fx(j)-fxi
                 fy(j)=fy(j)-fyi
                 fz(j)=fz(j)-fzi
               endif
             endif
c  3.3 tete-queue (S-C):
             if (j.gt.(nwat+nlip)) then
                ifx=rsq/drr
                epot=epot+eps*vscr(ifx)
                fxi=-eps*rxij*fscr(ifx)
                fyi=-eps*ryij*fscr(ifx)
                fzi=-eps*rzij*fscr(ifx)
                k1=zn(i)/delc+1
                k2=zn(j)/delc+1
                nbbox=abs(k1-k2)+1
                if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                endif
                if(k1.ge.k2)then
                  kmx=k1
                  kmn=k2
                endif
                if(k1.lt.k2)then
                  kmx=k2
                  kmn=k1
                endif
                do kind=kmn,kmx
                  idum=kind
                  if(kind.gt.nsc)idum=idum-nsc
                  pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                  pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                  pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                enddo
                fx(i)=fx(i)+fxi
                fy(i)=fy(i)+fyi
                fz(i)=fz(i)+fzi
                fx(j)=fx(j)-fxi
                fy(j)=fy(j)-fyi
                fz(j)=fz(j)-fzi
              endif
c end r cutoff:
            endif
c On a le prochain indice,
17          j=list(j)
            goto 26
c end j>0...
          endif

C ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************

          imap=13*(m-1)
          do 28 nbh=1,13
            loc=imap+nbh
            jc=map(loc)
            j=ihead(jc)
210         if (j.ne.0) then
              rxij=xi-xn(j)
              ryij=yi-yn(j)
              rzij=zi-zn(j)
              if(rxij.gt.sboxx2) rxij=rxij-sboxx
              if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
              if(ryij.gt.sboxy2) ryij=ryij-sboxy
              if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
              if(rzij.gt.sboxz2) rzij=rzij-sboxz
              if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
              rsq=rxij*rxij+ryij*ryij+rzij*rzij
              if(rsq.lt.rcut2) then
c 2.2) tete-eau (L-J):
                if (j.le.nwat) then
                  ifx=rsq/drr
                  epot=epot+eps*vlj(ifx)
                  fxi=-eps*rxij*flj(ifx)
                  fyi=-eps*ryij*flj(ifx)
                  fzi=-eps*rzij*flj(ifx)
                  k1=zn(i)/delc+1
                  k2=zn(j)/delc+1
                  nbbox=abs(k1-k2)+1
                  if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                  endif
                  if(k1.ge.k2)then
                    kmx=k1
                    kmn=k2
                  endif
                  if(k1.lt.k2)then
                    kmx=k2
                    kmn=k1
                  endif
                  do kind=kmn,kmx
                    idum=kind
                    if(kind.gt.nsc)idum=idum-nsc
                    pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                    pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                    pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                  enddo
                  fx(i)=fx(i)+fxi
                  fy(i)=fy(i)+fyi
                  fz(i)=fz(i)+fzi
                  fx(j)=fx(j)-fxi
                  fy(j)=fy(j)-fyi
                  fz(j)=fz(j)-fzi
               endif
c 3.2) tete-tete (L-J):
               if (j.le.(nwat+nlip)) then
                 if (j.gt.nwat) then
                   ifx=rsq/drr
                   epot=epot+eps*vlj(ifx)
                   fxi=-rxij*eps*flj(ifx)
                   fyi=-ryij*eps*flj(ifx)
                   fzi=-rzij*eps*flj(ifx)
                   k1=zn(i)/delc+1
                   k2=zn(j)/delc+1
                   nbbox=abs(k1-k2)+1
                   if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                   endif
                   if(k1.ge.k2)then
                     kmx=k1
                     kmn=k2
                   endif
                   if(k1.lt.k2)then
                     kmx=k2
                     kmn=k1
                   endif
                   do kind=kmn,kmx
                     idum=kind
                     if(kind.gt.nsc)idum=idum-nsc
                     pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                     pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                     pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                   enddo
                   fx(i)=fx(i)+fxi
                   fy(i)=fy(i)+fyi
                   fz(i)=fz(i)+fzi
                   fx(j)=fx(j)-fxi
                   fy(j)=fy(j)-fyi
                   fz(j)=fz(j)-fzi
                 endif
               endif
c  3.3 tete-queue (S-C):
               if (j.gt.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
                 fxi=-eps*rxij*fscr(ifx)
                 fyi=-eps*ryij*fscr(ifx)
                 fzi=-eps*rzij*fscr(ifx)
                 k1=zn(i)/delc+1
                 k2=zn(j)/delc+1
                 nbbox=abs(k1-k2)+1
                 if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 endif
                 if(k1.ge.k2)then
                   kmx=k1
                   kmn=k2
                 endif
                 if(k1.lt.k2)then
                   kmx=k2
                   kmn=k1
                 endif
                 do kind=kmn,kmx
                   idum=kind
                   if(kind.gt.nsc)idum=idum-nsc
                   pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                   pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                   pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                 enddo
                 fx(i)=fx(i)+fxi
                 fy(i)=fy(i)+fyi
                 fz(i)=fz(i)+fzi
                 fx(j)=fx(j)-fxi
                 fy(j)=fy(j)-fyi
                 fz(j)=fz(j)-fzi
               endif
c r cutoff end if:
             endif
214          j=list(j)
             goto 210
c end j>0 (2ieme fois)
           endif
c do sur les neighboors
28       continue
c  Les end if qui ferme la condition sur i:
        endif
       endif



c----------------------------------------------
c 3) Interactions avec queue: i=> queue
C----------------------------------------------
       if (i.gt.(nwat+nlip)) then
         xi=xn(i)
         yi=yn(i)
         zi=zn(i)
         j=list(i)
36       if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then
c  2.3 Queue-eau:
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
               fxi=-eps*rxij*fscr(ifx)
               fyi=-eps*ryij*fscr(ifx)
               fzi=-eps*rzij*fscr(ifx)
               k1=zn(i)/delc+1
               k2=zn(j)/delc+1
               nbbox=abs(k1-k2)+1
               if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
               endif
               if(k1.ge.k2)then
                 kmx=k1
                 kmn=k2
               endif
               if(k1.lt.k2)then
                 kmx=k2
                 kmn=k1
               endif
               do kind=kmn,kmx
                 idum=kind
                 if(kind.gt.nsc)idum=idum-nsc
                 pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                 pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                 pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
               enddo
               fx(i)=fx(i)+fxi
               fy(i)=fy(i)+fyi
               fz(i)=fz(i)+fzi
               fx(j)=fx(j)-fxi
               fy(j)=fy(j)-fyi
               fz(j)=fz(j)-fzi
             endif
c  3.3 queue-tete
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
                 fxi=-eps*rxij*fscr(ifx)
                 fyi=-eps*ryij*fscr(ifx)
                 fzi=-eps*rzij*fscr(ifx)
                 k1=zn(i)/delc+1
                 k2=zn(j)/delc+1
                 nbbox=abs(k1-k2)+1
                 if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 endif
                 if(k1.ge.k2)then
                   kmx=k1
                   kmn=k2
                 endif
                 if(k1.lt.k2)then
                   kmx=k2
                   kmn=k1
                 endif
                 do kind=kmn,kmx
                   idum=kind
                   if(kind.gt.nsc)idum=idum-nsc
                   pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                   pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                   pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                 enddo
                 fx(i)=fx(i)+fxi
                 fy(i)=fy(i)+fyi
                 fz(i)=fz(i)+fzi
                 fx(j)=fx(j)-fxi
                 fy(j)=fy(j)-fyi
                 fz(j)=fz(j)-fzi
               endif
             endif
c  4.3 queue-queue (L-J attractif):
             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+vlj(ifx)
               fxi=-eps*rxij*flj(ifx)
               fyi=-eps*ryij*flj(ifx)
               fzi=-eps*rzij*flj(ifx)
               k1=zn(i)/delc+1
               k2=zn(j)/delc+1
               nbbox=abs(k1-k2)+1
               if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
               endif
               if(k1.ge.k2)then
                 kmx=k1
                 kmn=k2
               endif
               if(k1.lt.k2)then
                 kmx=k2
                 kmn=k1
               endif
               do kind=kmn,kmx
                 idum=kind
                 if(kind.gt.nsc)idum=idum-nsc
                 pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                 pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                 pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
               enddo
               fx(i)=fx(i)+fxi
               fy(i)=fy(i)+fyi
               fz(i)=fz(i)+fzi
               fx(j)=fx(j)-fxi
               fy(j)=fy(j)-fyi
               fz(j)=fz(j)-fzi
             endif
c end r cutoff:
           endif
c On a le prochain indice,
37      j=list(j)
      goto 36
c end j>0...
      endif



c ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************
       imap=13*(m-1)
       do 38 nbh=1,13
         loc=imap+nbh
         jc=map(loc)
         j=ihead(jc)

310      if (j.ne.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if (rsq.lt.rcut2) then
c  2.3 Queue-eau:
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
               fxi=-eps*rxij*fscr(ifx)
               fyi=-eps*ryij*fscr(ifx)
               fzi=-eps*rzij*fscr(ifx)
               k1=zn(i)/delc+1
               k2=zn(j)/delc+1
               nbbox=abs(k1-k2)+1
               if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
               endif
               if(k1.ge.k2)then
                 kmx=k1
                 kmn=k2
               endif
               if(k1.lt.k2)then
                 kmx=k2
                 kmn=k1
               endif
               do kind=kmn,kmx
                 idum=kind
                 if(kind.gt.nsc)idum=idum-nsc
                 pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                 pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                 pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
               enddo
               fx(i)=fx(i)+fxi
               fy(i)=fy(i)+fyi
               fz(i)=fz(i)+fzi
               fx(j)=fx(j)-fxi
               fy(j)=fy(j)-fyi
               fz(j)=fz(j)-fzi
             endif
c  3.3 queue-tete
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
                 fxi=-eps*rxij*fscr(ifx)
                 fyi=-eps*ryij*fscr(ifx)
                 fzi=-eps*rzij*fscr(ifx)
                 k1=zn(i)/delc+1
                 k2=zn(j)/delc+1
                 nbbox=abs(k1-k2)+1
                 if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 endif
                 if(k1.ge.k2)then
                   kmx=k1
                   kmn=k2
                 endif
                 if(k1.lt.k2)then
                   kmx=k2
                   kmn=k1
                 endif
                 do kind=kmn,kmx
                   idum=kind
                   if(kind.gt.nsc)idum=idum-nsc
                   pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                   pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                   pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
                 enddo
                 fx(i)=fx(i)+fxi
                 fy(i)=fy(i)+fyi
                 fz(i)=fz(i)+fzi
                 fx(j)=fx(j)-fxi
                 fy(j)=fy(j)-fyi
                 fz(j)=fz(j)-fzi
               endif
             endif
c  4.3 queue-queue (L-J attractif):
             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+vlj(ifx)
               fxi=-eps*rxij*flj(ifx)
               fyi=-eps*ryij*flj(ifx)
               fzi=-eps*rzij*flj(ifx)
               k1=zn(i)/delc+1
               k2=zn(j)/delc+1
               nbbox=abs(k1-k2)+1
               if(nbbox.ge.(nsc/2))then
                 kdum1=k1
                 kdum2=k2
                 if(kdum1.gt.kdum2)then
                   k2=k2+nsc
                   nbbox=abs(k1-k2)+1
                 endif
                 if(kdum1.lt.kdum2)then
                   k1=k1+nsc
                   nbbox=abs(k1-k2)+1
                 endif
               endif
               if(k1.ge.k2)then
                 kmx=k1
                 kmn=k2
               endif
               if(k1.lt.k2)then
                 kmx=k2
                 kmn=k1
               endif
               do kind=kmn,kmx
                 idum=kind
                 if(kind.gt.nsc)idum=idum-nsc
                 pxx(idum)=pxx(idum)+rxij*fxi/(nbbox*1.0)
                 pyy(idum)=pyy(idum)+ryij*fyi/(nbbox*1.0)
                 pzz(idum)=pzz(idum)+rzij*fzi/(nbbox*1.0)
               enddo
               fx(i)=fx(i)+fxi
               fy(i)=fy(i)+fyi
               fz(i)=fz(i)+fzi
               fx(j)=fx(j)-fxi
               fy(j)=fy(j)-fyi
               fz(j)=fz(j)-fzi
             endif
c r cutoff end if:
           endif
47         j=list(j)
           goto 310
c end j>0 (2ieme fois)
         endif
38       continue
c  end condition i=queue:
       endif
       i=list(i)
       go to 5
c end loop i>0
       endif
4      continue
       epot=epot/(ntot*1.0)
       return
       end


       subroutine moveac(xn,yn,zn,vx,vy,vz,fx,fy,fz,dt,sboxx,sboxy
     $,sboxz,ntot)
      implicit none
C Variable Global:
      real*8 xn(16384),yn(16384),zn(16384),fx(16384),fy(16384)
     $,fz(16384),vx(16384),vy(16384),vz(16384),dt,sboxx,sboxy,sboxz
       integer ntot,i
C Variable Local
       real*8 dtsq2
       dtsq2=dt*dt/2.d0
       do 100 i=1,ntot
        xn(i)=xn(i)+dt*vx(i)+dtsq2*fx(i)
        yn(i)=yn(i)+dt*vy(i)+dtsq2*fy(i)
        zn(i)=zn(i)+dt*vz(i)+dtsq2*fz(i)
        if(yn(i).gt.sboxy) yn(i)=yn(i)-sboxy
        if(yn(i).le.0.) yn(i)=yn(i)+sboxy
        if(xn(i).gt.sboxx) xn(i)=xn(i)-sboxx
        if(xn(i).le.0) xn(i)=xn(i)+sboxx
        if(zn(i).gt.sboxz) zn(i)=zn(i)-sboxz
        if(zn(i).le.0.) zn(i)=zn(i)+sboxz
        vx(i)=vx(i)+.5d0*dt*fx(i)
        vy(i)=vy(i)+.5d0*dt*fy(i)
        vz(i)=vz(i)+.5d0*dt*fz(i)

100    continue
       return
       end




       subroutine movebc(ekx,eky,ekz,vx,vy,vz,fx,fy,fz,ntot,dt,vcorxx,
     $vcoryy,vcorzz)
       implicit none
       real*8 xn(16384),yn(16384),zn(16384),fx(16384),fy(16384)
     $,fz(16384),vx(16384),vy(16384),vz(16384),dtekx,eky,ekz,vcorxx
     $,vcoryy,vcorzz,ekx,dt,vcmx,vcmy,vcmz
       integer ntot,i
       ekx=0.d0
       eky=0.d0
       ekz=0.d0
       vcmx=0.
       vcmy=0.
       vcmz=0.
       vcorxx=0.d0
       vcoryy=0.d0
       vcorzz=0.d0
       do 200 i=1,ntot
         vx(i)=vx(i)+.5d0*dt*fx(i)
         vy(i)=vy(i)+.5d0*dt*fy(i)
         vz(i)=vz(i)+.5d0*dt*fz(i)
        vcmx=vcmx+vx(i)
        vcmy=vcmy+vy(i)
        vcmz=vcmz+vz(i)
        ekx=ekx+vx(i)**2
        eky=eky+vy(i)**2
        ekz=ekz+vz(i)**2
200    continue
       vcmx=vcmx/ntot
       vcmy=vcmy/ntot
       vcmz=vcmz/ntot
       do i=1,ntot
        vcorxx=vcorxx+(vx(i)-vcmx)**2
        vcoryy=vcoryy+(vy(i)-vcmy)**2
        vcorzz=vcorzz+(vz(i)-vcmz)**2
       enddo
       ekx=ekx/(2.d0*ntot)
       eky=eky/(2.d0*ntot)
       ekz=ekz/(2.d0*ntot)
       return
       end

       subroutine tempscale(vx,vy,vz,temp,ntot)
       implicit none
       real*8 TA,Temp,vx(16384),vy(16384),vz(16384),dum,superdum
       integer i,ntot
       TA=0.0
       do 11 i=1,ntot
         TA=TA+(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
11     continue
       superdum=1.0*ntot
       TA=(1.0/(3.0*superdum))*TA
       if (TA.lt.1E-3) then
         dum=1.0
       else
         dum=sqrt(temp/TA)
       endif
       do 13 i=1,ntot
        vx(i)=dum*vx(i)
        vy(i)=dum*vy(i)
        vz(i)=dum*vz(i)
13     continue

       return
       end

       subroutine ELevel1(Elipitt,Emax,Elevel,nEi)
       implicit none
       real*8 Elipitt(200),Emax
       integer Elevel(10000),nEi,i
       do i=1,200
          Elipitt(i)=0.0
       enddo
       do i=1,10000
         Elevel(i)=0
       enddo
       Emax=500.0
       nEi=0
       return
       end


       subroutine reajusteZ(zn,nlip,ntot,nwat,sboxz)
       implicit none
       real*8 zn(16384),sboxz,sboxxy,memCM
       integer i,ntot,nwat,nlip
       memCM=0.0
       do i=nwat+1,ntot
          memCM=memCM+zn(i)
       enddo
       memCM=memCM/(nlip*5.0)
C  Add to every particle the displacement of the membrane with respect
C  to the center of the cell...
       do i=1,ntot
         zn(i)=zn(i)+((sboxz/2.0)-memCM)
         if(zn(i).gt.sboxz) zn(i)=zn(i)-sboxz
         if(zn(i).le.0.) zn(i)=zn(i)+sboxz
       enddo
       return
       end




       
       subroutine Pres_1(Jxx,Jyy,Jzz,Jxxs,Jyys,Jzzs,PressC,sboxx,sboxy,
     $ sboxz,volume,nsc,cond,delz,gamma,gamAve)
       implicit none
       real*8 Jxx(50),Jyy(50),Jzz(50),Jxxs(50),Jyys(50),Jzzs(50),sboxx,
     $  sboxy,sboxz,volume,delz,cond(50),gamma(50),gamAve(50)
       integer PressC,nsc,i

       do i=1,50
         cond(i)=0.0
       enddo
       nsc=45
       cond(1)=0.0
       do i=1,nsc
         cond(i+1)=(1.0*i)*(sboxz/(nsc*1.0))
       enddo
       write(*,*)cond(1),cond(2)
       volume=sboxx*sboxy*(cond(2)-cond(1))
       delz=sboxz/(1.0*nsc)
       do i=1,50
         Jxx(i)=0.0
         Jyy(i)=0.0
         Jzz(i)=0.0
         Jxxs(i)=0.0
         Jyys(i)=0.0
         Jzzs(i)=0.0
         gamma(i)=0.0
         gamAve(i)=0.0
       enddo
       PressC=0
       return
       end


       subroutine Pres_2(Jxx,Jyy,Jzz,pxx,vcorxx,pyy,vcoryy,pzz,vcorzz,
     $    cond,delc,nsc,gamma)
       implicit none
       real*8 Jxx(50),Jyy(50),Jzz(50),pxx(50),vcorxx,pyy(50),vcoryy,
     $  pzz(50),vcorzz,cond(50),delc,gamma(50)
       integer i,nsc
       do i=1,nsc
         Jxx(i)=pxx(i)+Jxx(i)
         Jyy(i)=pyy(i)+Jyy(i)
         Jzz(i)=pzz(i)+Jzz(i)
         gamma(i)=gamma(i)+0.5*(pxx(i)+pyy(i))-pzz(i)
       enddo
       return
       end


       subroutine Pres_3(Jxx,Jyy,Jzz,Jxxs,Jyys,Jzzs,PressC,nsc,gamma,
     $  gamAve)
       implicit none
       real*8 Jxx(50),Jyy(50),Jzz(50),Jxxs(50),Jyys(50),Jzzs(50),
     $  gamma(50),gamAve(50)
       integer PressC,i,nsc
       PressC=PressC+1
       do i=1,nsc
         Jxxs(i)=Jxxs(i)+Jxx(i)/5000.0
         Jyys(i)=Jyys(i)+Jyy(i)/5000.0
         Jzzs(i)=Jzzs(i)+Jzz(i)/5000.0
         gamAve(i)=gamAve(i)+gamma(i)/5000.0
         Jxx(i)=0.0
         Jyy(i)=0.0
         Jzz(i)=0.0
         gamma(i)=0.0
       enddo
       return
       end

       subroutine Pres_4(Jxxs,Jyys,Jzzs,PressC,volume,iloop,nsc,cond,
     $  gamAve)
       implicit none
       real*8 Jxxs(50),Jyys(50),Jzzs(50),Ptot,volume,Pz,Pxy,dum1,dum2,
     $        dum3,cond(50),gamAve(50),dum4
       integer PressC,iloop,nsc,i
       write(*,*)nsc,PressC
       do i=1,nsc
         dum1=-Jxxs(i)/(1.0*PressC)
         dum2=-Jyys(i)/(1.0*PressC)
         dum3=-Jzzs(i)/(1.0*PressC)
         dum4=-gamAve(i)/(2.0*PressC*volume)
         Ptot=(1.0/(2.0*volume))*(0.5*(dum1+dum2)-dum3)
         Pz=(1.0/volume)*dum3
         Pxy=(1.0/volume)*0.5*(dum2+dum1)
         write(*,*)(cond(i)+cond(i+1))/2.0,Ptot,Pz,Pxy,dum4
         write(23,*)(cond(i)+cond(i+1))/2.0,Ptot,Pz,Pxy,dum4
       enddo
       call flush(23)
       return
       end

       
       
       
       
       
       
       subroutine Monte_Carlo(eps,sig,xn,yn,zn,ntot,eold,vlip,vscr,vlj,
     $  lip,nlip,nwat,drr,rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,
     $  sboxy2,sboxz2,map,ihead,list,ncell,temp,rboxxy,rboxz,nbx,nby,
     $  random)
       implicit none
       real*8 xn(16384),yn(16384),zn(16384),xnew,ynew,znew,eold,enew,
     $ random(100000),xdum(16384),ydum(16384),zdum(16384),vlip(200001),
     $ vscr(200001),vlj(200001),drr,rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,
     $ sboxy2,sboxz2,eps,sig,temp,echeck1,echeck2,rboxxy,rboxz
       integer i,j,isel,map(17303),ihead(17303),list(17303),ncell,n,nwat
     $ ,nlip,ntot,lip(16384),nbx,nby
       n=1
       do j=1,ntot
C      Step 1: Choose a particle
         isel=random(n)*ntot+1
         n=n+1
C      Step two: Vary its position by some value dr (x+dx,y+dy,z+dz)
C            - Here, we take dx=(0.005,-0.005),dy=(...), dz=(...)

         xnew=xn(isel)+0.01*(1.0-2.0*random(n))
         n=n+1
         ynew=yn(isel)+0.01*(1.0-2.0*random(n))
         n=n+1
         znew=zn(isel)+0.01*(1.0-2.0*random(n))
         n=n+1
         do i=1,ntot
           xdum(i)=xn(i)
           ydum(i)=yn(i)
           zdum(i)=zn(i)
         enddo

C      Step three: Calculate the energy of the old/new configuration
         call Energy2(xn,yn,zn,vlip,vscr,vlj,lip,ntot,nlip,nwat,
     $    drr,rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2,echeck1,
     $    map,ihead,list,ncell,sig,eps,isel,rboxxy,rboxz,nbx,nby)
         echeck2=eold-echeck1
         xdum(isel)=xnew
         ydum(isel)=ynew
         zdum(isel)=znew
         if(ydum(isel).gt.sboxy) ydum(isel)=ydum(isel)-sboxy
         if(ydum(isel).le.0.) ydum(isel)=ydum(isel)+sboxy
         if(xdum(isel).gt.sboxx) xdum(isel)=xdum(isel)-sboxx
         if(xdum(isel).le.0) xdum(isel)=xdum(isel)+sboxx
         if(zdum(isel).gt.sboxz) zdum(isel)=zdum(isel)-sboxz
         if(zdum(isel).le.0.) zdum(isel)=zdum(isel)+sboxz
         call Energy2(xdum,ydum,zdum,vlip,vscr,vlj,lip,ntot,nlip,nwat,
     $    drr,rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2,echeck1,
     $    map,ihead,list,ncell,sig,eps,isel,rboxxy,rboxz,nbx,nby)
         enew=echeck2+echeck1
C      Step four: Do the Monte-Carlo test(t=1 ==> 324K)
           if(exp(-(enew-eold)).gt.random(n)) then
C      Step five: Update the new configuration
             do i=1,ntot
               xn(i)=xdum(i)
               yn(i)=ydum(i)
               zn(i)=zdum(i)
             enddo
             eold=enew
             n=n+1
           endif
       enddo

       return
       end






        subroutine Energy(xn,yn,zn,vlip,vscr,vlj,lip,ntot,nlip,nwat,drr,
     $  rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2,epot,map,
     $  ihead,list,ncell,sig,eps)

       implicit none
C Variable global:
       real*8 xn(16384),yn(16384),zn(16384),vlip(200001),vscr(200001)
     $  ,vlj(200001)
       integer lip(16384),ntot,nlip,nwat
       real*8 sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2
       real*8 drr,rcut,eps,sig,rcut2,epot
       integer map(17303),ihead(17303),list(17303),ncell
C Variable local:
       real*8 rxij,ryij,rzij,rsq,fxi,fyi,fzi,xi,yi,zi
       integer m,n,i,j,ifx,imap,nbh,loc,jc
C ----------------------------------------------------------------------
C On initialise les force:
C ----------------------------------------------------------------------
       epot=0.0
c***********************************************************************
c**********              ENERGIE DANS LES CHAINES           ************
c***********************************************************************
       do 11 i=1,nlip
         do 110 j=1,4
         m=lip((i-1)*5+j)
         n=lip((i-1)*5+j+1)
         rxij=xn(n)-xn(m)
         ryij=yn(n)-yn(m)
         rzij=zn(n)-zn(m)
C  Condition frontiere:
         if(rxij.gt.sboxx2) rxij=rxij-sboxx
         if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
         if(ryij.gt.sboxy2) ryij=ryij-sboxy
         if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
         if(rzij.gt.sboxz2) rzij=rzij-sboxz
         if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
         rsq=rxij*rxij+ryij*ryij+rzij*rzij
         ifx=rsq/drr
         epot=epot+vlip(ifx)
110     continue
11      continue


c***********************************************************************
c**********          INTERACTIONS DEHORS DES CHAINES           *********
c***********************************************************************

       do 4 m=1,ncell
c Interactions dans la cellule (iptc avec jptc):

       i=ihead(m)
5      if (i.gt.0) then
C----------------------------------------------
C      1) Interactions avec eau: on veut que i=> eau
C----------------------------------------------
       if (i.le.nwat) then
         xi=xn(i)
         yi=yn(i)
         zi=zn(i)
         j=list(i)

C L'indice de la prochaine particule de cellulue:
6        if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then

c 2.1) Eau-eau (L-J):
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vlj(ifx)
             endif

c 2.2) Eau-tete (L-J):

             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
               endif
             endif

c  2.3 Eau-queue (S-C repulsif):

             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
             endif
c end r cutoff:
           endif
c On a le prochain indice,
           j=list(j)
           goto 6
c end j>0...
         endif
C ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************
         imap=13*(m-1)
         do 8 nbh=1,13
           loc=imap+nbh
           jc=map(loc)
           j=ihead(jc)
10         if (j.gt.0) then
             rxij=xi-xn(j)
             ryij=yi-yn(j)
             rzij=zi-zn(j)
             if(rxij.gt.sboxx2) rxij=rxij-sboxx
             if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
             if(ryij.gt.sboxy2) ryij=ryij-sboxy
             if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
             if(rzij.gt.sboxz2) rzij=rzij-sboxz
             if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
             rsq=rxij*rxij+ryij*ryij+rzij*rzij
             if(rsq.lt.rcut2) then
c 2.1) Eau-eau:
               if (j.le.nwat) then
                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
               endif

c 2.2) Eau-tete (L-J):

               if (j.gt.nwat) then
                 if (j.le.(nwat+nlip)) then
                   ifx=rsq/drr
                   epot=epot+eps*vlj(ifx)
                 endif
               endif

c  2.3 Eau-queue (S-C repulsif):
               if (j.gt.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
c r cutoff end if:
             endif
14           j=list(j)
             goto 10
c end j>0 (2ieme fois)
           endif
8        continue
c  fermeture de la condition sur i:
       endif

c----------------------------------------------
c 3) Interactions avec tete: on veut que i=> tete
C----------------------------------------------
       if (i.gt.nwat) then
         if (i.le.(nwat+nlip)) then
         xi=xn(i)
         yi=yn(i)
         zi=zn(i)
         j=list(i)

26       if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then
c 2.2) tete-eau (L-J):
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vlj(ifx)
             endif
c 3.2) tete-tete (L-J):
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
               endif
             endif
c  3.3 tete-queue (S-C):
             if (j.gt.(nwat+nlip)) then
                ifx=rsq/drr
                epot=epot+eps*vscr(ifx)
              endif
c end r cutoff:
            endif
c On a le prochain indice,
17          j=list(j)
            goto 26
c end j>0...
          endif

C ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************

          imap=13*(m-1)
          do 28 nbh=1,13
            loc=imap+nbh
            jc=map(loc)
            j=ihead(jc)
210         if (j.ne.0) then
              rxij=xi-xn(j)
              ryij=yi-yn(j)
              rzij=zi-zn(j)
              if(rxij.gt.sboxx2) rxij=rxij-sboxx
              if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
              if(ryij.gt.sboxy2) ryij=ryij-sboxy
              if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
              if(rzij.gt.sboxz2) rzij=rzij-sboxz
              if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
              rsq=rxij*rxij+ryij*ryij+rzij*rzij
              if(rsq.lt.rcut2) then
c 2.2) tete-eau (L-J):
                if (j.le.nwat) then
                  ifx=rsq/drr
                  epot=epot+eps*vlj(ifx)
               endif
c 3.2) tete-tete (L-J):
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                   ifx=rsq/drr
                   epot=epot+eps*vlj(ifx)
                 endif
               endif
c  3.3 tete-queue (S-C):
               if (j.gt.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
c r cutoff end if:
             endif
214          j=list(j)
             goto 210
c end j>0 (2ieme fois)
           endif
c do sur les neighboors
28       continue
c  Les end if qui ferme la condition sur i:
        endif
       endif



c----------------------------------------------
c 3) Interactions avec queue: i=> queue
C----------------------------------------------
       if (i.gt.(nwat+nlip)) then
         xi=xn(i)
         yi=yn(i)
         zi=zn(i)
         j=list(i)
36       if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then
c  2.3 Queue-eau:
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
             endif
c  3.3 queue-tete
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
             endif
c  4.3 queue-queue (L-J attractif):
             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+vlj(ifx)
             endif
c end r cutoff:
           endif
c On a le prochain indice,
37      j=list(j)
      goto 36
c end j>0...
      endif



c ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************
       imap=13*(m-1)
       do 38 nbh=1,13
         loc=imap+nbh
         jc=map(loc)
         j=ihead(jc)

310      if (j.ne.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if (rsq.lt.rcut2) then
c  2.3 Queue-eau:
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
             endif
c  3.3 queue-tete
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
             endif
c  4.3 queue-queue (L-J attractif):
             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+vlj(ifx)
             endif
c r cutoff end if:
           endif
47         j=list(j)
           goto 310
c end j>0 (2ieme fois)
         endif
38       continue
c  end condition i=queue:
       endif
       i=list(i)
       go to 5
c end loop i>0
       endif
4      continue
       return
       end


c       subroutine links2(xn,yn,zn,ihead,list,ntot,ncell,nbx,nby,nbz,
c     $ rboxxy,rboxz,Pnum)
c       implicit none
c       real*8 xn(16384),yn(16384),zn(16384),rboxxy,rboxz
c       integer ihead2(17303),list2(17303),ncell,ntot,ix,iy,iz,jcell
c       integer nbx,nby,nbz,i,Pnum
c       do 10 jcell=1,ncell
c        ihead2(jcell)=0
c10     continue
c        ihead2(jcell)=Pnum
c       do 20 i=1,ntot
c        iy=yn(i)/rboxxy
c        iz=zn(i)/rboxz
c        ix=xn(i)/rboxxy
c        jcell=1+ix+iy*nbx+iz*nby*nbx
c        list(i)=ihead2(jcell)
c
c20     continue
c       return
c       end




        subroutine Energy2(xn,yn,zn,vlip,vscr,vlj,lip,ntot,nlip,
     $  nwat,drr,rcut,rcut2,sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2,epot,
     $  map,ihead,list,ncell,sig,eps,Pnum,rboxxy,rboxz,nbx,nby)
       implicit none
C Variable global:
       real*8 xn(16384),yn(16384),zn(16384),vlip(200001),vscr(200001)
     $  ,vlj(200001),rboxxy,rboxz
       integer lip(16384),ntot,nlip,nwat
       real*8 sboxx,sboxy,sboxz,sboxx2,sboxy2,sboxz2
       real*8 drr,rcut,eps,sig,rcut2,epot
       integer map(17303),ihead(17303),list(17303),ncell,nbx,nby
C Variable local:
       real*8 rxij,ryij,rzij,rsq,fxi,fyi,fzi,xi,yi,zi
       integer m,n,i,j,ifx,imap,nbh,loc,jc,Pnum,ix,iy,iz,FstPart
C ----------------------------------------------------------------------
C On initialise les force:
C ----------------------------------------------------------------------

       epot=0.0
c***********************************************************************
c**********              FORCE DANS LES CHAINES             ************
c***********************************************************************
       if(Pnum.gt.nwat) then
       do 11 i=1,nlip
         do 110 j=1,4
         m=lip((i-1)*5+j)
         n=lip((i-1)*5+j+1)
         rxij=xn(n)-xn(m)
         ryij=yn(n)-yn(m)
         rzij=zn(n)-zn(m)
C  Condition frontiere:
         if(rxij.gt.sboxx2) rxij=rxij-sboxx
         if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
         if(ryij.gt.sboxy2) ryij=ryij-sboxy
         if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
         if(rzij.gt.sboxz2) rzij=rzij-sboxz
         if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
         rsq=rxij*rxij+ryij*ryij+rzij*rzij
         ifx=rsq/drr
         epot=epot+vlip(ifx)
110     continue
11      continue
        endif

c***********************************************************************
c**********          INTERACTIONS DEHORS DES CHAINES           *********
c***********************************************************************

       iy=yn(Pnum)/rboxxy
       iz=zn(Pnum)/rboxz
       ix=xn(Pnum)/rboxxy
       m=1+ix+iy*nbx+iz*nby*nbx
C      Interactions dans la cellule (iptc avec jptc):
       j=ihead(m)
       xi=xn(Pnum)
       yi=yn(Pnum)
       zi=zn(Pnum)
C----------------------------------------------
C      1) Interactions avec eau: on veut que i=> eau
C----------------------------------------------
       if (Pnum.le.nwat) then

C L'indice de la prochaine particule de cellulue:
6         if(Pnum.eq.j) j=list(j)
          if(j.ne.0.and.j.eq.list(j)) then
             write(*,*)j,list(j),list(list(j)),list(list(list(j))),Pnum,
     $       xn(j),yn(j),zn(j),ncell,m
             stop
           endif
          if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then

c 2.1) Eau-eau (L-J):
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vlj(ifx)
             endif

c 2.2) Eau-tete (L-J):

             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
               endif
             endif

c  2.3 Eau-queue (S-C repulsif):

             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
             endif
c end r cutoff:
           endif
c On a le prochain indice,
           j=list(j)
           goto 6
c end j>0...
         endif
C ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************
         imap=13*(m-1)
         do 8 nbh=1,13
           loc=imap+nbh
           jc=map(loc)
           j=ihead(jc)
10         if (j.gt.0) then
             rxij=xi-xn(j)
             ryij=yi-yn(j)
             rzij=zi-zn(j)
             if(rxij.gt.sboxx2) rxij=rxij-sboxx
             if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
             if(ryij.gt.sboxy2) ryij=ryij-sboxy
             if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
             if(rzij.gt.sboxz2) rzij=rzij-sboxz
             if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
             rsq=rxij*rxij+ryij*ryij+rzij*rzij
             if(rsq.lt.rcut2) then
c 2.1) Eau-eau:
               if (j.le.nwat) then

                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
               endif

c 2.2) Eau-tete (L-J):

               if (j.gt.nwat) then
                 if (j.le.(nwat+nlip)) then
                   ifx=rsq/drr
                   epot=epot+eps*vlj(ifx)
                 endif
               endif

c  2.3 Eau-queue (S-C repulsif):
               if (j.gt.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
c r cutoff end if:
             endif
14           j=list(j)
             goto 10
c end j>0 (2ieme fois)
           endif
8        continue
c  fermeture de la condition sur i:
       endif

c----------------------------------------------
c 3) Interactions avec tete: on veut que i=> tete
C----------------------------------------------
       if (Pnum.gt.nwat) then
         if (Pnum.le.(nwat+nlip)) then
26         if(Pnum.eq.j) j=list(j)
           if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then
c 2.2) tete-eau (L-J):
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vlj(ifx)
             endif
c 3.2) tete-tete (L-J):
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vlj(ifx)
               endif
             endif
c  3.3 tete-queue (S-C):
             if (j.gt.(nwat+nlip)) then
                ifx=rsq/drr
                epot=epot+eps*vscr(ifx)
              endif
c end r cutoff:
            endif
c On a le prochain indice,
17          j=list(j)
            goto 26
c end j>0...
          endif
C ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************

          imap=13*(m-1)
          do 28 nbh=1,13
            loc=imap+nbh
            jc=map(loc)
            j=ihead(jc)
210         if (j.ne.0) then
              rxij=xi-xn(j)
              ryij=yi-yn(j)
              rzij=zi-zn(j)
              if(rxij.gt.sboxx2) rxij=rxij-sboxx
              if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
              if(ryij.gt.sboxy2) ryij=ryij-sboxy
              if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
              if(rzij.gt.sboxz2) rzij=rzij-sboxz
              if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
              rsq=rxij*rxij+ryij*ryij+rzij*rzij
              if(rsq.lt.rcut2) then
c 2.2) tete-eau (L-J):
                if (j.le.nwat) then
                  ifx=rsq/drr
                  epot=epot+eps*vlj(ifx)
               endif
c 3.2) tete-tete (L-J):
               if (j.le.(nwat+nlip)) then
                 if (j.gt.nwat) then
                   ifx=rsq/drr
                   epot=epot+eps*vlj(ifx)
                 endif
               endif
c  3.3 tete-queue (S-C):
               if (j.gt.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
c r cutoff end if:
             endif
214          j=list(j)
             goto 210
c end j>0 (2ieme fois)
           endif
c do sur les neighboors
28       continue
c  Les end if qui ferme la condition sur i:
        endif
       endif


c----------------------------------------------
c 3) Interactions avec queue: i=> queue
C----------------------------------------------
       if (Pnum.gt.(nwat+nlip)) then
36       if(Pnum.eq.j) j=list(j)
         if (j.gt.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if(rsq.lt.rcut2) then
c  2.3 Queue-eau:
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
             endif
c  3.3 queue-tete
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
             endif
c  4.3 queue-queue (L-J attractif):
             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+vlj(ifx)
             endif
c end r cutoff:
           endif
c On a le prochain indice,
37      j=list(j)
       goto 36
c end j>0...
        endif

c ********************************************************
C On VERIFIE LES CELLULES VOISINANTES:
C ********************************************************
       imap=13*(m-1)
       do 38 nbh=1,13
         loc=imap+nbh
         jc=map(loc)
         j=ihead(jc)
310      if (j.ne.0) then
           rxij=xi-xn(j)
           ryij=yi-yn(j)
           rzij=zi-zn(j)
           if(rxij.gt.sboxx2) rxij=rxij-sboxx
           if(rxij.le.(-sboxx2)) rxij=rxij+sboxx
           if(ryij.gt.sboxy2) ryij=ryij-sboxy
           if(ryij.le.(-sboxy2)) ryij=ryij+sboxy
           if(rzij.gt.sboxz2) rzij=rzij-sboxz
           if(rzij.le.(-sboxz2)) rzij=rzij+sboxz
           rsq=rxij*rxij+ryij*ryij+rzij*rzij
           if (rsq.lt.rcut2) then
c  2.3 Queue-eau:
             if (j.le.nwat) then
               ifx=rsq/drr
               epot=epot+eps*vscr(ifx)
             endif
c  3.3 queue-tete
             if (j.gt.nwat) then
               if (j.le.(nwat+nlip)) then
                 ifx=rsq/drr
                 epot=epot+eps*vscr(ifx)
               endif
             endif
c  4.3 queue-queue (L-J attractif):
             if (j.gt.(nwat+nlip)) then
               ifx=rsq/drr
               epot=epot+vlj(ifx)
             endif
c r cutoff end if:
           endif
47         j=list(j)
           goto 310
c end j>0 (2ieme fois)
         endif
38       continue


c  end condition i=queue:
       endif
       return
       end


c Random number generator....
c   --------------------------------------------

       subroutine seeds(jseed)

 ! This routine sets the seeds for the random number generator randv.

 ! If the argument jseed > 0, this argument is used to calculate the

 ! five integers needed for randv. Otherwise the time of the day is used.

       double precision dseed,rng69069

       parameter (dseed=2d0**31)

       integer time,i1,i2,j1,j2,j3

       integer jseed

       common /rng_integers/ i1,i2,j1,j2,j3

       save /rng_integers/

       if (jseed .gt. 0) then

         iseed=jseed

       else

         it=time()

 ! it is now the number of seconds since Jan. 1, 1970. Since the last digits

 ! are most random, the number is flipped around.

         lt=log10(real(it))

         it_backward=0

         do i=0,lt

           idigit=it-(it/10)*10

           it_backward=it_backward+idigit*10**(lt-i)

           it=(it-idigit)/10

         end do

         iseed=dble(it_backward)*dseed/10d0**(lt+1)

       end if

       i1=dseed*rng69069(iseed)

       if ((i1/2)*2 .eq. i1)  i1=i1+1

       i2=dseed*rng69069(iseed)

       if ((i2/2)*2 .eq. i2)  i2=i2+1

       j1=dseed*rng69069(iseed)

       j2=dseed*rng69069(iseed)

       j3=dseed*rng69069(iseed)

       return

       end

c  -----------------------------------------------

       function rng69069(iseed)

 ! The linear congruential random number generator i=i*69069+1

 ! This is a "bad" random number generator, use it only to generate seeds,

 ! not in a real calculation!

       double precision const,randt,rng69069

       integer ia,iexp

       parameter (ia=69069,const=1048576d0)

       integer ir(2)

       equivalence (randt,ir)

       data iexp/z'41300000'/

       ir(2)=iexp

       iseed=iseed*ia+1

       ir(1)=iseed

       rng69069=randt-const

       return

       end

c  ------------------------------------------------

       subroutine randv(r,n)

 !*****

 !*    subroutine randv implements the following algorithm:

 !*

 !*    x(i) = (x(i-1)*x(i-2)) mod (2**32)

 !*    y(i)=  (y(i-3)-y(i-1)) mod (2**30 - 35)

 !*    z(i) = (x(i)-y(i)) mod (2**32)

 !*    r(i) = z(i)/(2**32)+0.5

 !*    x are odd integers in the range (-2**31,2**31)

 !*    y are integers in the range [0,2**30-35)

 !*    z are integers in the range [-2**31,2**31)

 !*    r is real [0,1)

 !*

 !*    Literature: G. A. Marsaglia, in "Computer Science and Statistics:

 !*       The Interface", L. Billard (ed.), (Elsevier, Amsterdam, 1985), p.3

 !*

 !*    randv calculates an real array r of lenth n of random numbers in [0,1).

 !*    At the first call of randv 5 integers i1,i2,j1,j2,j3 must be provided.

 !*    The first two, i1 and i2, must be odd integers.

 !*

 !* !! Warning: the implementation below works only on machines that     !!

 !* !! conform to the IEEE standard for real and double precision        !!

 !* !! numbers. Furthermore, it assumes that the equivalence statement   !!

 !* !! of a double precision number with an integer field of length 2    !!

 !* !! packs the first element of the integer field into the first word  !!

 !* !! of the double precision number and the second element into the    !!

 !* !! second word. This is the case for almost all machines. But on,    !!

 !* !! e.g., Intel processors it is the other way round. In that case    !!

 !* !! line 61 must be changed to ir(2)=iexp and ir(2) has to be         !!

 !* !! replaced by ir(1) in all occurences between line 71 and line 125. !!

 !*

 !

 !     If you want real random numbers (instead of double precision), uncomment

 !     the following statement and comment out the statement right after that.

 !     furthermore replace each line that assigns a value to ir(2) by

 !     the line that directly follows that line.

 !

 !     real r(n)

       double precision r(n)

       double precision randt,const

       integer ir(2),i1,i2,j1,j2,j3

       equivalence (randt,ir)

       common /rng_integers/ i1,i2,j1,j2,j3

       save /rng_integers/

 ! this module holds the five integers i1,i2,j1,j2,j3 that must be

 ! initialized, e.g., using the routine seeds

       parameter (icut=-256)

       data mody/1073741789/,const/1048576d0/,iexp/z'41300000'/

 ! iexp loads the double precision representation of the exponent value

 ! 2**(e-bias)=2**20. for bias=1023, e=1043 (this is the IEEE standard).

 ! icut corresponds to z'ffffff00', i.e., iand(ir(2),icut) sets the last 8

 ! bits of ir(2) to zero so that the conversion of randt to real still leaves

 ! a number within [0,1).

       data mask/z'7fffffff'/

 ! ensure j1, j2 and j3 are in the range [0,mody-1]

       j1=mod(iand(j1,mask),mody)

       j2=mod(iand(j2,mask),mody)

       j3=mod(iand(j3,mask),mody)

 ! the loop is unrolled by a factor of 6(=[2 x par]*[3 y par])

 ! to speed things up. it reduces the number of executed stmts per

 ! random number from 11 to 6. (cf. the 2nd do loop which does the

 ! last few r-numbers in a slightly less efficient way)

       lo=mod(n,6)

       n6=n-lo

 ! pack the exponent and the mantissa for the floating point number

       ir(2)=iexp

       do i=1,n6,6

 !  begin 1st of 6 iterations

         i2=i1*i2

         j3=j3-j1

 ! for j3 and j1 in the range [0,mody-1] this is the same as doing the

 ! operation y mod mody

         if (j3 .lt. 0) j3=j3+mody

 ! this should be z = (x-y) mod 2**32.  machine overflow should

 ! look after the mod  2**32 stuff.

         ir(1)=i2-j3

 !       ir(1)=iand(i2-j3,icut)

 ! randt is now in the form : (2**20)*(1+z31/2**21+...+z0/2**44)

 ! where zn is the value of the nth bit in z

 ! subtracting const=2**20 leaves a number in the range [0,1);

 ! setting the last 8 bits of randt to zero ensures that randt-const

 ! is in the range [0,1) even after the conversion to real.

         r(i)=randt-const

 !  begin 2nd of 6 iterations

         i1=i2*i1

         j2=j2-j3

         if (j2 .lt. 0) j2=j2+mody

         ir(1)=i1-j2

 !       ir(1)=iand(i1-j2,icut)

         r(i+1)=randt-const

 !  begin 3rd of 6 iterations

         i2=i1*i2

         j1=j1-j2

         if (j1 .lt. 0) j1=j1+mody

         ir(1)=i2-j1

 !       ir(1)=iand(i2-j1,icut)

         r(i+2)=randt-const

 !  begin 4th of 6 iterations

         i1=i2*i1

         j3=j3-j1

         if (j3 .lt. 0) j3=j3+mody

         ir(1)=i1-j3

 !       ir(1)=iand(i1-j3,icut)

         r(i+3)=randt-const

 !  begin 5th of 6 iterations

         i2=i1*i2

         j2=j2-j3

         if (j2 .lt. 0) j2=j2+mody

         ir(1)=i2-j2

 !       ir(1)=iand(i2-j2,icut)

         r(i+4)=randt-const

 !  begin 6th of 6 iterations

         i1=i2*i1

         j1=j1-j2

         if (j1 .lt. 0) j1=j1+mody

         ir(1)=i1-j1

 !       ir(1)=iand(i1-j1,icut)

         r(i+5)=randt-const

       end do

       do i=n6+1,n

         ix=i1*i2

         iy=j3-j1
          if (iy .lt. 0) iy=iy+mody
           i2=i1
           i1=ix
           j3=j2
           j2=j1
           j1=iy
           ir(1)=ix-iy
 !       ir(1)=iand(ix-iy,icut)
          r(i)=randt-const
        end do
       return
       end
c ***********************************************
c ******** end of subroutines *******************