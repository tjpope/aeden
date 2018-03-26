!---------------------------------------------------------------------!
      module scf
!---------------------------------------------------------------------!
      contains
      subroutine eigensolve(n,h,c)
      use rundata, only: nx,xmat,ne
      implicit none
      integer::i,n,info,np
      double precision,dimension(n)::c
      double precision,dimension(nx)::e
      double precision,dimension(3*nx)::wrk
      double precision,dimension(n,n)::h,v
      double precision,dimension(nx,nx)::h1
      intent(in)n,h
      h1=matmul(transpose(xmat(:,1:nx)),matmul(h,xmat(:,1:nx)))
      call dsyev('V','U',nx,h1,nx,e,wrk,3*nx,info)
      if(info/=0)stop'H-matrix diagonalization failed '
      v=matmul(xmat(:,1:nx),h1);
      c=v(:,1)*sqrt(dble(ne))
c      do i=1,n; c(i)=sum(v(i,1:ne)); enddo; c=c/sqrt(dble(ne))
      end subroutine eigensolve
!---------------------------------------------------------------------!      
      subroutine kohnsham
      use rundata, only: n0,n1,n2,nh,ne,smat,kmat,nmat,qmat,c,
     .                               mu,etot,hartfck,cputime,tol
      implicit none
      integer::i,j,k,l,cnt
      double precision::eold,enew,enuc,converge,mix
      double precision,dimension(n0,n0)::H0,S,RHO
      double precision,dimension(nh,nh)::H,P
      call cpu_time(cputime(5))
      converge=tol
      H0=kmat+nmat; S=smat; enuc=ionrep()
      P=fullden(c)
      if(hartfck) then
       RHO=P(1:n0,1:n0)
      else
       RHO=P(n1:n2,n1:n2)+P(1:n0,1:n0)
      endif
      H=focker(kmat,nmat,qmat,c,P,RHO)
c      open(62,file="conv.dat")
c      open(63,file="convcof.dat")
c      open(64,file="convrho.dat")
c      do i=1,nh; write(63,*) -1,i,c(i); enddo; write(63,*)
c      do i=1,n0; write(64,*) -1,i,RHO(i,i); enddo; write(64,*)       
c      write(62,*) -1,enew,enuc,mu
      call eigensolve(nh,H,c)
      P=fullden(c); 
      if(hartfck) then
       RHO=P(1:n0,1:n0)
      else
       RHO=P(n1:n2,n1:n2)+P(1:n0,1:n0)
      endif
      H=focker(kmat,nmat,qmat,c,P,RHO) 
      enew=sum(P*H); eold=2*enew
c      do i=1,nh; write(63,*) 0,i,c(i); enddo; write(63,*)
c      do i=1,n0; write(64,*) 0,i,RHO(i,i); enddo; write(64,*)       
c      write(62,*) 0,enew,enuc,mu
      cnt=0; do while(abs(enew-eold)>converge)
       cnt=cnt+1
c       do i=1,nh; write(63,*) cnt,i,c(i); enddo; write(63,*)
c       do i=1,n0; write(64,*) cnt,i,RHO(i,i); enddo; write(64,*)       
       call eigensolve(nh,H,c)
       P=fullden(c); 
       if(hartfck) then
        RHO=P(1:n0,1:n0)
       else
        RHO=P(n1:n2,n1:n2)+P(1:n0,1:n0)
       endif
       H=focker(kmat,nmat,qmat,c,P,RHO); 
       eold=enew; enew=sum(P*H)
       mu=enew/ne+enuc
c       write(62,*) cnt,enew,enuc,mu
c       if(cnt.gt.1e5)eold=enew
       if(cnt.gt.1e5) then; write(*,100); stop; endif
      enddo
c      close(62)
c      close(63)
c      close(64)
      etot=mu
100   format("#",t63,"#",/,"#",t6,"This isn't working. I'm giving up",
     .           t63,"#",/,"#",t63,"#",/,63("#"))
      end subroutine kohnsham
!---------------------------------------------------------------------!
      function focker(T,V,Q,c,P,RHO) result(F)
      use rundata, only: n0,n1,n2,nh,hartfck
      implicit none
      integer::i,j
      double precision,dimension(nh)::c
      double precision,dimension(n0,n0)::T,V,RHO
      double precision,dimension(nh,nh)::P,F
      double precision,dimension(n0,n0,n0,n0)::Q
      intent(in)T,V,Q,c,P,RHO
      F(1:n0,1:n0)=(T+V)+coulomb(RHO,Q)
      if(.not.hartfck) then
       F(n1:n2,n1:n2)=F(1:n0,1:n0)
       F(n1:n2,1:n0)=bivector(c,T,RHO)
       F(1:n0,n1:n2)=-transpose(F(n1:n2,1:n0))
      endif
      end function focker
!---------------------------------------------------------------------!
      function coulomb(RHO,Q) result(VH)
      use rundata, only: n0,hartfck,ne
      implicit none
      integer::i,j
      double precision,dimension(n0,n0)::RHO,VH
      double precision,dimension(n0,n0,n0,n0)::Q
      intent(in)RHO,Q
      if(hartfck)then 
       do i=1,n0; do j=1,n0
        VH(i,j)=sum(RHO*(2*Q(i,j,:,:)-Q(i,:,j,:)))
       enddo; enddo
      else
       do i=1,n0; do j=1,n0
        VH(i,j)=sum(RHO*Q(i,j,:,:))
       enddo; enddo
      endif
      VH=0.25*VH*(ne-1)/ne
      end function coulomb
!---------------------------------------------------------------------!
      function bivector(c,T,RHO) result(B) 
      use rundata, only: n0,n1,n2
      implicit none
      integer::i,j
      double precision,dimension(n2)::c
      double precision,dimension(n0,n0)::T,RHO,B
      intent(in)c,T,RHO
      do i=1,n0; B(i,:)=(c(i)*c(n1:n2)-c(1:n0)*c(n0+i))*T(i,:); enddo
      B=B/(RHO+1.e-10)
      end function bivector
!---------------------------------------------------------------------!
      function fullden(c) result(P)
      use rundata, only: n0,n1,n2,nh,ne,hartfck
      implicit none
      integer::i
      double precision,dimension(nh)::c
      double precision,dimension(nh,nh)::p
      intent(in)c
      do i=1,n0; P(i,1:n0)=c(i)*c(1:n0); enddo
      if(.not.hartfck)then
       P(1:n0,n1:n2)=0.d0; P(n1:n2,1:n0)=P(1:n0,n1:n2)
       do i=1,n0; P(i,1:n0)=c(i)*c(1:n0); enddo
       do i=n1,n2; P(i,n1:n2)=c(i)*c(n1:n2); enddo
      endif
      P=P*ne
      end function fullden
!---------------------------------------------------------------------!
      function ionrep() result(E)
      use rundata
      implicit none
      integer::i,j
      double precision::E,d
      E=0
      do i=1,nion; do j=1,nion
       d=sqrt(sum((ion(i,:,iani)-ion(j,:,iani))**2))
       if(d>1e-8)E=E+zion(i)*zion(j)/d
      enddo; enddo
      end function
!---------------------------------------------------------------------!
      end module scf
!---------------------------------------------------------------------!
