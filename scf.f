!---------------------------------------------------------------------!
      module scf
!---------------------------------------------------------------------!
      contains
      subroutine eigensolve(n,h,c)
      use rundata, only: ne,nx,xmat,rt2
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
      v=matmul(xmat(:,1:nx),h1); c=v(:,1)
      end subroutine eigensolve
!---------------------------------------------------------------------!      
      subroutine kohnsham
      use rundata, only: nb,smat,kmat,nmat,qmat,c,mu,rt2,np,etot,n1,n2
      implicit none
      integer::i,j,k,l,cnt
      double precision::eold,enew,enuc,converge,mix,lam
      double precision,dimension(nb,nb)::H0,S,RHO
      double precision,dimension(np,n2)::cold
      double precision,dimension(n2,n2)::H,P
      converge=1e-9
      H0=kmat+nmat; S=smat; enuc=ionrep()
      P=fullden(c,nb); RHO=P(n1:n2,n1:n2)+P(1:nb,1:nb)
      H=focker(kmat,nmat,qmat,c,P,RHO)
      call eigensolve(n2,H,c)
      P=fullden(c,nb); RHO=P(n1:n2,n1:n2)+P(1:nb,1:nb)
      H=focker(kmat,nmat,qmat,c,P,RHO); enew=sum(P*H); eold=2*enew
      open(62,file="conv.dat")
      open(63,file="convcof.dat")
      do i=1,nb*2; write(63,*) 0,i,c(i); enddo; write(63,*)
      write(62,*) 0,enew,enuc,mu
      cnt=0; do while(abs(enew-eold)>converge)
       cnt=cnt+1
       do i=1,nb*2; write(63,*) cnt,i,c(i); enddo; write(63,*)
       call eigensolve(n2,H,c)
       P=fullden(c,nb); RHO=P(n1:n2,n1:n2)+P(1:nb,1:nb)
       eold=enew; H=focker(kmat,nmat,qmat,c,P,RHO); enew=sum(P*H)
       mu=enew+enuc
       write(62,*) cnt,enew,enuc,mu,lam; flush(62)
c       if(cnt.gt.1e5)eold=enew
       if(cnt.gt.1e5)stop"WHAT YOU'RE DOING!!! Ran out of steps..."
      enddo
      close(62)
      close(63)
      etot=mu
      end subroutine kohnsham
!---------------------------------------------------------------------!
      function focker(T,V,Q,c,P,RHO) result(F)
      use rundata, only: nb,n1,n2
      implicit none
      integer::i,j
      double precision,dimension(n2)::c
      double precision,dimension(nb,nb)::T,V,RHO
      double precision,dimension(n2,n2)::P,F
      double precision,dimension(nb,nb,nb,nb)::Q
      intent(in)T,V,Q,c,P,RHO
      F(1:nb,1:nb)=2*(T+V)+coulomb(RHO,Q)
      F(n1:n2,n1:n2)=F(1:nb,1:nb)
      F(n1:n2,1:nb)=bivector(c,T,RHO)
      F(1:nb,n1:n2)=-transpose(F(n1:n2,1:nb))
      end function focker
!---------------------------------------------------------------------!
      function coulomb(RHO,Q) result(VH)
      use rundata, only: nb
      implicit none
      integer::i,j
      double precision,dimension(nb,nb)::RHO,VH
      double precision,dimension(nb,nb,nb,nb)::Q
      intent(in)RHO,Q
      do i=1,nb; do j=1,nb
       VH(i,j)=sum(RHO*Q(i,j,:,:))
      enddo; enddo
      end function coulomb
!---------------------------------------------------------------------!
      function bivector(c,T,RHO) result(B) 
      use rundata, only: nb,n1,n2
      implicit none
      integer::i,j
      double precision,dimension(n2)::c
      double precision,dimension(nb,nb)::T,RHO,B
      intent(in)c,T,RHO
      do i=1,nb; B(i,:)=(c(i)*c(n1:n2)-c(1:nb)*c(nb+i))*T(i,:); enddo
      B=B/RHO
      end function bivector
!---------------------------------------------------------------------!
      function fullden(c,n) result(P)
      implicit none
      integer::i,n,n1,n2
      double precision,dimension(n*2)::c
      double precision,dimension(n*2,n*2)::p
      intent(in)c,n
      n1=n+1; n2=n*2; P(1:n,n1:n2)=0.d0; P(n1:n2,1:n)=P(1:n,n1:n2)
      do i=1,n; P(i,1:n)=c(i)*c(1:n); enddo
      do i=n1,n2; P(i,n1:n2)=c(i)*c(n1:n2); enddo
      end function fullden
!---------------------------------------------------------------------!
      function ionrep() result(E)
      use rundata
      implicit none
      integer::i,j
      double precision::E,d
      E=0
      do i=1,nion; do j=1,nion
       d=sqrt(sum((ion(i,:)-ion(j,:))**2))
       if(d>1e-8)E=E+zion(i)*zion(j)/d
      enddo; enddo
      end function
!---------------------------------------------------------------------!
      end module scf
!---------------------------------------------------------------------!
