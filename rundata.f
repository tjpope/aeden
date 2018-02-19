!---------------------------------------------------------------------!
      module mytypes
!---------------------------------------------------------------------!
      type base
       integer,allocatable,dimension(:)::i
       integer,allocatable,dimension(:,:)::c
       double precision,allocatable,dimension(:)::a,d
      end type base
      type maxbase
       integer,dimension(6)::i
       integer,dimension(6,3)::c
       double precision,dimension(6)::a,d
      end type maxbase
      type gauss
       integer,dimension(3)::cart
       double precision::expo,coef
       double precision,dimension(3)::cent
      end type gauss
      type nuke
      integer::la,lb
      double precision::ba,ga,bb,gb,bp,br
      end type nuke
      type nrecval
      integer,dimension(3)::a,b,c,d
      double precision::ga,gb,gc,gd,gp,gq,gw
      double precision,dimension(3)::RA,RB,RC,RD,RP,RQ,RW
      end type nrecval
!---------------------------------------------------------------------!
      end module mytypes
!---------------------------------------------------------------------!
      module rundata
!---------------------------------------------------------------------!
      use mytypes
!---------------------------------------------------------------------!
!     input                                                           !
!---------------------------------------------------------------------!
      integer::ng,n(3),nion,nb,ne,np,n1,n2
      double precision::xmin,xmax,dx(3)
      double precision,allocatable,dimension(:)::zion,slat
      double precision,allocatable,dimension(:,:)::ion
      character(100)::syslab
      logical::verbose
!---------------------------------------------------------------------!
!     gauss-jacobi                                                    !
!---------------------------------------------------------------------!
      integer::ngj
      double precision,allocatable,dimension(:)::xgj,wgj
!---------------------------------------------------------------------!
!     running                                                         !
!---------------------------------------------------------------------!
      type(base) sto
      integer::nx,no
      double precision,allocatable,dimension(:)::c
      double precision,allocatable,dimension(:,:)::smat,kmat,nmat,xmat
      double precision,allocatable,dimension(:,:,:,:)::qmat
!---------------------------------------------------------------------!
!     output                                                          !
!---------------------------------------------------------------------!
      double precision::mu,etot
      double precision,allocatable,dimension(:,:,:)::den,mden,sden
!---------------------------------------------------------------------!
!       simple constants                                              !
!---------------------------------------------------------------------!
       double precision::zero,pi,rt2,pi52
       parameter(zero=1e-25,pi=3.1415926535897932385,rt2=0.707107,
     .                            pi52=17.49341833)
!---------------------------------------------------------------------!
!       composite constants for F- and T-matrices                     !
!---------------------------------------------------------------------!
       double precision::hcon,pi2on6
       parameter(hcon=1.128379167,pi2on6=1.644934067)
!---------------------------------------------------------------------!
!       exchange corraltion parameters                                !
!---------------------------------------------------------------------!
       double precision::chachiyo_a,chachiyo_b,chachiyo_c,lda_x
       parameter(chachiyo_a=-0.01554534534,chachiyo_b=20.4562557,
     .                chachiyo_c=1.611991954,lda_x=-0.7385587664)
!---------------------------------------------------------------------!
!       conversion units                                              !
!---------------------------------------------------------------------!
       double precision::eh2ev,bor2A,aF2nN,ry2ev
       parameter(eh2ev=27.211385,bor2A=0.52917772109,aF2nN=82.387,
     .   ry2ev=13.605693)
!---------------------------------------------------------------------!
      contains
      subroutine allocations
      implicit none
      if(verbose)then
       allocate(den(n(1),n(2),n(3)))!total density matrix
       allocate(mden(n(1),n(2),n(3)))!mass density matrix
       allocate(sden(n(1),n(2),n(3)))!spin density matrix
      endif
      allocate(sto%a(nb),sto%d(nb),sto%i(nb),sto%c(nb,3))!basis set
      allocate(smat(nb,nb),xmat(n2,n2))!gaussian overlap matrix 
      allocate(kmat(nb,nb))!gaussian kinetic matrix 
      allocate(nmat(nb,nb))!gaussian nuclear matrix 
      allocate(qmat(nb,nb,nb,nb))!Q matrix - BE VERY CAREFUL!!!!! 
      allocate(c(n2))!coefficients for density
      allocate(xgj(ngj),wgj(ngj))
      end subroutine allocations
!---------------------------------------------------------------------!
      end module rundata
!---------------------------------------------------------------------!
