!---------------------------------------------------------------------!
!   aeden                                                             !
!---------------------------------------------------------------------!
!                                                                     !
!   Program for calculating all-electron ground state density and     !
!    energy in small atoms/molecules                                  !
!   We use the Extended Electron Fock matrix:                         !
!                                                                     !
!         | H    -B |                                                 !
!      F= |         | ,                                               !
!         | B     H |                                                 !
!                                                                     !
!   where                                                             !
!      H = Kinetic + Nuclear + Coulomb                                !
!      B = Bivector                                                   !
!   Details on the Extended Electron model can be found at:           !
!                                                                     !
!   - Werner A Hofer. Unconventional approach to orbital-free density !
!     functional theory derived from a model of extended electrons.   !
!     Foundations of Physics, 41(4):754–791, 2011.                    !
!   - Thomas Pope and Werner Hofer. Spin in the extended electron     !
!     model. Frontiers of Physics, 12(3):128503, 2017.                !
!   - Thomas Pope and Werner Hofer. An Extended Electron Approach to  !
!     the General Many-Body Problem. arXiv preprint arXiv:1801.06242, !
!     2018.                                                           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!   The overlap, kinetic and nuclear integrals are generated using    !
!    the recursive algorithms presented repsectively in:              !
!                                                                     !
!   - Minhhuy Ho and Julio Manuel Hernandex-Peresz. Evaluation of     !
!     Gaussian Molecular Integrals I. Overlap Integrals. The          !
!     Mathematica Journal, 14, 2012.                                  !
!   - Minhhuy Ho and Julio Manuel Hernandex-Peresz. Evaluation of     !
!     Gaussian Molecular Integrals II. Kinetic Energy Integrals. The  !
!     Mathematica Journal, 15, 2013.                                  !
!   - Minhhuy Ho and Julio Manuel Hernandex-Peresz. Evaluation of     !
!     Gaussian Molecular Integrals III. Nuclear-Electron Attraction   !
!     Integrals. The Mathematica Journal, 16, 2014.                   !
!                                                                     !
!   The coulomb matrix is evaluated using the recursive algorithm     !
!    presneted in:                                                    !
!                                                                     !
!   - Shigeru Obara and A Saika. Efficient recursive computation of   !
!     molecular integrals over Cartesian Gaussian functions. The      !
!     Journal of chemical physics, 84(7):3963–3974, 1986.             !
!                                                                     !
!   Numerical integration for the nuclear matrix is performed using   !
!    Chebyshev-Gauss quadrature as implemented in:                    !
!                                                                     !
!   - William H Press, Saul A Teukolsky, William T Vetterling, and    !
!     Brian P Flannery. Numerical recipes in Fortran 77: the art of   !
!     scientific computing, 1992.                                     !
!                                                                     !
!   The basis sets available are the STO-nG type, where n=3,4,5 and   !
!    6. We use the values for orbital exponents and initial           !
!    coefficients presented in:                                       !
!                                                                     !
!   - Warren J Hehre, Robert F Stewart, and John A Pople. self-       !
!     consistent molecular-orbital methods. i. use of gaussian        !
!     expansions of Slater-type atomic orbitals. The Journal of       !
!     Chemical Physics, 51(6):2657–2664, 1969                         !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!   The LAPACK library is used - specifically, the routines:          !
!                                                                     !
!   - dsyev                                                           !
!   - dgetrf                                                          !
!   - dgetri                                                          !
!                                                                     !
!   Details found at:                                                 !
!                                                                     !
!   - Edward Anderson, Zhaojun Bai, Christian Bischof, Susan          !
!     Blackford, James Demmel, Jack Dongarra, Jeremy Du Croz, Anne    !
!     Greenbaum, S Hammerling, Alan McKenney, et al. LAPACK Users'    !
!     guide, volume 9. Siam, 1999.                                    !
!                                                                     !
!---------------------------------------------------------------------!
      program aeden
!---------------------------------------------------------------------!
      use rundata       !global variables
      use io            !input/output routines
      use basis         !basis set routines
      use scf           !scf routines
      implicit none
      call moses        !get input from command line
      call input        !get input from files
      call allocations  !allocate arrays
      call initbase     !initialize the chosen basis set
      if(animate)write(*,100)nani; iani=1
      call smatrix      !calculate overlap matrix
c      call output("smat.dat",n0,n0,1,smat)
      call kmatrix      !calculate kinetic matrix
c      call output("kmat.dat",n0,n0,1,kmat)
      call nmatrix      !calculate nuclear matrix
c      call output("nmat.dat",n0,n0,1,nmat)
      call qmatrix      !calculate coulomb matrix
      call kohnsham     !perform scf cycle
      call energies     !calculate and output energies
      if(animate) call interate !interate if requested
      call output("coef.dat",nh,1,1,c)
      if(verbose)then   !generate mass and spin densities on grid
       if(hartfck)then 
        call output("dens.dat",n(1),n(2),n(3),plotbase(c))
       else
        call output("mass.dat",n(1),n(2),n(3),plotbase(c(1:n0)))
        call output("spin.dat",n(1),n(2),n(3),plotbase(c(n1:n2)))
       endif
      endif
100   format("#",t63,"#",/,"#",t6,"Running on ",i0," Geometries",
     .       t63,"#",/,"#",t63,"#")
      end program aeden
!---------------------------------------------------------------------!
      subroutine interate
      use rundata       !global variables
      use io            !input/output routines
      use basis         !basis set routines
      use renew_basis   !basis renewal routines
      use scf           !scf routines
      implicit none
      double precision::steps,cutoff,cut
      cutoff=nani/63.0; steps=2; cut=cutoff
      do iani=2,nani
       if(uout.ne.6) then
        steps=steps+1; if(steps.ge.cut) then
         call system("echo -n '#'"); cut=cut+cutoff
        endif
       endif
c       call renew_smatrix !calculate overlap matrix again
c       call renew_kmatrix !calculate kinetic matrix again
c       call renew_nmatrix !calculate nuclear matrix again
       call smatrix !calculate overlap matrix again
c      call output("smat.dat",n0,n0,1,smat)
       call kmatrix !calculate kinetic matrix again
c      call output("kmat.dat",n0,n0,1,kmat)
       call nmatrix !calculate nuclear matrix again
c      call output("nmat.dat",n0,n0,1,nmat)
       call qmatrix       !calculate coulomb matrix again
       call kohnsham      !perform scf cycle
       call energies      !calculate and output energies
c       stop
      enddo
      if(uout.ne.6)then
       close(65); write(*,200); write(*,300)anitime
      else
       write(*,300)anitime
      endif
200   format(/,"#",t63,"#",/,"#",t6,
     .       "The rest of the ouput is in output.aeden",t63,"#",/,
     .       "#",t63,"#")
300   format(63("#"),/,"#",t63,"#",/,"#",t6,
     .       "Overall CPU Times:",t63,"#",/,"#"t63,"#",/,
     .  "#",t6,"Overlap Matrix: ",t26,f10.3,t36," seconds",t63,"#",/,
     .  "#",t6,"Kinetic Matrix: ",t26,f10.3,t36," seconds",t63,"#",/,
     .  "#",t6,"Nuclear Matrix: ",t26,f10.3,t36," seconds",t63,"#",/,
     .  "#",t6,"Coulomb Matrix: ",t26,f10.3,t36," seconds",t63,"#",/,
     .  "#",t6,"SCF: ",t26,f10.3,t36," seconds",t63,"#",/,
     .  "#"t63,"#",/,63("#"))
      end subroutine interate
!---------------------------------------------------------------------!
