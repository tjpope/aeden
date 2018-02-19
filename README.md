# aeden
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
