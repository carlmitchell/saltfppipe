SUBROUTINE voigtfit(num,wave,flux,sigma,guess,flag,fit,err,chisq)
IMPLICIT NONE

INTEGER, INTENT(IN) :: num
REAL, DIMENSION(num), INTENT(IN) :: wave, flux, sigma
REAL, DIMENSION(5), INTENT(IN) :: guess
LOGICAL, DIMENSION(5), INTENT(IN) :: flag

REAL, DIMENSION(5), INTENT(OUT) :: fit, err
REAL, INTENT(OUT) :: chisq

fit = guess

CALL evfit(wave,flux,sigma,num,fit,flag,err,chisq)

END SUBROUTINE voigtfit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE voigtfit2(num,wave,flux,sigma,guess,flag,fit,err,chisq)
IMPLICIT NONE

INTEGER, INTENT(IN) :: num
REAL, DIMENSION(num), INTENT(IN) :: wave, flux, sigma
REAL, DIMENSION(6), INTENT(IN) :: guess
LOGICAL, DIMENSION(6), INTENT(IN) :: flag

REAL, DIMENSION(7), INTENT(OUT) :: fit, err
REAL, INTENT(OUT) :: chisq

LOGICAL, DIMENSION(7) :: newflag
REAL, DIMENSION(7) :: newguess

newguess(1:6) = guess(1:6)
newflag(1:6) = flag(1:6)

newguess(7) = 0
newflag(7) = .FALSE.

fit = newguess

CALL evfitm(wave,flux,sigma,num,fit,newflag,err,chisq)

END SUBROUTINE voigtfit2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE voigtguess(wave,flux,num,guess)
IMPLICIT NONE

INTEGER, INTENT(IN) :: num
REAL, DIMENSION(num), INTENT(IN) :: wave, flux
REAL, DIMENSION(5), INTENT(OUT) :: guess

CALL evinit(wave,flux,num,guess)

END SUBROUTINE voigtguess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE voigtguess2(wave,flux,num,guess)
IMPLICIT NONE

INTEGER, INTENT(IN) :: num
REAL, DIMENSION(num), INTENT(IN) :: wave, flux
REAL, DIMENSION(6), INTENT(OUT) :: guess

CALL evinit2(wave,flux,num,guess)

END SUBROUTINE voigtguess2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE voigtprof2(wave,a,npar,vgt,dvda)
IMPLICIT NONE

INTEGER, INTENT(IN) :: npar
REAL, INTENT(IN) :: wave
REAL, DIMENSION(npar), INTENT(IN) :: a
REAL, INTENT(OUT) :: vgt
REAL, DIMENSION(npar), INTENT(OUT) :: dvda

CALL evoigt(wave, a, npar, vgt, dvda)

END SUBROUTINE voigtprof2

