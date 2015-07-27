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
SUBROUTINE voigtguess(wave,flux,num,guess)
IMPLICIT NONE

INTEGER, INTENT(IN) :: num
REAL, DIMENSION(num), INTENT(IN) :: wave, flux
REAL, DIMENSION(5), INTENT(OUT) :: guess

CALL evinit(wave,flux,num,guess)

END SUBROUTINE voigtguess
