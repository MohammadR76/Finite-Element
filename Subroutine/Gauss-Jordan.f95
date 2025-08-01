!===================================================================================
! Subroutine: GAUSS
! Author    : Mohammadreza seyfpour
! Email		: Mohammadrezaseyfpour@gmail.com
! Description:
!   Solves a system of linear equations Ax = b using Gauss elimination
!   with both forward and backward elimination (Gauss-Jordan style).
!
! INPUT:
!   N   - Number of equations/unknowns
!   A   - Augmented matrix of size (N, N+1) where last column is vector b
!
! OUTPUT:
!   B   - Solution vector of size (N, 1)
!===================================================================================

SUBROUTINE GAUSS(N, A, B)
    IMPLICIT NONE

    ! ---------- Declarations ----------
    INTEGER, INTENT(IN) :: N				   ! N is Dimension
    REAL, DIMENSION(N, N+1) :: A               ! Augmented matrix [A | b]
    REAL, DIMENSION(N, 1), INTENT(OUT) :: B    ! Solution vector
    INTEGER :: I, J, K                         ! Loop counters
    REAL :: FACTOR                             ! Multiplier for elimination

    ! ---------- Forward Elimination ----------
    DO K = 1, N-1
        WRITE(*,*) 'Forward step: K =', K
        DO I = K+1, N
            IF (A(K,K) == 0.0) CYCLE           ! Avoid division by zero
            IF (K == I) CYCLE                  ! Skip diagonal
            FACTOR = -A(I,K) / A(K,K)
            DO J = 1, N+1
                A(I,J) = FACTOR * A(K,J) + A(I,J)
            END DO
        END DO
    END DO

    ! ---------- Backward Elimination ----------
    DO K = N, 1, -1
        WRITE(*,*) 'Backward step: K =', K
        DO I = N, 1, -1
            IF (A(K,K) == 0.0) CYCLE           ! Avoid division by zero
            IF (K == I) CYCLE                  ! Skip diagonal
            FACTOR = -A(I,K) / A(K,K)
            DO J = 1, N+1
                A(I,J) = FACTOR * A(K,J) + A(I,J)
            END DO
        END DO
    END DO

    ! ---------- Extract Solution ----------
    DO I = 1, N
        IF (A(I,I) == 0.0) THEN
            B(I,1) = 0.0    
        ELSE
            B(I,1) = A(I,N+1) / A(I,I)
        END IF
    END DO

    RETURN
END SUBROUTINE GAUSS
