module utils_kahan
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    public :: kahan_add

contains
    !> -------------------------------------------------------------------------
    !> Helper subroutine for Kahan Summation Algorithm (Compensated Summation)
    !> Reduces numerical error when adding a sequence of finite precision floating point numbers.
    !> -------------------------------------------------------------------------
    pure elemental subroutine kahan_add(sum, c, input)
        implicit none
        real(real64), intent(inout) :: sum ! Current sum
        real(real64), intent(inout) :: c ! Compensation (carry) variable
        real(real64), intent(in) :: input ! Value to add

        real(real64) :: y, t

        y = input - c ! So far, so good: c is zero.
        t = sum + y ! Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y ! (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
        sum = t ! Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    end subroutine kahan_add
end module utils_kahan
