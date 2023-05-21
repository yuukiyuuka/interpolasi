program readmatrix

implicit none
integer, parameter :: dpr=kind(1.0D0)
integer :: i,n,k,j,row,col
real(dpr), allocatable :: D(:,:),X(:),Y(:),C(:,:),B(:),A(:)

!jumlah variabel yang tidak diketahui
write(*,*) 'masukkan ordo polinomial dan jumlah data'
read(*,*) n,row

allocate(D(1:row,1:2),X(1:row),Y(1:row),C(1:n,1:n),B(1:n),A(1:n))

!baca data matriks
open(unit=1,file='input_matrix.txt',action='read')
do i=1,row
    read(1,*,end=2) D(i,:)
    X(i)=D(i,1)
    Y(i)=D(i,2)
end do

2 close(1)

!tampilkan data X
write(*,*) 'Data X'
do i=1,n
    write(*,*) X(i)
end do

!tampilkan data Y
write(*,*) 'Matriks Y'
do i=1,n
    write(*,*) Y(i)
end do

!Inisialisasi matriks C dengan nilai 0
do i=1,n
   C(i,:) = 0.0
end do

!menghitung matriks C
do k=1,n
   do j=1,n 
      do i=1,row
         C(k,j)=C(k,j)+X(i)**(k+j-2)
      end do
   end do
end do

!tampilkan matriks C
write(*,*) 'Matriks C'
do i=1,n
    write(*,*) (C(i,j),j=1,n)
end do

!menghitung matriks B
do i=1,row 
   do k=1,n
      B(k)=Y(i)*X(i)**(k-1)
   end do
end do

!tampilkan matriks B
write(*,*) 'Matriks B'
do i=1,n
   write(*,*) B(i)
end do

!mencari matriks koefisien A
A=lu_solve(C,B,n)

!tampilkan matriks A
write(*,*) 'Matriks A'
do i=1,n
   write(*,*) A(i)
end do

!save matriks C,B, dan A ke file eksternal
write(*,*) 'cek file eksternal'
open(unit=3,file='matriks_CB.txt',status='replace',action='write')
write(3,*) 'Matriks C'
do i=1,n
   write(3,*) (C(i,j),j=1,n)
end do
write(3,*) 'Matriks B'
do i=1,n
   write(3,*) B(i)
end do
write(3,*) 'Matriks koefisien A'
do i=1,n
   write(3,*) A(i)
end do
close(3)

deallocate(D,X,Y,C,B)

stop

contains

function lu_solve(AA, BB, n) result(X)
  implicit none
  integer, intent(in) :: n
  real(dpr), intent(in) :: AA(n,n)
  real(dpr), intent(in) :: BB(n)
  real(dpr) :: X(n)
  integer :: i, j, k, maxind
  real(dpr) :: maxval, temp

  !copy AA and BB to temporary arrays
  real(dpr) :: A(n,n), B(n)
  A = AA
  B = BB

  ! Perform LU decomposition with partial pivoting
  do k = 1, n-1
    ! Find pivot element
    maxind = k
    maxval = abs(A(k,k))
    do i = k+1, n
      if (abs(A(i,k)) > maxval) then
        maxind = i
        maxval = abs(A(i,k))
      end if
    end do
    ! Swap rows if necessary
    if (maxind /= k) then
      do j = 1, n
        temp = A(k,j)
        A(k,j) = A(maxind,j)
        A(maxind,j) = temp
      end do
      temp = B(k)
      B(k) = B(maxind)
      B(maxind) = temp
    end if
    ! Perform elimination
    do i = k+1, n
      temp = A(i,k) / A(k,k)
      A(i,k) = temp
      do j = k+1, n
        A(i,j) = A(i,j) - temp * A(k,j)
      end do
      B(i) = B(i) - temp * B(k)
    end do
  end do

  ! Solve the system by back substitution
  X(n) = B(n) / A(n,n)
  do i = n-1, 1, -1
    X(i) = B(i)
    do j = i+1, n
      X(i) = X(i) - A(i,j) * X(j)
    end do
    X(i) = X(i) / A(i,i)
  end do

end function lu_solve

end program readmatrix