module mpi_context
    use mpi
    implicit none

    integer :: comm
    integer :: proc_id
    integer :: num_procs
    integer :: ierr

contains

    subroutine init_mpi()
      comm = MPI_COMM_WORLD
      call MPI_Init(ierr)
      call MPI_Comm_rank(comm, proc_id, ierr)
      call MPI_Comm_size(comm, num_procs, ierr)
    end subroutine init_mpi

    subroutine finalize_mpi()
      call MPI_Finalize(ierr)
    end subroutine finalize_mpi

    !add abort

end module mpi_context

