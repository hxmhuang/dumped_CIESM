module cfio
implicit none

integer, parameter :: CFIO_PROC_CLIENT = 1
integer, parameter :: CFIO_PROC_SERVER = 2
integer, parameter :: CFIO_PROC_BLANK = 3

integer :: char_len = 128
integer :: array_len = 128

integer, parameter :: cfio_byte   = 1
integer, parameter :: cfio_char   = 2
integer, parameter :: cfio_short  = 3
integer, parameter :: cfio_int	  = 4
integer, parameter :: cfio_float  = 5
integer, parameter :: cfio_double = 6

interface cfio_put_att
    module procedure cfio_put_att_str
    module procedure cfio_put_att_int
    module procedure cfio_put_att_real
    module procedure cfio_put_att_double
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
    module procedure cfio_put_atts_int
    module procedure cfio_put_atts_real
    module procedure cfio_put_atts_double
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end interface

interface cfio_put_vara
    module procedure cfio_put_vara_real
    module procedure cfio_put_vara_double
    module procedure cfio_put_vara_int
    module procedure cfio_put_vara_text
end interface

contains

integer(4) function cfio_init(inter_comm, x_proc_num, y_proc_num, ratio)
    implicit none
    integer(4), intent(in) :: inter_comm, x_proc_num, y_proc_num, ratio

    call cfio_init_c(inter_comm, x_proc_num, y_proc_num, ratio, cfio_init)

end function

integer(4) function cfio_finalize()
    implicit none

    call cfio_finalize_c(cfio_finalize)

end function

integer(4) function cfio_proc_type()
    implicit none

    call cfio_proc_type_c(cfio_proc_type)

end function

integer(4) function cfio_create(path, cmode, ncid)
    implicit none
    character(len=*), intent(in) :: path
    integer(4) :: cmode, ncid, length

    length = len(trim(path))

    call cfio_create_c(trim(path), length, cmode, ncid, cfio_create)

end function

integer(4) function cfio_inq_dimid(ncid, name, dimid)
    implicit none
    integer(4), intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer(4), intent(in) ::  dimid 
    integer(4) name_length

    name_length = len(trim(name))

    call cfio_inq_dimid_c(ncid, trim(name), name_length, dimid, cfio_inq_dimid)

end function

integer(4) function cfio_inq_varid(ncid, name, varid)
    implicit none
    integer(4), intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer(4), intent(in) ::  varid 
    integer(4) name_length

    name_length = len(trim(name))

    call cfio_inq_varid_c(ncid, trim(name), name_length, varid, cfio_inq_varid)

end function

integer(4) function cfio_inq_att(ncid, varid, name)
    implicit none
    integer(4), intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer(4), intent(in) ::  varid 
    integer(4):: name_length

    name_length = len(trim(name))

    call cfio_inq_att_c(ncid, varid, trim(name), name_length, cfio_inq_att)

end function

integer(4) function cfio_inq_varndims(ncid, varid, ndims)
    implicit none
    integer(4), intent(in) :: ncid, varid, ndims

    call cfio_inq_varndims_c(ncid, varid, ndims, cfio_inq_varndims)

end function

integer(4) function cfio_inq_vardimid(ncid, varid, dimids)
    implicit none
    integer(4), intent(in) :: ncid, varid
    integer(4), dimension(*), intent(in) :: dimids
    call cfio_inq_vardimid_c(ncid, varid, dimids, cfio_inq_vardimid)
end function

integer(4) function cfio_def_dim(ncid, name, length, dimid)
    implicit none
    integer(4), intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer(4), intent(in) :: length, dimid 
    integer(4) name_length

    name_length = len(trim(name))

    call cfio_def_dim_c(ncid, trim(name), name_length, length, dimid, cfio_def_dim)

end function

integer(4) function cfio_def_var_old(ncid, name, xtype, ndims, dimids, start, &
	count, varid)
    implicit none
    integer(4), intent(in) :: ncid, xtype, ndims, varid
    character(len=*), intent(in) :: name
    integer(4), dimension(*), intent(in) :: dimids, start, count
    integer(4) name_length

    name_length = len(trim(name))

    call cfio_def_var_c(ncid, trim(name), name_length, xtype, ndims, dimids, &
	start, count, varid, cfio_def_var_old)

end function

integer(4) function cfio_def_var(ncid, name, xtype, ndims, dimids, varid)
    implicit none
    integer(4), intent(in) :: ncid, xtype, ndims, varid
    character(len=*), intent(in) :: name
    integer(4), dimension(*), intent(in) :: dimids
    integer(4) name_length
    integer(4), dimension(:), allocatable :: start, count
    integer(4) i

    !allocate(start(ndims))
    !allocate(count(ndims))
    !do i = 1, ndims
    !    start(i) = 0
    !    count(i) = 0
    !end do

    name_length = len(trim(name))

    call cfio_def_var_c(ncid, trim(name), name_length, xtype, ndims, dimids, &
	varid, cfio_def_var)

    !deallocate(start)
    !deallocate(count)

end function

integer function cfio_put_att_str(ncid, varid, name, values)
    implicit none
    integer(4), intent(in) :: ncid, varid
    character(len=*), intent(in) :: name, values
    integer(4) name_length
    
    name_length = len(trim(name))

    call cfio_put_att_c(ncid, varid, trim(name), name_length, cfio_char, len(values), values, &
	cfio_put_att_str)

end function

integer function cfio_put_att_int(ncid, varid, name, values)
    implicit none
    integer(4), intent(in) :: ncid, varid
    character(len=*), intent(in) :: name
    integer(4), intent(in) :: values
    integer(4) name_length
    
    name_length = len(trim(name))

    call cfio_put_att_c(ncid, varid, trim(name), name_length, cfio_int, 1, values, &
	cfio_put_att_int)

end function

integer function cfio_put_att_real(ncid, varid, name, values)
    implicit none
    integer(4), intent(in) :: ncid, varid
    character(len=*), intent(in) :: name
    real(4), intent(in) :: values
    integer(4) name_length
    
    name_length = len(trim(name))

    call cfio_put_att_c(ncid, varid, trim(name), name_length, cfio_float, 1, values, &
	cfio_put_att_real)

end function

integer function cfio_put_att_double(ncid, varid, name, values)
    implicit none
    integer(4), intent(in) :: ncid, varid
    character(len=*), intent(in) :: name
    real(8), intent(in) :: values
    integer(4) name_length
    
    name_length = len(trim(name))

    call cfio_put_att_c(ncid, varid, trim(name), name_length, cfio_double, 1, values, &
	cfio_put_att_double)

end function

integer function cfio_put_atts_int(ncid, varid, name, att_len, values)
    implicit none
    integer(4), intent(in) :: ncid, varid, att_len
    character(len=*), intent(in) :: name
    integer(4), intent(in) :: values(:)
    integer(4) name_length
    
    name_length = len(trim(name))

    call cfio_put_att_c(ncid, varid, trim(name), name_length, cfio_int, att_len, values, &
	cfio_put_atts_int)

end function

integer function cfio_put_atts_real(ncid, varid, name, att_len, values)
    implicit none
    integer(4), intent(in) :: ncid, varid, att_len
    character(len=*), intent(in) :: name
    real(4), intent(in) :: values(:)
    integer(4) name_length
    
    name_length = len(trim(name))

    call cfio_put_att_c(ncid, varid, trim(name), name_length, cfio_float, att_len, values, &
	cfio_put_atts_real)

end function

integer function cfio_put_atts_double(ncid, varid, name, att_len, values)
    implicit none
    integer(4), intent(in) :: ncid, varid, att_len
    character(len=*), intent(in) :: name
    real(8), intent(in) :: values(:)
    integer(4) name_length
    
    name_length = len(trim(name))

    call cfio_put_att_c(ncid, varid, trim(name), name_length, cfio_double, att_len, values, &
	cfio_put_atts_double)

end function

integer function cfio_enddef(ncid)
    implicit none
    integer(4), intent(in) :: ncid

    call cfio_enddef_c(ncid, cfio_enddef)

end function

integer function cfio_put_vara_real(ncid, varid, ndims, start, count, fp)
    implicit none
    integer(4), intent(in) :: ncid, varid, ndims
    integer(4), dimension(*), intent(in) :: start, count 
    real(4), dimension(*), intent(in) :: fp

    call cfio_put_vara_float_c(ncid, varid, ndims, start, count, fp, &
	cfio_put_vara_real)

end function

integer function cfio_put_vara_double(ncid, varid, ndims, start, count, fp)
    implicit none
    integer(4), intent(in) :: ncid, varid, ndims
    integer(4), dimension(*), intent(in) :: start, count 
    real(8), dimension(*), intent(in) :: fp

    call cfio_put_vara_double_c(ncid, varid, ndims, start, count, fp, &
	cfio_put_vara_double)

end function

integer function cfio_put_vara_int(ncid, varid, ndims, start, count, fp)
    implicit none
    integer(4), intent(in) :: ncid, varid, ndims
    integer(4), dimension(*), intent(in) :: start, count 
    integer(4), dimension(*), intent(in) :: fp

    call cfio_put_vara_int_c(ncid, varid, ndims, start, count, fp, &
	cfio_put_vara_int)
end function

integer function cfio_put_vara_text(ncid, varid, ndims, start, count, fp)
    implicit none
    integer(4), intent(in) :: ncid, varid, ndims
    integer(4), dimension(*), intent(in) :: start, count 
    character, intent(in) :: fp(:)

    call cfio_put_vara_text_c(ncid, varid, ndims, start, count, fp, &
	cfio_put_vara_text)
end function

integer function cfio_io_end()

    call cfio_io_end_c(cfio_io_end)

end function

integer function cfio_close(ncid)
    integer(4), intent(in) :: ncid

    call cfio_close_c(ncid, cfio_close)

end function

end module

