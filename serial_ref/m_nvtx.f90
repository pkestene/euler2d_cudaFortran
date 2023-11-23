#define nop(x) associate( x => x ); end associate

module m_nvtx

  use, intrinsic :: iso_c_binding
  implicit none

  integer,private,parameter :: nbcol=7
  integer(kind=C_INT32_T),private :: col(nbcol) = [ &
       & int(Z'0000ff00',kind=C_INT32_T), &
       & int(Z'000000ff',kind=C_INT32_T), &
       & int(Z'00ffff00',kind=C_INT32_T), &
       & int(Z'00ff00ff',kind=C_INT32_T), &
       & int(Z'0000ffff',kind=C_INT32_T), &
       & int(Z'00ff0000',kind=C_INT32_T), &
       & int(Z'00ffffff',kind=C_INT32_T) ]

  character,private,target :: tempName(256)

  type, bind(C):: nvtxEventAttributes
     integer(C_INT16_T):: version=1
     integer(C_INT16_T):: size=48 !
     integer(C_INT):: category=0
     integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT):: color
     integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT):: reserved0
     integer(C_INT64_T):: payload   ! union uint,int,double
     integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1
     type(C_PTR):: message  ! ascii char
  end type nvtxEventAttributes

  interface nvtxRangePush
     ! push range with custom label and standard color
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=C_CHAR) :: name(256)
     end subroutine nvtxRangePushA

     ! push range with custom label and custom color
     subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
       use iso_c_binding
       import:: nvtxEventAttributes
       type(nvtxEventAttributes):: event
     end subroutine nvtxRangePushEx
  end interface nvtxRangePush

  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop

contains

  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
#ifdef EULER2D_NVTX_ANNOTATION_ENABLED
    type(nvtxEventAttributes):: event
    character(kind=c_char,len=256) :: trimmed_name
    integer:: i

    trimmed_name=trim(name)//c_null_char

    ! move scalar trimmed_name into character array tempName
    do i=1,LEN(trim(name)) + 1
       tempName(i) = trimmed_name(i:i)
    enddo


    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
#else
    nop(name)
    nop(id)
    nop(col)
    nop(tempname)
#endif
  end subroutine nvtxStartRange

  subroutine nvtxEndRange
#ifdef EULER2D_NVTX_ANNOTATION_ENABLED
    call nvtxRangePop
#endif
  end subroutine nvtxEndRange

end module m_nvtx
