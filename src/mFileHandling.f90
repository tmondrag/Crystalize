! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 07 April 2017

! This is a more OOP version of my previous safe file handling module
! as suggested by Wil England to help protect the file descriptor integer
!
! functions:
!   filename = FILEobj%closeKeepScratchfile()
!   fileunit = FILEobj%getFunit()
! subroutines
!   CALL FILEobj%init()
!   CALL FILEobj%openReadOnly(filename)
!   CALL FILEobj%openReadWrite(filename)
!   CALL FILEobj%openWriteNew(filename)
!   CALL FILEobj%openWriteAppend(filename)
!   CALL FILEobj%openWriteReplace(filename)
!   CALL FILEobj%openScratchfile()
!   CALL FILEobj%close()
!   CALL FILEobj%closeDeleteFile()
!   CALL find_IU_info(fileunit)

MODULE mFileHandling
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  IMPLICIT NONE

  INTEGER, PARAMETER  :: stdout = OUTPUT_UNIT
  INTEGER, PARAMETER  :: stdin  = INPUT_UNIT
  INTEGER, PARAMETER  :: stderr = ERROR_UNIT

  TYPE,PUBLIC :: FILE
    INTEGER, PRIVATE      :: fd ! file unit or file descriptor number
    CHARACTER(LEN=256)    :: filename

  CONTAINS
    PROCEDURE :: init                     => init_FILE
    PROCEDURE :: openReadOnly
    PROCEDURE :: openReadWrite
    PROCEDURE :: openWriteNew
    PROCEDURE :: openWriteAppend
    PROCEDURE :: openWriteReplace
    PROCEDURE :: openScratchfile
    PROCEDURE :: close                    => closeDefault
    PROCEDURE :: closeDeleteFile
    PROCEDURE :: closeKeepScratchfile

    PROCEDURE :: getFunit
  END TYPE FILE

  PUBLIC  ::  openReadOnly, openReadWrite, openWriteNew, openWriteAppend, openScratchfile
  PUBLIC  ::  openWriteReplace, closeDeleteFile, closeKeepScratchfile, getFunit
  PRIVATE ::  init_FILE,closeDefault
CONTAINS
  ! Initalize FILE object to a safe, well known state.
  ! USAGE: CALL this%init()
  SUBROUTINE init_FILE(this)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this

    this%fd = stderr
  END SUBROUTINE init_FILE

  ! Open a file in read only mode
  ! USAGE: CALL this%openReadOnly(filename)
  SUBROUTINE openReadOnly(this,filename)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    CHARACTER(LEN=*),INTENT(IN)                   :: filename
    CHARACTER(LEN=256)                            :: io_message = ""
    INTEGER                                       :: fd
    INTEGER                                       :: io_err
    LOGICAL                                       :: opened, exists

    INQUIRE(FILE=TRIM(filename), EXIST=exists, OPENED=opened)
    IF(.NOT. exists) THEN
      WRITE(stderr,'(A,A,A)') "ERROR: Could not open ",TRIM(filename)," for reading; file doesn't exist."
      STOP
    ELSE IF (opened) THEN
      WRITE(stderr,'(A,A,A)') "ERROR: Refused to open ",TRIM(filename)," for reading; file already opened."
      STOP
    END IF

    OPEN(NEWUNIT=fd, FILE=TRIM(filename), ACTION='READ', IOSTAT=io_err, IOMSG=io_message)
    IF(io_err .NE. 0) THEN
      WRITE(stderr,'(A)') io_message
      STOP
    END IF

    this%fd = fd
    this%filename = filename
  END SUBROUTINE openReadOnly

  ! Open a file in readwrite mode
  ! USAGE: CALL this%openReadWrite(filename)
  SUBROUTINE openReadWrite(this,filename)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                      :: this
    CHARACTER(LEN=*),INTENT(IN)                   :: filename
    INTEGER                                       :: fd
    CHARACTER(LEN=256)                            :: io_message = ""
    INTEGER                                       :: io_err
    LOGICAL                                       :: opened, exists

    INQUIRE(FILE=TRIM(filename), EXIST=exists, OPENED=opened)
    IF(.NOT. exists) THEN
      WRITE(stderr,'(A,A,A)') "Warning: ",TRIM(filename)," doesn't exist; new empty file will be created."
    ELSE IF (opened) THEN
      WRITE(stderr,'(A,A,A)') "ERROR: Refused to open ",TRIM(filename),"; file already opened."
      STOP
    END IF

    OPEN(NEWUNIT=fd, FILE=TRIM(filename), ACTION='READWRITE', IOSTAT=io_err, IOMSG=io_message)
    IF(io_err .NE. 0) THEN
      WRITE(stderr,'(A)') io_message
      STOP
    END IF

    this%fd = fd
    this%filename = filename
  END SUBROUTINE openReadWrite

  ! Open a new file for writing
  ! USAGE: CALL this%openWriteNew(filename)
  SUBROUTINE openWriteNew(this,filename)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    CHARACTER(LEN=*),INTENT(IN)                   :: filename
    INTEGER                                       :: fd
    CHARACTER(LEN=256)                            :: io_message = ""
    INTEGER                                       :: io_err
    LOGICAL                                       :: opened, exists

    INQUIRE(FILE=TRIM(filename), EXIST=exists, OPENED=opened)
    IF(exists) THEN
      WRITE(stderr,'(A,A,A)') "ERROR: ",TRIM(filename)," already exists; refused to open in case of overwrite."
      STOP
    ELSE IF (opened) THEN
      WRITE(stderr,'(A,A,A)') "ERROR: ",TRIM(filename)," doesn't exist but is already open;", &
                              "something is seriously wrong with your OS."
      STOP
    END IF

    OPEN(NEWUNIT=fd, FILE=TRIM(filename), ACTION='WRITE', STATUS='NEW', IOSTAT=io_err, IOMSG=io_message)
    IF(io_err .NE. 0) THEN
      WRITE(stderr,'(A)') io_message
      STOP
    END IF

    this%fd = fd
    this%filename = filename
  END SUBROUTINE openWriteNew

  ! Open a file for writing, appending new content to the file's previous content
  ! USAGE: CALL this%openWriteAppend(filename)
  SUBROUTINE openWriteAppend(this,filename)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    CHARACTER(LEN=*),INTENT(IN)                   :: filename
    INTEGER                                       :: fd
    CHARACTER(LEN=256)                            :: io_message = ""
    INTEGER                                       :: io_err
    LOGICAL                                       :: opened, exists

    INQUIRE(FILE=TRIM(filename), EXIST=exists, OPENED=opened)
    IF(.NOT.exists) THEN
      WRITE(stderr,'(A,A,A)') "Warning: ",TRIM(filename)," doesn't exist; new empty file will be created."
    ELSE IF (opened) THEN
      WRITE(stderr,'(A,A,A)') "ERROR: ",TRIM(filename)," is already open; refused to open to prevent ", &
                              "unintentional writing."
      STOP
    END IF

    IF (exists) THEN
      OPEN(NEWUNIT=fd, FILE=TRIM(filename), STATUS="old", POSITION="append", ACTION="write", IOSTAT=io_err, &
           IOMSG=io_message)
      IF(io_err .NE. 0) THEN
          WRITE(stderr,'(A)') io_message
          STOP
      END IF
    ELSE
      OPEN(NEWUNIT=fd, FILE=TRIM(filename), STATUS="new", ACTION="write", IOSTAT=io_err, IOMSG=io_message)
      IF(io_err .NE. 0) THEN
          WRITE(stderr,'(A)') io_message
      STOP
      END IF
    END IF

    this%fd = fd
    this%filename = filename
  END SUBROUTINE openWriteAppend

  ! Open a file for writing, overwriting its previous content
  ! USAGE: CALL this%openWriteReplace(filename)
  SUBROUTINE openWriteReplace(this,filename)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    INTEGER                                       :: fd
    CHARACTER(LEN=*),INTENT(IN)                   :: filename
    CHARACTER(LEN=256)                            :: io_message = ""
    INTEGER                                       :: io_err
    LOGICAL                                       :: opened, exists

    INQUIRE(FILE=TRIM(filename), EXIST=exists, OPENED=opened)
    IF(.NOT.exists) THEN
      WRITE(stderr,'(A,A,A)') "Warning: ",TRIM(filename)," doesn't exist; new empty file will be created."
    ELSE IF (opened) THEN
      WRITE(stderr,'(A,A,A)') "ERROR: ",TRIM(filename)," is already open; refused to open to prevent ", &
                              "unintentional writing."
      STOP
    END IF

    OPEN(NEWUNIT=fd, FILE=TRIM(filename), STATUS="REPLACE", ACTION="write", IOSTAT=io_err, IOMSG=io_message)
    IF(io_err .NE. 0) THEN
      WRITE(stderr,'(A)') io_message
      STOP
    END IF

    this%fd = fd
    this%filename = filename
  END SUBROUTINE openWriteReplace

  ! Open a scratch file, which will be deleted upon closing by default
  ! USAGE: fd = safeopen_scratchfile()
  SUBROUTINE openScratchfile(this)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    INTEGER                                       :: fd
    CHARACTER(LEN=256)                            :: io_message = "",filename
    INTEGER                                       :: io_err

    OPEN(NEWUNIT=fd, STATUS="SCRATCH", ACTION="READWRITE", IOSTAT=io_err, IOMSG=io_message)
    IF(io_err .NE. 0) THEN
        WRITE(stderr,'(A)') io_message
        STOP
    END IF

    INQUIRE(UNIT=fd,NAME=filename)
    this%fd = fd
    this%filename = filename
  END SUBROUTINE openScratchfile

  ! Default file close. Will close and delete only if the file is a scratchfile.
  ! USAGE: CALL this%close()
  SUBROUTINE closeDefault(this)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    LOGICAL                                       :: opened
    LOGICAL                                       :: exists
    CHARACTER(LEN=256)                            :: filename

    INQUIRE(UNIT=this%fd, OPENED=opened, EXIST=exists, NAME=filename)
    IF (.NOT.opened) THEN
      IF (exists) THEN
        WRITE(stderr,'(A,A,A)') "ERROR: Closing file ",TRIM(filename),". File is not open."
        this%filename = filename
        STOP
      ELSE
        WRITE(stderr,'(A,i4,A)') "ERROR:  Closing file descriptor unit ",this%fd,". File does not exist."
        STOP
      END IF
    END IF
    CLOSE(UNIT=this%fd)
  END SUBROUTINE closeDefault

  ! Special close for the special occasion when a file should be deleted after being closed
  ! don't bother using for files opened with the status SCRATCH; fortran deletes those on close already
  ! USAGE: CALL closeDeleteFilelete()
  SUBROUTINE closeDeleteFile(this)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    LOGICAL                                       :: opened
    LOGICAL                                       :: exists
    CHARACTER(LEN=256)                            :: filename

    INQUIRE(UNIT=this%fd, OPENED=opened, EXIST=exists, NAME=filename)
    IF (.NOT.opened) THEN
      IF (exists) THEN
        WRITE(stderr,'(A,A,A)') "ERROR: Closing file ",TRIM(filename),". File is not open; refused to delete ", &
                                 "to prevent unintentional deletion."
        this%filename = filename
        STOP
      ELSE
        WRITE(stderr,'(A,i4,A)') "ERROR: Closing file descriptor unit ",this%fd,". File does not exist; refused to delete ", &
                                 "to prevent unintentional deletion."
        STOP
      END IF
    END IF
    CLOSE(UNIT=this%fd,STATUS='DELETE')
  END SUBROUTINE closeDeleteFile

  ! Special close function for the rare occasion you want to keep a scratchfile
  ! USAGE: newfilename = closeKeepScratchfile()
  FUNCTION closeKeepScratchfile(this) RESULT(filename)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    CHARACTER(LEN=256)                            :: filename
    LOGICAL                                       :: named, opened

    INQUIRE(UNIT=this%fd, NAMED=named, OPENED=opened)
    IF (opened) THEN
      IF (named) THEN
        INQUIRE(UNIT=this%fd,NAME=filename)
        WRITE(stdout,'(A,i6,A,A,A)') "INFO: Closing scatchfile descriptor unit ",this%fd," will be named ",TRIM(filename),"."
        this%filename = filename
      ELSE
        WRITE(stderr,'(A,i6,A)') "ERROR: Closing scratchfile descriptor unit ",this%fd," is unnamed. Naughty computer!"
        STOP
      END IF
        CLOSE(UNIT=this%fd,STATUS="KEEP")
    ELSE
      WRITE(stderr,'(A,i6,A)') "ERROR: Closing scratchfile descriptor unit ",this%fd," does not exist. Naughty programmer!"
      STOP
    END IF
  END FUNCTION closeKeepScratchfile

  ! A getter for the file's fileunit integer
  ! USAGE: funit = this%getFunit()
  FUNCTION getFunit(this) RESULT(fd)
    IMPLICIT NONE
    CLASS(FILE),INTENT(INOUT)                     :: this
    INTEGER                                       :: fd

    fd = this%fd
  END FUNCTION getFunit

  ! This is just an informative printer of information about a given IO file descriptor unit
  ! USAGE: CALL find_IU_info(fd)
  SUBROUTINE find_IU_info(fd)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                            :: fd
    CHARACTER(LEN=256)                            :: filename, omode, pos
    LOGICAL                                       :: opened, exists, named

    INQUIRE(UNIT=fd, OPENED=opened, EXIST=exists, NAMED=named, ACTION=omode, POSITION=pos)
    IF (opened) THEN
      WRITE(stdout,'(A,i6,A,A,A,A,A)') "INFO: File descriptor unit ",fd," is opened &
                                       &in mode ",TRIM(omode), " at position ",TRIM(pos),"."
    ELSE
      WRITE(stdout,'(A,i6,A)') "INFO: File descriptor unit ",fd," is not open."
    END IF
    IF (exists) THEN
      WRITE(stdout,'(A,i6,A)') "INFO: File descriptor unit ",fd," is in the range of values &
                                 &allowed by the compiler."
    ELSE
      WRITE(stdout,'(A,i6,A)') "INFO: File descriptor unit ",fd," not allowed by the compiler."
    END IF
    IF (named) THEN
      INQUIRE(UNIT=fd,NAME=filename)
      WRITE(stdout,'(A,i6,A,A,A)') "INFO: File descriptor unit ",fd," is named ",TRIM(filename),"."
    ELSE
      WRITE(stdout,'(A,i6,A)') "INFO: File descriptor unit ",fd," is unnamed."
    END IF
  END SUBROUTINE find_IU_info
END MODULE mFileHandling
