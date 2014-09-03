module testParser_mod
  use pFUnit
  use Parser_mod
  implicit none
  private

  ! legacy
  public :: testStripComment_noComment
  public :: testStripComment
  public :: testStripComment_1stChar

  public :: testStripComment_new
  public :: testSkipHeader
  public :: testNoBeginData
  public :: testIsEndData

  public :: testGetValueType
  public :: testGetCommonValueType

  public :: testSplitTokensA
  public :: testSplitTokensB
  public :: testSplitTokensC
  public :: testSplitTokensD
  public :: testSplitEmbeddedComma
  public :: testBadFirstSeparator

  public :: testParseEOF
  public :: testParseSimple
  public :: testParse
  public :: testParseNoValue

  public :: testWriteFormatted
  public :: testWriteTextExample


!TO DO - create tests for these:
!!$  public :: testMissingEndQuote
!!$  public :: testFirstSep ! '='
!!$  public :: testIgnoreEmbeddedComments

contains
  subroutine testStripComment_noComment()
    character(len=30) :: expectedString
    character(len=30) :: string

    string = 'has no comment characters'
    expectedString = string
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment_noComment

  subroutine testStripComment()
    character(len=30) :: expectedString
    character(len=30) :: string

    string = 'comment starts here ! ignore this'
    expectedString = 'comment starts here'
    
    call assertEqual(expectedString, strip_comment(string))

    string = 'comment starts #'
    expectedString = 'comment starts'
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment

  subroutine testStripComment_1stChar()
    character(len=30) :: expectedString
    character(len=30) :: string

    string = '! all comment'
    expectedString = ''
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment_1stChar

  subroutine testStripComment_new()
    type (Parser_type) :: parser
    character(len=40) :: expectedString
    character(len=40) :: string

    call setCommentCharacters(parser, ':')
    string = 'has no comment characters? ! comment'
    expectedString = string
    call assertEqual(expectedString, stripComment(parser, string))

    call setCommentCharacters(parser, ' ,')
    expectedString = 'has'
    call assertEqual(expectedString, stripComment(parser, string))

    call setCommentCharacters(parser, '#!')
    expectedString = 'has no comment characters?'
    call assertEqual(expectedString, stripComment(parser, string))

  end subroutine testStripComment_new

  subroutine testSkipHeader()
    use FileManager
    type (Parser_type) :: parser
    integer :: unit
    integer :: status
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=256) :: line

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) 'a'
    write(unit,*) 'b'
    write(unit,*) BEGIN_DATA
    write(unit,*) 'first line after header'
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call skipHeader(parser, unit, status)
    read(unit,'(a257)') line
    call assertEqual('first line after header', line)

    close(unit, status='delete')

  end subroutine testSkipHeader

  subroutine testNoBeginData()
    use FileManager
    type (Parser_type) :: parser
    integer :: unit
    integer :: status
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=256) :: line

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) 'a'
    write(unit,*) 'b'
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call skipHeader(parser, unit, status)
    call assertTrue(status /= 0)

    close(unit, status='delete')
  end subroutine testNoBeginData

  subroutine testGetValueType()
    use Dictionary_mod

    call assertEqual(INTEGER_TYPE, getValueType('3'), 'case A')
    call assertEqual(INTEGER_TYPE, getValueType('123456'), 'case B')

    call assertEqual(REAL64_TYPE, getValueType('1.'), 'case C')
    call assertEqual(REAL64_TYPE, getValueType('1.2'), 'case D')
    call assertEqual(REAL64_TYPE, getValueType('1.2e+10'), 'case E')
    call assertEqual(REAL64_TYPE, getValueType('234.6d-12'), 'case F')

    call assertEqual(LOGICAL_TYPE, getValueType('.true.'), 'case G')
    call assertEqual(LOGICAL_TYPE, getValueType('.false.'), 'case H')
    call assertEqual(LOGICAL_TYPE, getValueType('T'), 'case I')
    call assertEqual(LOGICAL_TYPE, getValueType('F'), 'case J')
    call assertEqual(LOGICAL_TYPE, getValueType('t'), 'case K')
    call assertEqual(LOGICAL_TYPE, getValueType('f'), 'case L')
    call assertEqual(LOGICAL_TYPE, getValueType('t'), 'case M')
    call assertEqual(LOGICAL_TYPE, getValueType('false'), 'case N')
    call assertEqual(LOGICAL_TYPE, getValueType('true'), 'case O')

    call assertEqual(STRING_TYPE, getValueType('hello'), 'case P')
    call assertEqual(STRING_TYPE, getValueType('truely'), 'case Q')
    call assertEqual(STRING_TYPE, getValueType('Fa'), 'case R')
    call assertEqual(STRING_TYPE, getValueType('''.true.'''), 'case S') ! has quotes => string
    call assertEqual(STRING_TYPE, getValueType('".true."'), 'case T') ! has quotes => string

  end subroutine testGetValueType

  subroutine testGetCommonValueType()
!@sum Test that a mixed set of tokens are treated as a category that
!@+ sufficiently broad.  E.g mix of integer and string should be
!@+ treated as string.
    use Dictionary_mod

    call assertEqual(INTEGER_TYPE, getValueType(['3','4']))
    call assertEqual(REAL64_TYPE, getValueType(['3.','4 ']))
    call assertEqual(REAL64_TYPE, getValueType(['3 ','4.', '5 ']))
    call assertEqual(REAL64_TYPE, getValueType(['3 ','4.']))
    call assertEqual(REAL64_TYPE, getValueType(['3.','4.']))

    call assertEqual(LOGICAL_TYPE, getValueType(['T     ','F     ','.true.']))
    call assertEqual(STRING_TYPE,  getValueType(['3     ','F     ','.true.']))

  end subroutine testGetCommonValueType



  subroutine testIsEndData()
    use FileManager
    type (Parser_type) :: parser
    character(len=*), parameter :: END_DATA = '*** end params ***'

    call setEndData(parser, END_DATA)
    call assertTrue(isEndData(parser, END_DATA))
    call assertFalse(isEndData(parser, 'some other string'))
  end subroutine testIsEndData

  subroutine testSplitTokensA()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, '  ')
    call assertEqual(0, size(tokens))
    deallocate(tokens)

  end subroutine testSplitTokensA

  subroutine testSplitTokensB()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, '  hello  ')
    call assertEqual(1, size(tokens))
    if (size(tokens) == 1) call assertEqual('hello', tokens(1))
    deallocate(tokens)

  end subroutine testSplitTokensB

  subroutine testSplitTokensC()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, ' var = 3  ')
    call assertEqual(2, size(tokens))
    if (size(tokens) == 2) then
      call assertEqual('var', tokens(1))
      call assertEqual('3', tokens(2))
    end if
    deallocate(tokens)

  end subroutine testSplitTokensC

  subroutine testSplitTokensD()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, ' var = "hello", 3, 7.891e+10, .true., 123456  ')
    call assertEqual(6, size(tokens))
    if (size(tokens) == 6) then
      call assertEqual('var', tokens(1))
      call assertEqual('"hello"', tokens(2))
      call assertEqual('3', tokens(3))
      call assertEqual('7.891e+10', tokens(4))
      call assertEqual('.true.', tokens(5))
      call assertEqual('123456', tokens(6))
    end if
    deallocate(tokens)

  end subroutine testSplitTokensD

  subroutine testSplitEmbeddedComma()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, " var = 'short name', 'name with , in it'")
    call assertEqual(3, size(tokens))
    if (size(tokens) == 3) then
      call assertEqual('var', tokens(1))
      call assertEqual("'short name'", tokens(2))
      call assertEqual("'name with , in it'", tokens(3))
    end if
    deallocate(tokens)

  end subroutine testSplitEmbeddedComma

  subroutine testBadFirstSeparator()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '= ,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, " my var = 3")
    call assertFailedAssert('Parser_mod: Illegal syntax.  "=" not first separator.', &
         & 'Failed to detect illegal syntax.')
    deallocate(tokens)

  end subroutine testBadFirstSeparator

  subroutine testParseEOF()
    use FileManager
    use Dictionary_mod
    type (Parser_type) :: parser
    type (Dictionary_type) :: aDictionary

    integer :: status
    integer :: unit

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    rewind(unit)

    aDictionary = parse(parser, unit, status)
    call assertTrue(status /= 0) ! EOF
    close(unit, status='delete')

  end subroutine testParseEOF

  subroutine testParseSimple()
    use FileManager
    use Dictionary_mod
    use GenericType_mod
    type (Parser_type) :: parser

    type (Dictionary_type) :: aDictionary
    integer :: unit
    integer :: i(3)
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'
    integer :: status

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) BEGIN_DATA
    write(unit,*) 'key = 1, 2, 3'
    write(unit,*) END_DATA
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aDictionary = parse(parser, unit, status)
    close(unit, status='delete')

    call assertEqual(1, size(getKeys(aDictionary)))
    call assertTrue(hasKey(aDictionary, 'key'), 'There should not be a "key".')
    i = lookup(aDictionary, 'key')
    call assertEqual([1,2,3], i, 'Wrong value for "key".')

    call clean(aDictionary)

  end subroutine testParseSimple

  subroutine testParse()
    use FileManager
    use Dictionary_mod
    use GenericType_mod
    type (Parser_type) :: parser

    type (Dictionary_type) :: aDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'

    integer :: status

    integer :: betaValue
    real(kind=r64) :: gammaValue
    logical :: deltaValues(4)

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) 'alpha = 0' ! ignore - in header
    write(unit,*) 'b'
    write(unit,*) BEGIN_DATA
    write(unit,*) '! alpha = 1' ! ignore - only a comment
    write(unit,*) 'beta = 2'
    write(unit,*) ! empty line
    write(unit,*) '!' ! effectively empty line
    write(unit,*) 'gamma = 3.'
    write(unit,*) '! alpha = 2' ! ignore - only a comment
    write(unit,*) 'delta = T, F, T, T'
    write(unit,*) END_DATA
    write(unit,*) '! alpha = 3' ! ignore - after header
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aDictionary = parse(parser, unit, status)
    close(unit, status='delete')

    call assertEqual(3, size(getKeys(aDictionary)))
    call assertFalse(hasKey(aDictionary, 'alpha'), 'There should not be an "alpha".')
    call assertTrue(hasKey(aDictionary, 'beta'), 'There should be a "beta".')
    call assertTrue(hasKey(aDictionary, 'gamma'), 'There should be a "gamma".')
    call assertTrue(hasKey(aDictionary, 'delta'), 'There should be a "delta".')

    betaValue = lookup(aDictionary, 'beta')
    call assertEqual(2, betaValue, 'Wrong value for "beta".')

    gammaValue = lookup(aDictionary, 'gamma')
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "gamma".')

    gammaValue = lookup(aDictionary, 'Gamma')
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "Gamma".')

    deltaValues = lookup(aDictionary, 'delta')
    call assertTrue(all([.true.,.false.,.true.,.true.] .eqv. deltaValues), 'Wrong value for "delta".')

    call clean(aDictionary)

  end subroutine testParse

  subroutine testParseNoValue()
!@sum Tests that an exception is thrown if there is no value to associate with a key.
    use FileManager
    use Dictionary_mod
    type (Parser_type) :: parser

    type (Dictionary_type) :: aDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'
    character(len=256) :: line

    integer :: status

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) BEGIN_DATA
    write(unit,*) 'beta'
    write(unit,*) END_DATA
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, ' =,')
    call setCommentCharacters(parser, '!#')
    aDictionary = parse(parser, unit, status)
    call assertFailedAssert('Parser_mod: syntax error in input unit.', &
         & 'Failed to detect syntax error in parse')

    close(unit, status='delete')
  end subroutine testParseNoValue

  subroutine testWriteFormatted()
    use FileManager
    use Dictionary_mod
    use GenericType_mod
    type (Parser_type) :: parser

    type (Dictionary_type) :: aDictionary
    type (Dictionary_type) :: bDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'
    character(len=256) :: line
    character(len=MAX_LEN_KEY), pointer :: keys(:)

    integer :: status

    integer :: i(3)

    aDictionary = Dictionary()
    call insert(aDictionary,'key',[1,2,3])

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')
    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    call writeFormatted(parser, unit, aDictionary)
    close(unit)

    call openUnit('testParser.txt', unit, qold=.true., qbin=.false.)

    bDictionary = parse(parser, unit, status)

    keys => getKeys(bDictionary)
    i = lookup(bDictionary, 'key')

    call assertTrue(aDictionary == bDictionary)

    call clean(bDictionary)
    call clean(aDictionary)

    close(unit, status='delete')

  end subroutine testWriteFormatted

  subroutine testWriteTextExample()
    use FileManager
    use Dictionary_mod
    type (Parser_type) :: parser

    type (Dictionary_type) :: aDictionary
    type (Dictionary_type) :: bDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '{'
    character(len=*), parameter :: END_DATA = '}'
    character(len=256) :: line
    character(len=MAX_LEN_KEY), pointer :: keys(:)

    integer :: status

    integer :: scale
    real*8 :: mass
    character(len=MAX_LEN_TOKEN) :: name

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')
    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)

    aDictionary = Dictionary()
    call insert(aDictionary,'name', 'Air')
    call insert(aDictionary,'logScale', -2)
    call insert(aDictionary,'molecularMass', 28.9655d0)
    call writeFormatted(parser, unit, aDictionary)
    
    rewind(unit)
    bDictionary = parse(parser, unit, status)

    keys => getKeys(bDictionary)
    call assertTrue(aDictionary == bDictionary, 'Dictionaries differ.')

    call clean(bDictionary)
    call clean(aDictionary)

    close(unit, status='delete')

  end subroutine testWriteTextExample

end module testParser_mod
