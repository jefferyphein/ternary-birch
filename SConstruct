# scons makefile.
# To compile, simply run "scons" in the same directory as this file.
# If compiling on windows, you may need to modify the libpath list so that you
# can properly load external libraries.
#

env = Environment()
debug = ARGUMENTS.get('debug', 0)

srcs = Split('''
    main.cpp
    AutomorphismZZ.cpp
    CharacterZZ.cpp
    Eigenvector.cpp
    GenusZZ.cpp
    IsometryZZ.cpp
    MathZZ.cpp
    NeighborIteratorZZ.cpp
    QuadFormZZ.cpp
    SetCover.cpp
    SparseMatrix.cpp
''')

libs = Split('''
    gmpxx
    gmp
    m
    pthread
''')

libpath = Split('''
    .
    /usr/lib
    /usr/local/lib
''')

ccflags = Split('''
    -Wall
    -O3
    -Wextra
    -Werror
    -pedantic
    -std=c++11
''')

if int(debug):
    ccflags += Split('''
        -g
        -DDEBUG
    ''')

Program(
    target = 'birch',
    source = srcs,
    LIBS = libs,
    LIBPATH = libpath,
    CCFLAGS = ccflags
)
