# scons makefile.
# To compile, simply run "scons" in the same directory as this file.
# If compiling on windows, you may need to modify the libpath list so that you
# can properly load external libraries.

srcs = Split('''
    main.cpp
    GenusZZ.cpp
    Ideal.cpp
    NeighborIterator.cpp
    QuadFormZZ.cpp
''')

libs = Split('''
    gmpxx
    gmp
    m
''')

libpath = Split('''
    .
    /usr/lib
    /usr/local/lib
''')

ccflags = Split('''
    -g
    -Wall
    -O3
    -Wextra
    -Werror
    -pedantic
    -std=c++11
''')

Program(
    target = 'birch',
    source = srcs,
    LIBS = libs,
    LIBPATH = libpath,
    CCFLAGS = ccflags
)
