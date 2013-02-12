INCLUDES = $(GSL_CFLAGS)
AM_LDFLAGS  = 
AM_CFLAGS = -DUSE_LAPACK
AM_CPPFLAGS = -DUSE_LAPACK -Wno-deprecated

lib_LIBRARIES = libms.a
libms_a_SOURCES = abrand.cc glm.cc matrix.cc poly.cc timeseries.cc svm.cc garch.cc
# libms_a_CXXFLAGS = -ggdb -O0

EXTRA_DIST = *.h

