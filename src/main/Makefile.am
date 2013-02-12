INCLUDES = -I$(top_srcdir)/src/mslib
AM_CPPFLAGS = $(GTKMM_CFLAGS) $(GSL_CFLAGS) $(MUPARSER_CFLAGS) -Wno-deprecated
EXTRA_DIST = *.h *.xpm

bin_PROGRAMS = cronos
cronos_SOURCES = cronos.cc mainwin.cc xyplot.cc transform.cc otherplots.cc tsmrender.cc dialogs.cc epsdc.cc cubespline.cc otherpopups.cc textrender.cc tsinterface.cc
cronos_LDADD = ../mslib/libms.a
cronos_LDFLAGS = $(MUPARSER_LIBS) $(GTKMM_LIBS) $(GSL_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -llapack
# cronos_CXXFLAGS = -ggdb -O0
