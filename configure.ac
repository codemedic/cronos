AC_INIT([cronos],[1.40],[Anthony Brockwell <abrock@stat.cmu.edu>])
AM_INIT_AUTOMAKE


#########################################################################
#  Dependency checks
#########################################################################
gtkmm_min_glibmm_version=2.4.0
gtkmm_min_gtk_version=2.4.0


GLIBMM_PREFIX=`pkg-config --variable=prefix glibmm-2.4`
GMMPROC_DIR=$GLIBMM_PREFIX/lib/glibmm-2.4/proc
AC_SUBST(GMMPROC_DIR)
GMMPROC=$GMMPROC_DIR/gmmproc
AC_SUBST(GMMPROC)


PKG_CHECK_MODULES(ATKMM, glibmm-2.4 >= ${gtkmm_min_glibmm_version} atk >= 1.6.0)
AC_SUBST(ATKMM_CFLAGS)
AC_SUBST(ATKMM_LIBS)

PKG_CHECK_MODULES(PANGOMM, glibmm-2.4 >= ${gtkmm_min_glibmm_version} pango >= 1.4.0)
AC_SUBST(PANGOMM_CFLAGS)
AC_SUBST(PANGOMM_LIBS)

PKG_CHECK_MODULES(GDKMM, glibmm-2.4 >= ${gtkmm_min_glibmm_version} gtk+-2.0 >= ${gtkmm_min_gtk_version})
AC_SUBST(GDKMM_CFLAGS)
AC_SUBST(GDKMM_LIBS)

PKG_CHECK_MODULES(GTKMM, glibmm-2.4 >= ${gtkmm_min_glibmm_version} gtk+-2.0 >= ${gtkmm_min_gtk_version} gtkmm-2.4 >= 2.4.0 )
AC_SUBST(GTKMM_CFLAGS)
AC_SUBST(GTKMM_LIBS)

AM_PATH_GSL(1.6,,[AC_MSG_ERROR([could not find required version of GSL])])

PKG_CHECK_MODULES(MUPARSER, muparser >= 1.25)

#AC_CHECK_LIB(lapack, dsyev_,,[AC_MSG_ERROR([Unable to locate LAPACK.])])

AC_CONFIG_FILES([Makefile src/Makefile src/main/Makefile src/mslib/Makefile])
AC_PROG_RANLIB
AC_PROG_CXX
AC_OUTPUT
