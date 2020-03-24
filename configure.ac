AC_INIT([mapseq],[1.2.6])

: ${CXXFLAGS=-O3}

AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
AC_LANG([C++])
AC_PROG_CXX
AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_INSTALL

PKG_PROG_PKG_CONFIG
LT_INIT


AX_PTHREAD()

LIBS="$PTHREAD_LIBS $LIBS"
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"

AX_CXX_COMPILE_STDCXX_11(,mandatory)
#AC_CANONICAL_SYSTEM

case "${host_os}" in
  darwin*)
    echo ===========================================================
    echo Setting up build environment for ${target_cpu}${target_os}
    echo ===========================================================
    host_os=osx
    LDLIBS="${LDLIBS} -framework Carbon"
    LDFLAGS="${LDFLAGS} -headerpad_max_install_names"
    ;;
  *linux*)
    echo ===========================================================
    echo Setting up build environment for ${target_cpu}${target_os}
    echo ===========================================================
    host_os=xlib
    ;;
  *)
    echo ===========================================================
    echo Setting up build environment for ${target_cpu}${target_os}
    echo ===========================================================
    AC_MSG_ERROR([ "Unsupported target operating system, please correct" ])
  esac

AC_ARG_WITH(syseutils,[AS_HELP_STRING([--with-syseutils],[Use system installed eutils.])],,with_syseutils=no)
if test "x$with_syseutils" = "xno"; then
  AC_CONFIG_SUBDIRS([libs/eutils])
else
  AM_PATH_EUTILS
  if test "x$no_eutils" = "xyes"; then
    AC_MSG_ERROR([
*** eutils library not found.
*** You can use the eutils version included in this package by
*** not using the --with-syseutils argument.])
  fi
  LIBS="${LIBS} ${EUTILS_LIBS}"
fi
AM_CONDITIONAL(USE_SYSEUTILS, test ! "x$with_syseutils" = "xno")


#AM_INIT_AUTOMAKE("mapseq","1.0.0")

#CPPFLAGS="${CPPFLAGS} -I/usr/ilog/cplex101/include"

# eutils, edlib
#AM_PATH_EUTILS(0.6.0)
#if test "$no_eutils" != yes; then
#  CXXFLAGS="$CXXFLAGS $EUTILS_CFLAGS"
#  CFLAGS="$CFLAGS $EUTILS_CFLAGS"
#  LIBS="$LIBS $EUTILS_LIBS"
#  AC_DEFINE_UNQUOTED(HAVE_EUTILS,1,[define when libeutils exists])
#else
# AC_MSG_ERROR([
#   Please install the EUTILS library.
#   If EUTILS is installed you might need to indicate where using the argument: --with-eutils-prefix=DIR
# ])
#fi

AC_ARG_ENABLE([makestatic], [  --enable-makestatic    compile statically linked executable. ],[enable_makestatic=$enableval],[enable_makestatic=no])
if test "x$enable_makestatic" = "xyes"; then
  AC_MSG_NOTICE([Enabling static compilation])
  LDFLAGS="${LDFLAGS} -all-static"
fi



AC_CONFIG_HEADERS(config.h)
AX_PREFIX_CONFIG_H

AC_CONFIG_FILES([Makefile])
#AC_CONFIG_SUBDIRS([eutils])

AC_OUTPUT
