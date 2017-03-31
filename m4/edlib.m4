# Configure path for the EasyC++ library
# Joao F. Matias Rodrigues <nyone@enyon.homeip.net>, September 2007


AC_DEFUN([AM_PATH_EDLIB],
[
AC_ARG_WITH(edlib-prefix,[  --with-edlib-prefix=PFX   Prefix where EDLIB is installed (optional)],
            edlib_prefix="$withval", edlib_prefix="")
AC_ARG_WITH(edlib-exec-prefix,[  --with-edlib-exec-prefix=PFX Exec prefix where EDLIB is installed (optional)],
            edlib_exec_prefix="$withval", edlib_exec_prefix="")
AC_ARG_ENABLE(edlibtest, [  --disable-edlibtest       Do not try to compile and run a test EDLIB program],
		    , enable_edlibtest=yes)

  if test "x${EDLIB_CONFIG+set}" != xset ; then
     if test "x$edlib_prefix" != x ; then
         EDLIB_CONFIG="$edlib_prefix/bin/edlib-config"
     fi
     if test "x$edlib_exec_prefix" != x ; then
        EDLIB_CONFIG="$edlib_exec_prefix/bin/edlib-config"
     fi
  fi

  AC_PATH_PROG(EDLIB_CONFIG, edlib-config, no)
  min_edlib_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for EDLIB - version >= $min_edlib_version)
  no_edlib=""
  if test "$EDLIB_CONFIG" = "no" ; then
    no_edlib=yes
  else
    EDLIB_CFLAGS=`$EDLIB_CONFIG --cxxflags`
    EDLIB_LIBS=`$EDLIB_CONFIG --libs`

    edlib_major_version=`$EDLIB_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${edlib_major_version}" = "x" ; then
       edlib_major_version=0
    fi

    edlib_minor_version=`$EDLIB_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${edlib_minor_version}" = "x" ; then
       edlib_minor_version=0
    fi

    edlib_micro_version=`$EDLIB_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${edlib_micro_version}" = "x" ; then
       edlib_micro_version=0
    fi

    if test "x$enable_edlibtest" = "xyes" ; then
      ac_save_CXXFLAGS="$CXXFLAGS"
      ac_save_LIBS="$LIBS"
      CXXFLAGS="$CXXFLAGS $EDLIB_CXXFLAGS"
      LIBS="$LIBS $EDLIB_LIBS"

      rm -f conf.edlibtest
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_RUN([
#include <edlib-2/edlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* my_strdup (const char *str);

char*
my_strdup (const char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = (char *)malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main (void)
{
  int major = 0, minor = 0, micro = 0;
  int n;
  char *tmp_version;

  system ("touch conf.edlibtest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_edlib_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_edlib_version");
     exit(1);
   }

   if (($edlib_major_version > major) ||
      (($edlib_major_version == major) && ($edlib_minor_version > minor)) ||
      (($edlib_major_version == major) && ($edlib_minor_version == minor) && ($edlib_micro_version >= micro)))
    {
      exit(0);
    }
  else
    {
      printf("\n*** 'edlib-config --version' returned %d.%d.%d, but the minimum version\n", $edlib_major_version, $edlib_minor_version, $edlib_micro_version);
      printf("*** of EDLIB required is %d.%d.%d. If edlib-config is correct, then it is\n", major, minor, micro);
      printf("*** best to upgrade to the required version.\n");
      printf("*** If edlib-config was wrong, set the environment variable EDLIB_CONFIG\n");
      printf("*** to point to the correct copy of edlib-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
      exit(1);
    }
}

],, no_edlib=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       AC_LANG_RESTORE
       CXXFLAGS="$ac_save_CXXFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_edlib" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$EDLIB_CONFIG" = "no" ; then
       echo "*** The edlib-config script installed by EDLIB could not be found"
       echo "*** If EDLIB was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the EDLIB_CONFIG environment variable to the"
       echo "*** full path to edlib-config."
     else
       if test -f conf.edlibtest ; then
        :
       else
          echo "*** Could not run EDLIB test program, checking why..."
          CXXFLAGS="$CXXFLAGS $EDLIB_CXXFLAGS"
          LIBS="$LIBS $EDLIB_LIBS"
          AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_TRY_LINK([
#include <stdio.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding EDLIB or finding the wrong"
          echo "*** version of EDLIB. If it is not finding EDLIB, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means EDLIB was incorrectly installed"
          echo "*** or that you have moved EDLIB since it was installed. In the latter case, you"
          echo "*** may want to edit the edlib-config script: $EDLIB_CONFIG" ])
          AC_LANG_RESTORE
          CXXFLAGS="$ac_save_CXXFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(EDLIB_CXXFLAGS)
  AC_SUBST(EDLIB_LIBS)
  rm -f conf.edlibtest
])


