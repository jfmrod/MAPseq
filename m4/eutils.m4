# Configure path for the EasyC++ library
# Joao F. Matias Rodrigues <nyone@enyon.homeip.net>, September 2007


AC_DEFUN([AM_PATH_EUTILS],
[
AC_ARG_WITH(eutils-prefix,[  --with-eutils-prefix=PFX   Prefix where EUTILS is installed (optional)],
            eutils_prefix="$withval", eutils_prefix="")
AC_ARG_WITH(eutils-exec-prefix,[  --with-eutils-exec-prefix=PFX Exec prefix where EUTILS is installed (optional)],
            eutils_exec_prefix="$withval", eutils_exec_prefix="")
AC_ARG_ENABLE(eutilstest, [  --disable-eutilstest       Do not try to compile and run a test EUTILS program],
		    , enable_eutilstest=yes)

  if test "x${EUTILS_CONFIG+set}" != xset ; then
     if test "x$eutils_prefix" != x ; then
         EUTILS_CONFIG="$eutils_prefix/bin/eutils-config"
     fi
     if test "x$eutils_exec_prefix" != x ; then
        EUTILS_CONFIG="$eutils_exec_prefix/bin/eutils-config"
     fi
  fi

  AC_PATH_PROG(EUTILS_CONFIG, eutils-config, no)
  min_eutils_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for EUTILS - version >= $min_eutils_version)
  no_eutils=""
  if test "$EUTILS_CONFIG" = "no" ; then
    no_eutils=yes
  else
    EUTILS_CXXFLAGS=`$EUTILS_CONFIG --cxxflags`
    EUTILS_LIBS=`$EUTILS_CONFIG --libs`

    eutils_major_version=`$EUTILS_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${eutils_major_version}" = "x" ; then
       eutils_major_version=0
    fi

    eutils_minor_version=`$EUTILS_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${eutils_minor_version}" = "x" ; then
       eutils_minor_version=0
    fi

    eutils_micro_version=`$EUTILS_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${eutils_micro_version}" = "x" ; then
       eutils_micro_version=0
    fi

    if test "x$enable_eutilstest" = "xyes" ; then
      ac_save_CXXFLAGS="$CXXFLAGS"
      ac_save_LIBS="$LIBS"
      CXXFLAGS="$CXXFLAGS $EUTILS_CXXFLAGS"
      LIBS="$LIBS $EUTILS_LIBS"

      rm -f conf.eutilstest
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_RUN([
#include <eutils/estr.h>

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

  system ("touch conf.eutilstest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_eutils_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_eutils_version");
     exit(1);
   }

   if (($eutils_major_version > major) ||
      (($eutils_major_version == major) && ($eutils_minor_version > minor)) ||
      (($eutils_major_version == major) && ($eutils_minor_version == minor) && ($eutils_micro_version >= micro)))
    {
      exit(0);
    }
  else
    {
      printf("\n*** 'eutils-config --version' returned %d.%d.%d, but the minimum version\n", $eutils_major_version, $eutils_minor_version, $eutils_micro_version);
      printf("*** of EUTILS required is %d.%d.%d. If eutils-config is correct, then it is\n", major, minor, micro);
      printf("*** best to upgrade to the required version.\n");
      printf("*** If eutils-config was wrong, set the environment variable EUTILS_CONFIG\n");
      printf("*** to point to the correct copy of eutils-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
      exit(1);
    }
}

],, no_eutils=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       AC_LANG_RESTORE
       CXXFLAGS="$ac_save_CXXFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_eutils" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$EUTILS_CONFIG" = "no" ; then
       echo "*** The eutils-config script installed by EUTILS could not be found"
       echo "*** If EUTILS was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the EUTILS_CONFIG environment variable to the"
       echo "*** full path to eutils-config."
     else
       if test -f conf.eutilstest ; then
        :
       else
          echo "*** Could not run EUTILS test program, checking why..."
          CXXFLAGS="$CXXFLAGS $EUTILS_CXXFLAGS"
          LIBS="$LIBS $EUTILS_LIBS"
          AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_TRY_LINK([
#include <stdio.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding EUTILS or finding the wrong"
          echo "*** version of EUTILS. If it is not finding EUTILS, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means EUTILS was incorrectly installed"
          echo "*** or that you have moved EUTILS since it was installed. In the latter case, you"
          echo "*** may want to edit the eutils-config script: $EUTILS_CONFIG" ])
          AC_LANG_RESTORE
          CXXFLAGS="$ac_save_CXXFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(EUTILS_CXXFLAGS)
  AC_SUBST(EUTILS_LIBS)
  rm -f conf.eutilstest
])


