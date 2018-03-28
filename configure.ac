dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[SCIPhI],
	[0.1],
	[jochen.singer@bsse.ethz.ch],
	[sciphi],
	[https://gitlab.com/singerj/SCITE])
dnl AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src/findBestTrees.cpp])
AC_LANG([C++])



dnl ========================================================================
dnl Check whether we want to set defaults for CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ========================================================================
AC_MSG_CHECKING([whether configure should try to set CXXFLAGS/CPPFLAGS/LDFLAGS])
AS_IF([test "x${CXXFLAGS+set}" = "xset" || test "x${CPPFLAGS+set}" = "xset" || test "x${LDFLAGS+set}" = "xset"],
	[enable_flags_setting=no],
	[enable_flags_setting=yes]
)
AC_MSG_RESULT([${enable_flags_setting}])

AX_CHECK_ENABLE_DEBUG
AS_IF([test "x${enable_debug}" = "xno"],
	[
		cxxflags_test="-O3"
		cppflags_test=""

		AC_CANONICAL_HOST
		AS_CASE([${host_os}],
			[darwin*],
			[ldflags_test="-Wl,-dead_strip_dylibs"],
			[linux*],
			[ldflags_test="-Wl,-O1 -Wl,--as-needed"]
		)
	]
)



dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
AX_COMPILER_VENDOR
AX_COMPILER_VERSION
AC_PROG_SED
AC_PROG_CXX
AC_CHECK_HEADERS([cstdint])
AX_PTHREAD



dnl ==================================
dnl Set CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ==================================
AS_IF([test "x${enable_flags_setting}" = "xyes" && test "x${enable_debug}" = "xno"],
	[
	AX_APPEND_COMPILE_FLAGS([${cxxflags_test} -Wall -Wno-redeclared-class-member], [CXXFLAGS])
	CXXFLAGS=$( echo ${CXXFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )
	
	AX_APPEND_COMPILE_FLAGS([${cppflags_test}], [CPPFLAGS])
	CPPFLAGS=$( echo ${CPPFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )
	
	AX_APPEND_LINK_FLAGS([${ldflags_test}], [LDFLAGS])
	LDFLAGS=$( echo ${LDFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )
	]
)



dnl ===============
dnl Check for C++14
dnl ===============
AX_CXX_COMPILE_STDCXX([14], [noext], [mandatory])


dnl ===================
dnl Initialise Automake
dnl ===================
AM_INIT_AUTOMAKE([1.15 foreign dist-bzip2 no-dist-gzip subdir-objects silent-rules])



dnl ===============
dnl Check for Boost
dnl ===============
AX_BOOST_BASE([1.50], [], AC_MSG_ERROR([could not find a suitable boost installation.]))
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_SYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_FILESYSTEM
AX_BOOST_THREAD
AX_BOOST_REGEX


dnl Test whether the unit test framework needs -DBOOST_TEST_DYN_LINK
AC_MSG_CHECKING([whether the Boost::Unit_Test_Framework library needs -DBOOST_TEST_DYN_LINK])
save_CPPFLAGS="$CPPFLAGS"
save_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
BOOST_UNIT_TEST_FRAMEWORK_CPPFLAGS=""
LIBS="$BOOST_LDFLAGS $BOOST_UNIT_TEST_FRAMEWORK_LIB"
AC_LINK_IFELSE([
    AC_LANG_SOURCE(
        [[
        #define BOOST_TEST_MODULE MyTest
        #include <boost/test/unit_test.hpp>
        BOOST_AUTO_TEST_CASE( my_test )
        {
            BOOST_CHECK(1 == 1);
        }
        ]]
    )],
    [AC_MSG_RESULT([no])],
    [BOOST_UNIT_TEST_FRAMEWORK_CPPFLAGS="-DBOOST_TEST_DYN_LINK"
    CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS $BOOST_UNIT_TEST_FRAMEWORK_CPPFLAGS"
    AC_LINK_IFELSE([
        AC_LANG_SOURCE(
            [[
            #define BOOST_TEST_MODULE MyTest
            #include <boost/test/unit_test.hpp>
            BOOST_AUTO_TEST_CASE( my_test )
            {
                BOOST_CHECK(1 == 1);
            }
            ]]
        )],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_ERROR([could not determine flags])]
    )]
)
CPPFLAGS="$save_CPPFLAGS"
LIBS="$save_LIBS"
AC_SUBST([BOOST_UNIT_TEST_FRAMEWORK_CPPFLAGS])



dnl ===============
dnl Check for Seqan
dnl ===============
AC_ARG_VAR([SEQAN_INCLUDEDIR], [contains the include path to the Seqan headers])
AC_MSG_CHECKING([whether seqan headers are available])
save_cppflags="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS -I$SEQAN_INCLUDEDIR"
AC_LINK_IFELSE([
    AC_LANG_SOURCE(
        [[
        #include <seqan/basic.h>
        int main() { seqan::Pair<int, int> pair; return 0; }
        ]]
    )],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_ERROR([could not find seqan headers! If the seqan headers are not located within the compiler's default search path, specify them with the SEQAN_INCLUDEDIR variable to the configure script.])]
)
CPPFLAGS="$save_cppflags"
PKG_CHECK_MODULES([ZLIB], [zlib])

dnl ==============
dnl Check for DLIB
dnl ==============
AC_ARG_VAR([DLIB_INCLUDEDIR], [contains the include path to the Seqan headers])
AC_MSG_CHECKING([whether dlib headers are available])
save_cppflags="$CPPFLAGS"
DLIB_CPPFLAGS="-DDLIB_NO_GUI_SUPPORT -lpthread "
AC_SUBST([DLIB_CPPFLAGS])
DLIB_LDFLAGS="-DDLIB_NO_GUI_SUPPORT $DLIB_INCLUDEDIR/dlib/all/source.cpp"
AC_SUBST([DLIB_LDFLAGS])
CPPFLAGS="-std=c++11 -I$DLIB_INCLUDEDIR $DLIB_CPPFLAGS $DLIB_LDFLAGS"
AC_LINK_IFELSE([
    AC_LANG_SOURCE(
        [[
        #include <dlib/optimization.h>
        #include <dlib/global_optimization.h>
        int main() { dlib::matrix<double,0,1> column_vector;; return 0; }
        ]]
    )],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_ERROR([could not find dlib headers! If the dlib headers are not located within the compiler's default search path, specify them with the DLIB_INCLUDEDIR variable to the configure script.])]
)
CPPFLAGS="$save_cppflags"

dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([Makefile])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	$PACKAGE_NAME $VERSION
	Using $ax_cv_cxx_compiler_vendor, $ax_cv_cxx_compiler_version

	CXX:        $CXX
	CXXFLAGS:   $CXXFLAGS
	CPPFLAGS:   $CPPFLAGS
	LDFLAGS:    $LDFLAGS
])
