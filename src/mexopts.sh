#
# mexopts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Note: For the version of system compiler supported with this release,
#       refer to the Supported and Compatible Compiler List at:
#       http://www.mathworks.com/support/compilers/current_release/
#
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building MEX-files
#
# Copyright 1984-2011 The MathWorks, Inc.
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    if [ "$ENTRYPOINT" = "mexLibrary" ]; then
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lmwservices -lut"
    else  
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat"
    fi
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        glnx86)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh glnx86 12
#----------------------------------------------------------------------------
            ;;
        glnxa64)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            RPATH="$RPATH -Wl,-rpath,./"
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            # CkeyLinkerName: GNU ld
            # CkeyLinkerVersion:
            CC='gcc-4.7'
            CFLAGS='-I../include -I$MKLROOT/include -D_GNU_SOURCE '
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS -fPIC -fno-omit-frame-pointer -m64"
            MKLLIBS="-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group"
            CLIBS="-L./ -leigpow $MKLLIBS"
            CLIBS="$CLIBS $RPATH $MLIBS -lgomp -lpthread -lm -ldl"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion:
            # C++keyLinkerName: GNU ld
            # C++keyLinkerVersion:  
            CXX='g++-4.7'
            CXXFLAGS='-I$MKLROOT/include -ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -fopenmp"
            CXXLIBS="-L$MKLROOT/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread"
            CXXLIBS="$CXXLIBS $RPATH $MLIBS -lm -ldl"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: gfortran
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion:
            # FortrankeyLinkerName: GNU ld
            # FortrankeyLinkerVersion:  
#
            FC='gfortran'
            FFLAGS='-fexceptions -fbackslash'
            FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexa64'
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol64)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh sol64 12
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh mac 12
#----------------------------------------------------------------------------
            ;;
        maci)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh maci 12
#----------------------------------------------------------------------------
            ;;
        maci64)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: Clang
            # CkeyManufacturer: Apple
            # CkeyLanguage: C
            # CkeyVersion:
            # CkeyLinkerName:
            # CkeyLinkerVersion:
            CC='xcrun  -sdk macosx10.8  clang'
## workaround clang defect temporarily use line below           SDKROOT='/Developer/SDKs/MacOSX10.6.sdk'
# compute SDK root on the fly
# target 10.8
            MW_SDKROOT_TMP="find `xcode-select -print-path` -name MacOSX10.8.sdk"
			MW_SDKROOT=`$MW_SDKROOT_TMP`
            MACOSX_DEPLOYMENT_TARGET='10.7'
            ARCHS='x86_64'
            CFLAGS="-fno-common -arch $ARCHS -isysroot $MW_SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CFLAGS="$CFLAGS  -fexceptions"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lc++"
            # C++keyName: Clang++
            # C++keyManufacturer: Apple
            # C++keyLanguage: C++
            # C++keyVersion:
            # C++keyLinkerName:
            # C++keyLinkerVersion:
            CXX='xcrun  -sdk macosx10.8  clang++'
            CXXFLAGS="-fno-common -fexceptions -arch $ARCHS -isysroot $MW_SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET -std=c++11 -stdlib=libc++"
	    CXXLIBS="-lc++"
            CXXLIBS="$CXXLIBS $MLIBS"
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: Intel Fortran
            # FortrankeyManufacturer: Intel
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            # FortrankeyLinkerName: 
            # FortrankeyLinkerVersion:
            FC='ifort'
            FFLAGS='-fexceptions -fpp -mp1 -fp-model source -assume bscc -D__LP64__'
            FC_LIBDIR="$IFORT_COMPILER14/compiler/lib"
            FLIBS="$MLIBS -L$FC_LIBDIR -lifcore -limf -lintlc -lirc -lsvml"
            FOPTIMFLAGS='-O2'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmaci64'
            LDFLAGS="-arch $ARCHS -Wl,-syslibroot,$MW_SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            LDFLAGS="$LDFLAGS -bundle -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
