#----------------------------------------------------------------------
#
# File:     CMakeList.txt
# Project:  libsscv
# Author:   Jeffrey Byrne
#
# $Id: CMakeLists.txt 51 2009-01-29 18:46:14Z jbyrne $
#
#----------------------------------------------------------------------
# MEX
IF (WIN32)
   SET(MEX_EXTENSION mexw32)
ENDIF (WIN32)
IF (UNIX)
   SET(MEX_EXTENSION mexglx)
ENDIF (UNIX)
set(MEX_ARGS -O)

# corrDn
IF (MEX_COMPILER)
  SET(MEX_BASE corrDn)
  SET(MEX_SRC ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/${MEX_BASE}.c) 
  SET(MEX_SRC ${MEX_SRC} ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/wrap.c) 
  SET(MEX_SRC ${MEX_SRC} ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/convolve.c) 
  SET(MEX_SRC ${MEX_SRC} ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/edges.c) 
  SET(MEX_OUT ${LIBSSCV_BINARY_DIR}/${MEX_BASE}.${MEX_EXTENSION}) 
  ADD_CUSTOM_TARGET(pyr-${MEX_BASE} ALL DEPENDS ${MEX_OUT})
  ADD_CUSTOM_COMMAND(
    OUTPUT    ${MEX_OUT}
    DEPENDS   ${MEX_SRC}
    COMMAND   ${MEX_COMPILER} 
    ARGS      ${MEX_ARGS}
	            ${MEX_SRC}
              -output ${MEX_OUT}
    COMMENT   "matlabPyrTools-1.3 - ${MEX_BASE} MEX"
  )
ENDIF(MEX_COMPILER)

# upConv
IF (MEX_COMPILER)
  SET(MEX_BASE upConv)
  SET(MEX_SRC ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/${MEX_BASE}.c) 
  SET(MEX_SRC ${MEX_SRC} ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/wrap.c) 
  SET(MEX_SRC ${MEX_SRC} ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/convolve.c) 
  SET(MEX_SRC ${MEX_SRC} ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/edges.c) 
  SET(MEX_OUT ${LIBSSCV_BINARY_DIR}/${MEX_BASE}.${MEX_EXTENSION}) 
  ADD_CUSTOM_TARGET(pyr-${MEX_BASE} ALL DEPENDS ${MEX_OUT})
  ADD_CUSTOM_COMMAND(
    OUTPUT    ${MEX_OUT}
    DEPENDS   ${MEX_SRC}
    COMMAND   ${MEX_COMPILER} 
    ARGS      ${MEX_ARGS}
	            ${MEX_SRC}
              -output ${MEX_OUT}
    COMMENT   "matlabPyrTools-1.3 - ${MEX_BASE} MEX"
  )
ENDIF(MEX_COMPILER)

# pointOp
IF (MEX_COMPILER)
  SET(MEX_BASE pointOp)
  SET(MEX_SRC ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/${MEX_BASE}.c) 
  SET(MEX_OUT ${LIBSSCV_BINARY_DIR}/${MEX_BASE}.${MEX_EXTENSION}) 
  ADD_CUSTOM_TARGET(pyr-${MEX_BASE} ALL DEPENDS ${MEX_OUT})
  ADD_CUSTOM_COMMAND(
    OUTPUT    ${MEX_OUT}
    DEPENDS   ${MEX_SRC}
    COMMAND   ${MEX_COMPILER} 
    ARGS      ${MEX_ARGS}
	            ${MEX_SRC}
              -output ${MEX_OUT}
    COMMENT   "matlabPyrTools-1.3 - ${MEX_BASE} MEX"
  )
ENDIF(MEX_COMPILER)

# histo
IF (MEX_COMPILER)
  SET(MEX_BASE histo)
  SET(MEX_SRC ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/${MEX_BASE}.c) 
  SET(MEX_OUT ${LIBSSCV_BINARY_DIR}/${MEX_BASE}.${MEX_EXTENSION}) 
  ADD_CUSTOM_TARGET(pyr-${MEX_BASE} ALL DEPENDS ${MEX_OUT})
  ADD_CUSTOM_COMMAND(
    OUTPUT    ${MEX_OUT}
    DEPENDS   ${MEX_SRC}
    COMMAND   ${MEX_COMPILER} 
    ARGS      ${MEX_ARGS}
	            ${MEX_SRC}
              -output ${MEX_OUT}
    COMMENT   "matlabPyrTools-1.3 - ${MEX_BASE} MEX"
  )
ENDIF(MEX_COMPILER)

# range2
IF (MEX_COMPILER)
  SET(MEX_BASE range2)
  SET(MEX_SRC ${LIBSSCV_SOURCE_DIR}/libs/matlabPyrTools-1.3/MEX/${MEX_BASE}.c) 
  SET(MEX_OUT ${LIBSSCV_BINARY_DIR}/${MEX_BASE}.${MEX_EXTENSION}) 
  ADD_CUSTOM_TARGET(pyr-${MEX_BASE} ALL DEPENDS ${MEX_OUT})
  ADD_CUSTOM_COMMAND(
    OUTPUT    ${MEX_OUT}
    DEPENDS   ${MEX_SRC}
    COMMAND   ${MEX_COMPILER} 
    ARGS      ${MEX_ARGS}
	            ${MEX_SRC}
              -output ${MEX_OUT}
    COMMENT   "matlabPyrTools-1.3 - ${MEX_BASE} MEX"
  )
ENDIF(MEX_COMPILER)

