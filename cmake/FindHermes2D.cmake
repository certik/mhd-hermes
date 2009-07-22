execute_process(
	COMMAND python -c "import hermes2d; print hermes2d.get_include()"
	OUTPUT_VARIABLE HERMES2D_INCLUDE
	)
string(STRIP ${HERMES2D_INCLUDE} HERMES2D_INCLUDE)

execute_process(
	COMMAND python -c "import hermes2d; print hermes2d.get_pxd_include()"
	OUTPUT_VARIABLE HERMES2D_PXD_INCLUDE
	)
string(STRIP ${HERMES2D_PXD_INCLUDE} HERMES2D_PXD_INCLUDE)
FIND_PATH(HERMES2D_PXD_INCLUDE hermes2d.pxd ${HERMES2D_PXD_INCLUDE})

execute_process(
	COMMAND python -c "import hermes2d; print hermes2d.get_lib()"
	OUTPUT_VARIABLE HERMES2D_LIB
	)
string(STRIP ${HERMES2D_LIB} HERMES2D_LIB)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HERMES2D DEFAULT_MSG HERMES2D_PXD_INCLUDE)
