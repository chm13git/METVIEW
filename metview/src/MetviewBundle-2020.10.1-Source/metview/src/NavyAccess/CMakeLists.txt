
ecbuild_add_executable( TARGET       NavyAccess
                        SOURCES      NavyAccess.cc NavyAccess.h
                        DEFINITIONS  ${METVIEW_EXTRA_DEFINITIONS}
                        INCLUDES     ${METVIEW_STANDARD_INCLUDES}
                        LIBS         ${STANDARD_METVIEW_LIBS}
                    )

metview_module_files(ETC_FILES ObjectSpec.NavyAccess
                               NavyAccessDef
                               NavyAccessRules
                     XPM_FILES NAVYACCESS.xpm)
