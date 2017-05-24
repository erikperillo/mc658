# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


SET(CPACK_ALL_INSTALL_TYPES "Full;Developer")
SET(CPACK_BINARY_7Z "")
SET(CPACK_BINARY_BUNDLE "")
SET(CPACK_BINARY_CYGWIN "")
SET(CPACK_BINARY_DEB "")
SET(CPACK_BINARY_DRAGNDROP "")
SET(CPACK_BINARY_IFW "")
SET(CPACK_BINARY_NSIS "")
SET(CPACK_BINARY_OSXX11 "")
SET(CPACK_BINARY_PACKAGEMAKER "")
SET(CPACK_BINARY_PRODUCTBUILD "")
SET(CPACK_BINARY_RPM "")
SET(CPACK_BINARY_STGZ "")
SET(CPACK_BINARY_TBZ2 "")
SET(CPACK_BINARY_TGZ "")
SET(CPACK_BINARY_TXZ "")
SET(CPACK_BINARY_TZ "")
SET(CPACK_BINARY_WIX "")
SET(CPACK_BINARY_ZIP "")
SET(CPACK_BUILD_SOURCE_DIRS "/home/erik/downloads/rand/lemon-1.3.1;/home/erik/downloads/rand/lemon-1.3.1/build")
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_COMPONENTS_ALL "headers;library;html_documentation;bin")
SET(CPACK_COMPONENTS_ALL_SET_BY_USER "TRUE")
SET(CPACK_COMPONENT_BIN_DESCRIPTION "Command line utilities")
SET(CPACK_COMPONENT_BIN_DISPLAY_NAME "Command line utilities")
SET(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "Components needed to develop software using LEMON")
SET(CPACK_COMPONENT_GROUP_DOCUMENTATION_DESCRIPTION "Documentation of LEMON")
SET(CPACK_COMPONENT_HEADERS_DEPENDS "library")
SET(CPACK_COMPONENT_HEADERS_DESCRIPTION "C++ header files")
SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ headers")
SET(CPACK_COMPONENT_HEADERS_GROUP "Development")
SET(CPACK_COMPONENT_HEADERS_INSTALL_TYPES "Developer;Full")
SET(CPACK_COMPONENT_HTML_DOCUMENTATION_DESCRIPTION "Doxygen generated documentation")
SET(CPACK_COMPONENT_HTML_DOCUMENTATION_DISPLAY_NAME "HTML documentation")
SET(CPACK_COMPONENT_HTML_DOCUMENTATION_GROUP "Documentation")
SET(CPACK_COMPONENT_HTML_DOCUMENTATION_INSTALL_TYPES "Full")
SET(CPACK_COMPONENT_LIBRARY_DESCRIPTION "DLL and import library")
SET(CPACK_COMPONENT_LIBRARY_DISPLAY_NAME "Dynamic-link library")
SET(CPACK_COMPONENT_LIBRARY_GROUP "Development")
SET(CPACK_COMPONENT_LIBRARY_INSTALL_TYPES "Developer;Full")
SET(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
SET(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
SET(CPACK_GENERATOR "NSIS")
SET(CPACK_INSTALL_CMAKE_PROJECTS "/home/erik/downloads/rand/lemon-1.3.1/build;LEMON;ALL;/")
SET(CPACK_INSTALL_PREFIX "/usr/local")
SET(CPACK_MODULE_PATH "/home/erik/downloads/rand/lemon-1.3.1/cmake")
SET(CPACK_NSIS_CONTACT "lemon-user@lemon.cs.elte.hu")
SET(CPACK_NSIS_CREATE_ICONS_EXTRA "
    CreateShortCut \"$SMPROGRAMS\\$STARTMENU_FOLDER\\Documentation.lnk\" \"$INSTDIR\\share\\doc\\index.html\"
    ")
SET(CPACK_NSIS_DELETE_ICONS_EXTRA "
    !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    Delete \"$SMPROGRAMS\\$MUI_TEMP\\Documentation.lnk\"
    ")
SET(CPACK_NSIS_DISPLAY_NAME "LEMON 1.3.1 LEMON")
SET(CPACK_NSIS_DISPLAY_NAME_SET "TRUE")
SET(CPACK_NSIS_HELP_LINK "http:\\\\lemon.cs.elte.hu")
SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\lemon.ico")
SET(CPACK_NSIS_INSTALLER_ICON_CODE "")
SET(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
SET(CPACK_NSIS_MUI_ICON "/home/erik/downloads/rand/lemon-1.3.1/cmake/nsis/lemon.ico")
SET(CPACK_NSIS_MUI_UNIICON "/home/erik/downloads/rand/lemon-1.3.1/cmake/nsis/uninstall.ico")
SET(CPACK_NSIS_PACKAGE_NAME "LEMON 1.3.1 LEMON")
SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\lemon.cs.elte.hu")
SET(CPACK_OUTPUT_CONFIG_FILE "/home/erik/downloads/rand/lemon-1.3.1/build/CPackConfig.cmake")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "/")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake-3.8/Templates/CPack.GenericDescription.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LEMON - Library for Efficient Modeling and Optimization in Networks")
SET(CPACK_PACKAGE_FILE_NAME "LEMON-1.3.1-Linux")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "LEMON 1.3.1")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "LEMON 1.3.1")
SET(CPACK_PACKAGE_NAME "LEMON")
SET(CPACK_PACKAGE_RELOCATABLE "true")
SET(CPACK_PACKAGE_VENDOR "EGRES")
SET(CPACK_PACKAGE_VERSION "1.3.1")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_RESOURCE_FILE_LICENSE "/home/erik/downloads/rand/lemon-1.3.1/LICENSE")
SET(CPACK_RESOURCE_FILE_README "/usr/share/cmake-3.8/Templates/CPack.GenericDescription.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "/usr/share/cmake-3.8/Templates/CPack.GenericWelcome.txt")
SET(CPACK_SET_DESTDIR "OFF")
SET(CPACK_SOURCE_7Z "")
SET(CPACK_SOURCE_CYGWIN "")
SET(CPACK_SOURCE_GENERATOR "TBZ2;TGZ;TXZ;TZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/erik/downloads/rand/lemon-1.3.1/build/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_RPM "OFF")
SET(CPACK_SOURCE_TBZ2 "ON")
SET(CPACK_SOURCE_TGZ "ON")
SET(CPACK_SOURCE_TXZ "ON")
SET(CPACK_SOURCE_TZ "ON")
SET(CPACK_SOURCE_ZIP "OFF")
SET(CPACK_SYSTEM_NAME "Linux")
SET(CPACK_TOPLEVEL_TAG "Linux")
SET(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/erik/downloads/rand/lemon-1.3.1/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()