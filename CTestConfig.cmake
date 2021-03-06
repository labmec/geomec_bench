## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "PZ-External_Projects")
set(CTEST_NIGHTLY_START_TIME "00:00:00 BRT")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "www.labmec.org.br")
set(CTEST_DROP_LOCATION "/pz/CDash/submit.php?project=PZ-External_Projects")
set(CTEST_DROP_SITE_CDASH TRUE)
