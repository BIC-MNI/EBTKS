LINK_LIBRARIES(EBTKS)

ADD_DEFINITIONS(-DUSE_COMPMAT)
ADD_DEFINITIONS(-DUSE_DBLMAT)
ADD_DEFINITIONS(-DUSE_FCOMPMAT)

ADD_EXECUTABLE(ebtks_test_fft ebtks_test_fft.cc)

#TODO: add more test cases!


ADD_TEST(ebtks_test_fft ebtks_test_fft)
