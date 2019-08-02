
CC=xcrun --sdk iphoneos clang -arch arm64 --sysroot=/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS12.4.sdk \
-I/Users/goodclass/Desktop/ftoc_temp \
-I/Users/goodclass/Downloads/PythonSite/PythonApple/build/macOS/python/lib/python3.7/site-packages/numpy/core/include \
-I/Users/goodclass/Downloads/PythonSite/PythonApple/build/macOS/python/include/python3.7m
CXX=xcrun --sdk iphoneos clang++ -arch arm64 --sysroot=/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS12.4.sdk 
CFLAGS=-mios-version-min=11.0

SRCDIR := .

CSOURCES := $(shell find $(SRCDIR) -name "*.c")
CXXSOURCES := $(shell find $(SRCDIR) -name "*.cxx")
OBJECTS  := $(patsubst %.c, %.o, $(CSOURCES))
CXXOBJECTS  := $(patsubst %.cxx, %.o, $(CXXSOURCES))

.PHONY: all

all: $(OBJECTS) $(CXXOBJECTS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cxx
	$(CC) $(CXXFLAGS) -c $< -o $@

