# Dan Shervheim
# danielshervheim.com

CXX = g++
CXXFLAGS = -std=c++11  # c++ 11 required for initializing non-static class members.
BIN = atmosphere
SRC = src
BLD = build
OPENEXR = -l IlmImf -l half -I /usr/local/include/OpenEXR/

# default target
all: $(BIN)

# make a folder to store executable.
$(BLD):
	mkdir -p $(BLD)

# cleanup the compiled executable and build directory.
clean:
	rm -rf $(BLD)

# compile and link the source code from the source directory
# into an executable, and store it in the build directory.
$(BIN): $(BLD) $(SRC)/*.cc
	$(CXX) $(CXXFLAGS) -o $(BLD)/$(BIN) $(OPENEXR) $(SRC)/*.cc
