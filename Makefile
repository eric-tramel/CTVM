INCLUDE_DIR=include
SRC_DIR=src
BIN_DIR=bin
TEST_DIR=test
LIB_DIR=lib
MK_BOOST_LIB=$(BOOST_LIBDIR)
MK_BOOST_INC=$(BOOST_INCDIR)

MAGICK_CFLAG=`Magick++-config --cppflags --cxxflags`
MAGICK_LDFLAG=`Magick++-config --ldflags --libs`

CXX=g++
CPPFLAGS=-I$(INCLUDE_DIR) -I$(MK_BOOST_INC) $(MAGICK_CFLAG)
LDLIBS=-lboost_system -lboost_random -lboost_date_time
LDFLAGS=-L$(MK_BOOST_LIB) $(MAGICK_LDFLAG) $(LDLIBS)
DEPS=$(INCLUDE_DIR)/ctvm.h $(INCLUDE_DIR)/ctvm_util.h

all: checkdir ctvmlib executable test1


windows:
	$(CXX) -shared $(CPPFLAGS) -o $(LIB_DIR)/libctvm.dll $(SRC_DIR)/ctvm.cpp $(LDFLAGS)
	$(CXX) -shared $(CPPFLAGS) -o $(LIB_DIR)/libctvm_util.dll $(SRC_DIR)/ctvm_util.cpp $(LDFLAGS)
	$(CXX) -Wall $(CPPFLAGS) -Llib $(LDFLAGS) -lctvm -lctvm_util -c -o $(TEST_DIR)/test1.o $(TEST_DIR)/test1.cpp
	$(CXX) -Wall -Llib $(LDFLAGS) -lctvm -lctvm_util -o $(BIN_DIR)/test1 $(TEST_DIR)/test1.o


# ctvmlib: $(DEPS)
# 	$(CXX) -dynamiclib -fPIC  $(CPPFLAGS) -o $(LIB_DIR)/libctvm.dylib $(SRC_DIR)/ctvm.cpp $(LDFLAGS)
# 	$(CXX) -dynamiclib -fPIC  $(CPPFLAGS) -o $(LIB_DIR)/libctvm_util.dylib $(SRC_DIR)/ctvm_util.cpp $(LDFLAGS)

ctvmlib: $(DEPS)
	# Compile both of the libraries to object files
	$(CXX) -Wall -fPIC $(CPPFLAGS) -o $(SRC_DIR)/ctvm.o -c $(SRC_DIR)/ctvm.cpp
	$(CXX) -Wall -fPIC $(CPPFLAGS) -o $(SRC_DIR)/ctvm_util.o -c $(SRC_DIR)/ctvm_util.cpp
	# Link object files together into shared libraries
	$(CXX) -shared $(LDFLAGS) -o $(LIB_DIR)/libctvm.dylib $(SRC_DIR)/ctvm.o
	$(CXX) -shared $(LDFLAGS) -o $(LIB_DIR)/libctvm_util.dylib $(SRC_DIR)/ctvm_util.o


$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp $(LIB_DIR)/libctvm.dylib $(LIB_DIR)/libctvm_util.dylib
	$(CXX) $(CPPFLAGS) -c -o $@ $< 

test1: $(TEST_DIR)/test1.o
	$(CXX) -Llib $(LDFLAGS) -lctvm -lctvm_util -o $(BIN_DIR)/test1 $(TEST_DIR)/test1.o 

executable: ctvmlib
	$(CXX) $(CPPFLAGS) -Llib $(LDFLAGS) -lctvm -lctvm_util -o $(BIN_DIR)/ctvm-recover $(SRC_DIR)/ctvm_recover.cpp

clean:
	rm -f $(BIN_DIR)/*
	rm -f $(SRC_DIR)/*.o
	rm -f $(TEST_DIR)/*.o
	rm -f $(LIB_DIR)/*

test: clean all
	$(BIN_DIR)/test1
	$(BIN_DIR)/ctvm-recover test/data/testSino.png test/data/testAngles.dat a3


checkdir: 
	if [ -d "$(LIB_DIR)" ]; then \
		echo " "; \
	else \
		mkdir $(LIB_DIR); \
	fi

	if [ -d "$(BIN_DIR)" ]; then \
		echo " "; \
	else \
		mkdir $(BIN_DIR); \
	fi