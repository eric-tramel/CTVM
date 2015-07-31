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

all: checkdir ctvmlib_static executable test1 test_cs

ctvmlib: $(DEPS)
		# Compile both of the libraries to object files
		$(CXX) -Wall $(CPPFLAGS) -o $(SRC_DIR)/ctvm.o -c $(SRC_DIR)/ctvm.cpp
		$(CXX) -Wall $(CPPFLAGS) -o $(SRC_DIR)/ctvm_util.o -c $(SRC_DIR)/ctvm_util.cpp
		# Link object files together into shared libraries
		$(CXX) -shared -fPIC $(LDFLAGS) -o $(LIB_DIR)/cygctvm_util.dll $(SRC_DIR)/ctvm_util.o
		$(CXX) -shared -fPIC $(LDFLAGS) -o $(LIB_DIR)/cygctvm.dll $(SRC_DIR)/ctvm.o

ctvmlib_static: $(DEPS)
		# Compile both of the libraries to their respective object files
		$(CXX) -Wall $(CPPFLAGS) -o $(SRC_DIR)/ctvm_util.o -c $(SRC_DIR)/ctvm_util.cpp
		$(CXX) -Wall $(CPPFLAGS) -o $(SRC_DIR)/ctvm.o -c $(SRC_DIR)/ctvm.cpp

		# Build static libraries
		ar rcs $(LIB_DIR)/libctvm_util.a $(SRC_DIR)/ctvm_util.o
		ar rcs $(LIB_DIR)/libctvm.a $(SRC_DIR)/ctvm.o

$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp $(LIB_DIR)/cygctvm.dll $(LIB_DIR)/cygctvm_util.dll
		$(CXX) $(CPPFLAGS) -c -o $@ $< 

test1: $(TEST_DIR)/test1.o
		$(CXX) -Llib $(LDFLAGS) -lctvm -lctvm_util -o $(BIN_DIR)/test1 $(TEST_DIR)/test1.o 

test_cs: $(TEST_DIR)/test_cs.o
		$(CXX) -Llib $(LDFLAGS) -lctvm -lctvm_util -o $(BIN_DIR)/test-cs $(TEST_DIR)/test_cs.o 	

executable: ctvmlib_static
		$(CXX) -Wall $(CPPFLAGS) -o $(SRC_DIR)/ctvm_recover.o -c $(SRC_DIR)/ctvm_recover.cpp
		$(CXX) -Llib -lctvm -lctvm_util $(LDFLAGS) -o $(BIN_DIR)/ctvm-recover $(SRC_DIR)/ctvm_recover.o

clean:
		rm -f $(BIN_DIR)/*
		rm -f $(SRC_DIR)/*.o
		rm -f $(TEST_DIR)/*.o
		rm -f $(LIB_DIR)/*

test: clean all
		$(BIN_DIR)/test1 test/data/peppers.jpg $(HOME)/tmp/output.jpg
		$(BIN_DIR)/ctvm-recover test/data/testSino.png test/data/testAngles.dat a3
		$(BIN_DIR)/test-cs 32 test/data/box.png ~/tmp/test-cs-out.png


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