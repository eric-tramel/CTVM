INCLUDE_DIR=include
SRC_DIR=src
BIN_DIR=bin
TEST_DIR=test
LIB_DIR=lib
MK_BOOST_LIB=$(BOOST_LIBDIR)
MK_BOOST_INC=$(BOOST_INCDIR)

CXX=g++
CPPFLAGS=-I$(INCLUDE_DIR) -I$(MK_BOOST_INC)
LDFLAGS=-L$(MK_BOOST_LIB)
LDLIBS=-lboost_system
DEPS=$(INCLUDE_DIR)/ctvm.h


all: checkdir ctvmlib test1

# $(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(DEPS)
# 	$(CXX) -c -o $@ $< $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

ctvmlib: $(DEPS)
	$(CXX) -dynamiclib -fPIC -o $(LIB_DIR)/libctvm.dylib $(SRC_DIR)/ctvm.cpp $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

# $(LIB_DIR)/%.dylib: $(SRC_DIR)/%.c $(DEPS)
# 	$(CXX) -dynamiclib -fPIC -o $@ $< $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

$(TEST_DIR)/%.o: $(TEST_DIR)/%.c $(LIB_DIR)/ctvm.dylib
	$(CXX) -c -o $@ $< $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

test1: $(TEST_DIR)/test1.o
	$(CXX) -o $(BIN_DIR)/test1 $(TEST_DIR)/test1.o $(LDFLAGS) $(LDLIBS) -Llib -lctvm

clean:
	rm -f $(BIN_DIR)/*
	rm -f $(SRC_DIR)/*.o
	rm -f $(TEST_DIR)/*.o
	rm -f $(LIB_DIR)/*

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