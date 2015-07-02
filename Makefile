INCLUDE_DIR=include
SRC_DIR=src
BIN_DIR=bin
TEST_DIR=test
LIBDIR=/usr/local/lib
BOOST_INC=/usr/local/include

CXX=g++
CPPFLAGS=-I$(INCLUDE_DIR) -I$(BOOST_INC)
LDFLAGS=-L$(LIBDIR)
LDLIBS=-lboost_system
DEPS=$(INCLUDE_DIR)/ctvm.h


all: test1

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

$(TEST_DIR)/%.o: $(TEST_DIR)/%.c $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

test1: $(TEST_DIR)/test1.o
	$(CXX) -o $(BIN_DIR)/test1 $(TEST_DIR)/test1.o $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(BIN_DIR)/*
	rm -f $(SRC_DIR)/*.o
	rm -f $(TEST_DIR)/*.o
