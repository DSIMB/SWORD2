OS:=$(shell uname)
CFLAGS  = -fopenmp -std=c++0x -Wall -Wextra -O3 -march=native
ifeq ($(OS),Darwin)
    CFLAGS  = -std=c++0x -Wall -Wextra -O3
endif
CC      = g++
SRC1    = src/training.cpp
SRC2    = src/scoring.cpp
SRC3    = src/scoring_omp.cpp
EXE1    = training
EXE2    = scoring
EXE3    = scoring_omp

all: $(EXE1) $(EXE2) $(EXE3)
	
$(EXE1): $(SRC1)
	$(CC) $(SRC1) $(CFLAGS) -o $(EXE1)

$(EXE2): $(SRC2)
	$(CC) $(SRC2) $(CFLAGS) -o $(EXE2)

$(EXE3): $(SRC3)
	$(CC) $(SRC3) $(CFLAGS) -o $(EXE3)
