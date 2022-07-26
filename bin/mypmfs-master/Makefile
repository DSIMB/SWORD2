CC      = g++
CFLAGS  = -std=c++0x -Wall -Wextra -O3 -march=native
SRC1    = src/training.cpp
SRC2    = src/scoring.cpp
EXE1    = training
EXE2    = scoring

all: $(EXE1) $(EXE2)
	
$(EXE1): $(SRC1)
	$(CC) $(SRC1) $(CFLAGS) -o $(EXE1)

$(EXE2): $(SRC2)
	$(CC) $(SRC2) $(CFLAGS) -o $(EXE2)

