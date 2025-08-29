CXX = g++
CXXFLAGS = -std=c++17 -O3 -g -fopenmp
TARGET = Staple_Insertion
SRC = main.cpp io.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) $(TARGET).o