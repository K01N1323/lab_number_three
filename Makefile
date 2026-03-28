CXX = clang++
CXXFLAGS = -std=c++17 -Wall -IBaseStructures -ISequences -IComplexSequences -ILinearAlgebra -IUI -ITests

# Автоматически находим все .cpp файлы в корне и во всех папках
SRCS = $(wildcard *.cpp) $(wildcard */*.cpp)
OBJS = $(SRCS:.cpp=.o)
TARGET = main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
