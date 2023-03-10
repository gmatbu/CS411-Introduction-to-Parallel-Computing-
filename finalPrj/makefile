CC = g++

CFLAGS = --std=c++11 -fopenmp -lm
  
  all: main.cpp
	$(CC) $(CFLAGS) -o main main.cpp
  
  debug: main.cpp
	$(CC) $(CFLAGS) -g -o main main.cpp

  clean: 
	rm main	