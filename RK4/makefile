CC = g++
PY = python
FLAGS = -std=c++11 -O3

all: main exec data

main: main.cpp
	$(CC) $(FLAGS) -o $< $@

exec: main
	time ./$<

data: plot.py main
	$(PY) $< 
