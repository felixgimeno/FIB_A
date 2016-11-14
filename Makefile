all: main

main: main.cpp
	g++ -Wall -Wextra -std=c++11  -o main main.cpp



clean:
	-rm ./main
