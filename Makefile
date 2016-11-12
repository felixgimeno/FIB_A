all: data getkshingles

data: data.cpp
	g++ -Wall -Wextra -std=c++11  -o data data.cpp

getkshingles: getkshingles.cpp
	g++ -Wall -Wextra -std=c++11  -o getkshingles getkshingles.cpp

clean:
	-rm ./data ./getkshingles ./data.exe ./getkshingles.exe