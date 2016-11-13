#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <string>
#include <sstream> 
#include <stdio.h>
using namespace std;



int main() {
	srand(time(NULL));
	int nwords = 50;
	
	vector<string> vwords;
	vector<int> vlength;
	ifstream file;
	file.open("500words_avglong3.txt");
	
	if (file.is_open()) {
	string s;
	//leemos del vector.
	while (!file.eof()) {
		file >> s;
		vwords.push_back(s);
		vlength.push_back(s.length());
		}
	}
	file.close();
	float avglen = 0;
	for(int i=0; i<vlength.size(); ++i){
		avglen += vlength[i];
	}
	cout << "longitud media de las palabras: " << (avglen/vlength.size()) << endl;
	//calculando las 20 permutaciones randoms, tomando como base 50words
	//cada vector vec_s, repsentanda el contenido de un archivo , habra
	// 20 txt's totales.
	
	cout << "calculando 20 permutaciones aleatorias" << endl;
	for(int i = 1; i <= 20; ++i) {
	    stringstream stream; 
        string palabra; 
	    stream << i; 
        palabra = stream.str();
		palabra +=".txt";
		
		ofstream outfile;
		outfile.open(palabra.c_str());
		
		for (int k = 0; k < nwords; ++k) {
				int num=rand()%nwords;
				string s = vwords[num];
				outfile << s << " ";
			}
		outfile.close();
		}
	 }
