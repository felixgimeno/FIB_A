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


int ndocs = 20; // numero de documentos
int nwords = 50; //numero de palabras por documentos

int main() {
	srand(time(NULL));
	vector<string> vwords;
	vector<int> vlength;
	ifstream file;
	file.open("50words.txt");
	
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
	for(int i=0; i<int(vlength.size()); ++i){
		avglen += vlength[i];
	}
	cout << "longitud media de las palabras: " << (avglen/vlength.size()) << endl;
	//calculando las 20 permutaciones randoms, tomando como base 50words
	//cada vector vec_s, repsentanda el contenido de un archivo , habra
	// 20 txt's totales.
	
	cout << "Numero de documentos generados: " << ndocs << endl;
	cout << "Numero de palabras por documentos generados: " << nwords << endl;
	for(int i = 1; i <= ndocs; ++i) {
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
