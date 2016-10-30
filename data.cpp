#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
using namespace std;


//leer vector
void leer_vec(vector<string> &v) {
	int n = v.size();
	for (int i = 0; i < n; ++i) {
		int n;
		cin >> n;
		v[i] = n; 
		}
	}



//escrirbir vector
void escribir_vec(vector<string> &v) {
	int n = v.size();
	for (int i = 0; i < n; ++i) {
		cout << v[i] << " "; 
		}
	cout << endl;
	}




int main() {
	
	ifstream file;
	
	vector<string> vwords;
	
	file.open("50words.txt");
	
	if (file.is_open()) {
	string s;
	//leemos del vector.
	while (!file.eof()) {
		file >> s;
		vwords.push_back(s);
		}
	}
	
	file.close();
	escribir_vec(vwords);
	cout << endl;
	
	
	//calculando las 20 permutaciones randoms, tomando como base 50words
	//cada vector vec_s, repsentanda el contenido de un archivo , habra
	// 20 txt's totales.
	cout << "calculando 20 permutaciones aleatorias" << endl;
	for(int i = 0; i < 20; ++i) {
	   cout << "calculando la permutacion " << i + 1 << endl << endl;
	   vector<string> vec_s;
	   for (int j = 0; j < 50; ++j) {
		int num=rand()%50;
		string s = vwords[num];
		vec_s.push_back(s);
		} 
		escribir_vec(vec_s);
		cout << endl << endl;
		
		}
	 }
