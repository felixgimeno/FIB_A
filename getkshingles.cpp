#include <functional>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <random>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <sstream> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>
using namespace std;

typedef set<string> k_shingles; //set de ksingles
typedef vector<k_shingles> collection_of_k_shingles;
typedef map<size_t, string> row_to_string;
typedef vector<vector<size_t> > signature_matrix;
typedef function<size_t (size_t) > hash_function;
typedef function<size_t (vector<size_t>) > hash_function_for_vectors;
typedef vector< hash_function > vector_of_hash_functions;
typedef vector< hash_function_for_vectors > vector_of_hash_function_for_vectors;

/*
 *  given k returns set of all substrings of length k that appear in str 
 *  str shouldn't have newlines or tabs or multiple consecutive spaces
 */ 
k_shingles get_k_shingles_from_string(const size_t& k, const string& str){
	k_shingles ret = k_shingles ();
	for (size_t i = 0; i < str.size() ; i+= 1){
		if (i >= k){
			ret.insert(str.substr(i-k,k));
		}
	}
	return ret;
}

/*
 * computes the jaccard similarity from two sets
 * jaccard similarity of two sets is the number of elements in both divided by the number of elements in one or the other
 */
double jaccard_similarity(const k_shingles& A, const k_shingles& B){
	k_shingles C = B;
	uint_fast64_t count_intersection = 0;
	uint_fast64_t count_union = 0;
	for (const auto& v : A){
		const auto it = C.find(v);
		const bool t = it == C.end(); 
		if (t){
			// v is in A, isn't in B
			count_intersection += 0;
			count_union += 1;
		}
		else {
			// v is in A, is in B
			count_intersection += 1;
			count_union += 1;
			// we remove from our copy of B 
			// the intersection with A
			C.erase(it);			
		}
	}
	count_union += C.size();
	C.clear();
	return double(count_intersection)/double(count_union);	
}

/*
 * (for debugging) prints the set of k-shingles
 */
void print(const k_shingles& A){for (const auto& v : A){cout << "\"" << v << "\", ";} cout << endl;}

/*
 * each row in the characteristic matrix is one of the k_shingles cointained in any of the sets in cks
 * this function returns a map of row number to k-shingle addigned to that row
 */
row_to_string get_rows(const collection_of_k_shingles& cks){
	row_to_string ret = row_to_string ();
	k_shingles rows = k_shingles ();
	size_t i = 0;
	for (const k_shingles& ks : cks){
		for (const string& s : ks){
			if (rows.find(s) == rows.end()){
				rows.insert(s);
				ret[i] = s;
				i += 1; 
			}
		}
	}
	return ret;	
}
/*
 * section 3.3.5 of the book
 * computes a signature matrix from which the jaccard similarity of two set can be approximated
 */
signature_matrix compute_signature_matrix(const collection_of_k_shingles& cks, const vector_of_hash_functions& vf){
	signature_matrix&& ret = signature_matrix (vf.size(), vector<size_t> (cks.size(), numeric_limits<size_t>::max()));
	row_to_string&& rows = get_rows(cks);
	const size_t number_of_rows = rows.size(); 
	for (const auto& kv : rows){
		vector<size_t> vfr = vector<size_t> (vf.size(), kv.first);
		transform(vfr.begin(), vfr.end(), vf.begin(), vfr.begin(), [number_of_rows](size_t row_number, function<size_t (size_t)> f){return f(row_number)%number_of_rows;}); 
		for (size_t i = 0; i < cks.size(); i += 1){
			const bool column_has_zero_in_this_row = cks.at(i).find(kv.second) ==  cks.at(i).end();
			if (not column_has_zero_in_this_row){
				for (size_t j = 0; j < vf.size(); j += 1){
					ret.at(j).at(i) = min(ret.at(j).at(i), vfr.at(j));
				}
			}
		}
	}
	return ret;
}


/*
 * jaccard similarity given a signature matrix and the indexes of two k_shingles
 */
double jaccard_similarity(const size_t index1, const size_t index2, signature_matrix& sm){
	size_t acc = 0; 
	for (size_t j = 0; j < sm.size(); j += 1){
		acc += int(sm.at(j).at(index1) == sm.at(j).at(index2));
	}
	return acc/double(sm.size());
}

/*
 * (for debugging) print signature matrix
 */
void print(const signature_matrix& sm){
	cout << "printing signature matrix:\n";
	for (auto& r : sm){
		for (auto& e : r){
			cout << e << " "; 
		}
		cout << "\n";
	}
	cout << "end of signature matrix\n";
}

/*
 * this is a universal hashing framework from which we can obtain hash functions
 * https://en.wikipedia.org/wiki/Universal_hashing#Hashing_integers
 * p should be prime >= m
 */
hash_function get_hash_function(size_t a, size_t c /* , size_t p */){
	return [a,c /* ,p */](size_t x) -> size_t {
		return (a*x+c)/* %p */;
	};
}

hash_function_for_vectors get_hash_function_for_vectors(size_t a , size_t c /* , size_t p */){
	return [a,c /* ,p */](vector<size_t> x) -> size_t {
		size_t ret = 0;
		for (size_t k : x){
			ret = ( ret * a + c + k ) /* %p */ ;  
			}
		return ret;
	};	
}

/*
 * locality sensitive hashing
 * returns a set of all candidate pairs of similar columns
 */
 
set<pair<size_t, size_t> > lsh(signature_matrix& sm, const vector_of_hash_function_for_vectors& vf, float t, size_t r){
	
	map< pair<size_t, size_t> , size_t> coincidencias; //si quereis cambiad el tipo
	size_t nbands = int(ceil(sm.size()/r));
	for (size_t i = 0; i < nbands; i += 1){
		map<size_t, set<size_t> > bucket = map<size_t, set<size_t> > ();
		for(size_t j = 0; j < sm.at(0).size(); j += 1){
			//i numero banda, j numero columna, vf.size() es b, t es el threshold
			vector<size_t> input_hash_function = vector<size_t> ();
			for(size_t q = 0; q < r; q += 1){
				input_hash_function.push_back(sm.at(i*r+q).at(j));
			}
			size_t output_hash_function = vf.at(i % vf.size())(input_hash_function);
			bucket[output_hash_function].insert(j);
		} 
		//miramos los "buckets" y si coinciden añadimos 1 a "coincidencias"
		for (auto itmap = bucket.begin(); itmap != bucket.end(); ++itmap){ //para cada set del bucket
			for (auto itseta = (itmap->second).begin(); itseta != (itmap->second).end(); ++itseta){ //fijamos el primer elemento del set
				for (auto itsetb = itseta; itsetb != (itmap->second).end(); ++itsetb){ //recorremos los demás elementos a partir del fijo
					if(*itseta < *itsetb){
						coincidencias[pair<size_t, size_t> (*itseta, *itsetb)]++;
					}
					if(*itseta == *itsetb){}
					if(*itseta > *itsetb){
						coincidencias[pair<size_t, size_t> (*itsetb, *itseta)]++;
					}
				}
			}
		}
	}
	//recorremos coincidencias y ponemos en un set las parejas que tengan más de ? coincidencias		
	set<pair<size_t, size_t> > pair_candidates = set<pair<size_t, size_t> > ();
	for(auto it: coincidencias){
		//determine whether the fraction of components in which they agree is at least t
		if(it.second >= t * nbands){
			pair_candidates.insert(it.first);
		}
	}
	return pair_candidates;
}

double absolute(const double s){return s > 0 ? s : -s;}

int main(void){
	const clock_t begin_time = clock();
    
    vector<string> textos;
    //vamos leyendo los  20 documentos y lo vamos metiendo en el documento
	for (int i = 1; i <= 20; ++i) {
		
		string str = "";
		vector<string> vwords;
		ifstream file;
		
		stringstream stream; 
        string palabra; 
	    stream << i; 
        palabra = stream.str();
		palabra +=".txt";
		
		
		file.open(palabra.c_str());
		if (file.is_open()) {
		string s;
		//leemos del vector.
		while (!file.eof()) {
			file >> s;
			vwords.push_back(s);
			}
		}
		file.close();
		//pasamos el word iessimo a un solo string con sus espacios.
		for (int j = 0; j < int(vwords.size()); ++j) {
			str += vwords[j];
			str += " "; //ojo con el ultimo espacio de la ultima word.
			}
		//print debugger
		cout << "el string of the file  " << i << " is " << endl;
		cout  << str << endl;
		textos.push_back(str);
	
		}
	
	//ahora creamos las funciones de hash aleatorias que simularan las permutaciones
	const size_t number_of_hash_functions = 100;
	vector_of_hash_functions vh = {[](size_t a){return a;}};
	srand(time(0));
	for (size_t q = 0; q < number_of_hash_functions ; q += 1){
		vh.push_back(get_hash_function(rand(),rand()));
	}
	const size_t number_of_vector_of_hash_function_for_vectors = 100;
	vector_of_hash_function_for_vectors vf = vector_of_hash_function_for_vectors ();
	for (size_t q = 0; q < number_of_vector_of_hash_function_for_vectors ; q += 1){
		vf.push_back(get_hash_function_for_vectors(rand(),rand()));
	}
	const vector<float>  t_values = {0.1,0.2,0.3,0.4,0.5 ,0.6,0.7,0.8,0.9,1};	
    for (size_t k = 1; k < 3; ++k){
		
		for (float t : t_values) {
			for (size_t r = 1; r < 3; ++r) {
				collection_of_k_shingles cks;
				
				for (int i = 0; i < 20; ++i){
					k_shingles A = get_k_shingles_from_string(k, textos[i]);
					//print(A);
					//lo metemos en la coleccion
					cks.push_back(A);
				}
				signature_matrix sm = compute_signature_matrix(cks, vh); //no tocar
				
				for (int i = 0; i < int(cks.size()); ++i) {
					for (int j = i+1;j < int(cks.size()); ++j) {
						//cout << "jaccard similarity " << i+1 << " y " << j+1 << " es: " << jaccard_similarity(cks.at(i),cks.at(j)) << endl; 
						}
					//cout << endl;
				}
				//cout << endl;
				//cout << "ahora calculamos la variante de jaccard_similarity(i,j,sm)" << endl;
				//cout << "de momento no funciona la jaccard_similarity(i,j,sm) por eso da zero's siempre" << endl;	
				
				for (int i = 0; i < int(cks.size()); ++i) {
					for (int j = i + 1; j < int(cks.size()); ++j) {
						//cout << "jaccard similarity " << i+1 << " y " << j+1 << " es " <<  jaccard_similarity(i,j,sm) << endl; 
						}
					//cout << endl;
				}
				//Mirate este trozo de codigo.
				//solo quiero obtener el set de lsh y imprimrlo , pero hace nada
				//compila y tal ejecuta , pero el lsh no lo esta haciendo
				cout << "Candidatos según LSH constantes k " << k << " r " << r << " t " << t ; // << endl;
				set<pair<size_t, size_t> > pairconcidence =  lsh(sm,vf,t,r);
				for (auto& k : pairconcidence){ 
					cout << " documentos " << k.first << " " << k.second; // << endl; 
					//cout << "estoy calculando la del lsh" << endl;
					const double js_sm = jaccard_similarity(k.first, k.second, sm);
					const double js = jaccard_similarity(cks.at(k.first), cks.at(k.second));
					cout << " jaccard sm " <<  js_sm;// << endl;
					cout << " jaccard directo " << js;
					cout << " error relativo " << absolute((js-js_sm)/js);
					cout << endl;
					
					}
		
				
			  }
		   }
		}
	
    cout << "elapsed time " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
	
  } 
