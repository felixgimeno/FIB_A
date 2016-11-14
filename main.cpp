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
#include <ctime>
#include <fstream>
#include <iomanip>
#include <thread>
#include <time.h> 
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

typedef set<string> k_shingles; //set de ksingles
typedef vector<k_shingles> collection_of_k_shingles;
typedef map<size_t, string> row_to_string;
typedef vector<vector<size_t> > signature_matrix;
typedef function<size_t (size_t) > hash_function;
typedef function<size_t (vector<size_t>) > hash_function_for_vectors;
typedef vector< hash_function > vector_of_hash_functions;
typedef vector< hash_function_for_vectors > vector_of_hash_function_for_vectors;

// GLOBAL VARS

vector<string> textos;
const int debug = 0; // jaccard results without times
const int debug_time = 1; // only times and similarty values
const int generation = true; // generations of documents
const int ndocs = 50; // number of documents
const int nwords = 50; //number of word/documents


void generating_docs() {
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
	/*clock_t time_start, time_end;
	time_start = clock();*/
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
	/*time_end = clock();
	cout << "Jaccard sim. elapsed " << ((float)(time_end - time_start))/CLOCKS_PER_SEC << " seconds" << endl;*/
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
	/*clock_t time_start, time_end;
	time_start = clock();*/
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
	/*time_end = clock();
	cout << "compute signature matrix elapsed " << ((float)(time_end - time_start))/CLOCKS_PER_SEC << " seconds" << endl;*/
	return ret;
}


/*
 * jaccard similarity given a signature matrix and the indexes of two k_shingles
 */
double jaccard_similarity(const size_t index1, const size_t index2, const signature_matrix& sm){
	/*clock_t time_start, time_end;
	time_start = clock();*/
	size_t acc = 0; 
	for (size_t j = 0; j < sm.size(); j += 1){
		acc += int(sm.at(j).at(index1) == sm.at(j).at(index2));
	}
	/*time_end = clock();
	cout << "sim (document) elapsed " << ((float)(time_end - time_start))/CLOCKS_PER_SEC << " seconds" << endl;*/
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
hash_function get_hash_function(const size_t a, const size_t c , const size_t p ){
	return [a,c,p ](size_t x) -> size_t {
		return (a*x+c) %p ;
	};
}

hash_function_for_vectors get_hash_function_for_vectors(const size_t a , const size_t c , size_t p){
	return [a,c,p ](vector<size_t> x) -> size_t {
		size_t ret = 0;
		for (size_t k : x){
			ret = ( ret * a + c + k ) %p ;  
			}
		return ret;
	};	
}

/*
 * locality sensitive hashing
 * returns a set of all candidate pairs of similar columns
 */
 
set<pair<size_t, size_t> > lsh(signature_matrix& sm, const vector_of_hash_function_for_vectors& vf, float t, size_t r){
	/*clock_t time_start, time_end;
	time_start = clock();*/
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
	/*time_end = clock();
	cout << "LSH elapsed " << ((float)(time_end - time_start))/CLOCKS_PER_SEC << " seconds" << endl;*/
	return pair_candidates;
}

double absolute(const double s){return s > 0 ? s : -s;}

double get_relative_error(const collection_of_k_shingles& cks, const signature_matrix& sm, const size_t index1, const size_t index2){
	const double js_sm = jaccard_similarity(index1, index2, sm);
	const double js = jaccard_similarity(cks.at(index1), cks.at(index2));
	const double epsilon = 0.0001;
	return absolute((js-js_sm)/(absolute(js)+epsilon));
}

/*
 * returns the square mean of the errors between the two methods for jaccard similarity 
 */
double test(const collection_of_k_shingles& cks, const signature_matrix& sm){
	double ret = 0;
	size_t count = 0;
	for (size_t u = 0; u < cks.size(); u += 1){
		for (size_t v = u+1; v < cks.size(); v += 1){
			const double err = get_relative_error(cks, sm, u, v);
			ret += err*err;
			count += 1;
		}
	}
	return sqrt(ret/count);
}





string get_text_file(string filename){
	string str = "";
	vector<string> vwords;
	ifstream file;
	file.open(filename);
	if (file.is_open()) {
		string s;
		while (!file.eof()) {
			file >> s;
			vwords.push_back(s);
		}
	}
	file.close();
	//CONVERTING TO STRING.
	for (int j = 0; j < int(vwords.size()); ++j) {
		str += vwords[j];
		str += " "; //space beetwen words.
		}

	//cout  << str << endl;	
	return str;
}



void reading_documents() {
	for (int i = 1; i <= ndocs; ++i) {
		stringstream stream; 
        string palabra; 
	    stream << i; 
        palabra = stream.str();
		palabra +=".txt";
		textos.push_back(get_text_file(palabra.c_str()));
	}
  }
/*
 * http://primes.utm.edu/prove/prove2_3.html 
 */
size_t get_prime(const size_t lower_bound){
	function<size_t(size_t,size_t,size_t)> pow = [&pow](size_t a, size_t b, size_t c)->size_t {
		if (b == 0) {return 1;}
		size_t res = pow(a,b/2,c);
		size_t res2 = res*res % c;
		if (b % 2 == 0){
			return res2;
			}
		return res2*a%c; };
	const vector<size_t> list = {2, 3, 5, 7, 11, 13, 17};
	auto atomic_test = [&pow](size_t a, size_t n, size_t d, size_t two_to_r)->bool{
		const size_t res = pow(a,d,n);
		if ( res == 1 ){return true;}
		if ((pow(res, two_to_r,n)+1)%n == 0){return true;}
		return false;
		};
	auto mytest = [&list,&atomic_test](size_t n) -> bool {
		const size_t m = n-1;
		size_t d = m;
		size_t two_to_r = 1;
		while (d % 2 == 0) {
			d = d/2;
			two_to_r *= 2;
			} 
		for(size_t v : list){
			if (not(atomic_test(v,n,d,two_to_r))){
				return false;
				}
			}
		return true;
		};
	for (size_t u = lower_bound; ; u += 1){
		if (u%2 == 1 and mytest(u)){
			return u;
			}
		}
	return 0;}

vector_of_hash_functions get_vector_of_hash_functions(const size_t number_of_hash_functions, const size_t lower_bound){
	vector_of_hash_functions vh = {};
	srand(time(0));
	for (size_t q = 0; q < number_of_hash_functions ; q += 1){
		vh.push_back(get_hash_function(rand(),rand(),get_prime(lower_bound/*+rand()*/)));
	}
	return vh;	
}
vector_of_hash_function_for_vectors get_vector_of_hash_function_for_vectors(const size_t number_of_vector_of_hash_function_for_vectors, const size_t lower_bound){
	vector_of_hash_function_for_vectors vf = {};
	srand(time(0));
	for (size_t q = 0; q < number_of_vector_of_hash_function_for_vectors ; q += 1){
		vf.push_back(get_hash_function_for_vectors(rand(),rand(),get_prime(lower_bound/*+rand()*/)));
	}
	return vf;	
}



int main(void) {
	clock_t time_start, time_end;
    const size_t number_of_hash_functions = 100; 
	const size_t number_of_vector_of_hash_function_for_vectors = 100;  
	const vector<float>  t_values = {0.1,0.2,0.3,0.4,0.5 ,0.6,0.7,0.8,0.9,1}; // Threshold values
	//Reading Documents
     if (generation) {generating_docs();}
     reading_documents();
	
	//ahora creamos las funciones de hash aleatorias que simularan las permutaciones
	const size_t lower_bound = 1299827; //https://primes.utm.edu/lists/small/100000.txt
	clock_t start_vector_of_hash_functions = clock();
	vector_of_hash_functions vh = get_vector_of_hash_functions(number_of_hash_functions, lower_bound);
	clock_t end_vector_of_hash_functions = clock();
	if (debug_time){cout << "time to compute vector_of_hash_functions "  << ((double)(end_vector_of_hash_functions - start_vector_of_hash_functions))/CLOCKS_PER_SEC << endl;}
	
	
	//ahora creamos las funciones que hashearan las bandas del LSH
	clock_t start_vector_of_hash_function_for_vectors = clock();
	vector_of_hash_function_for_vectors vf = get_vector_of_hash_function_for_vectors(number_of_vector_of_hash_function_for_vectors, lower_bound);
	clock_t end_vector_of_hash_function_for_vectors = clock();
	if (debug_time){cout << "time to compute start_vector_of_hash_function_for_vectors "  << ((double)(end_vector_of_hash_function_for_vectors - start_vector_of_hash_function_for_vectors))/CLOCKS_PER_SEC << endl;}
    
    for (size_t k = 5; k < 10; ++k){
		//ahora hacemos los calculos que solo dependen de k y number_of_hash_functions
		collection_of_k_shingles cks;
		for (string s : textos){
			cks.push_back(get_k_shingles_from_string(k, s));
		}
		
		time_start = clock();
		signature_matrix sm = compute_signature_matrix(cks, vh);
		time_end = clock();
		if(debug_time) cout << "csm	" << k << 0 << "	" << 0 << "	" << ((float)(time_end - time_start)/CLOCKS_PER_SEC) << endl;
		
		//if(debug) cout << "Para la k " << k << " la media cuadrada de los errores es: " << test(cks, sm) << endl;
		//fin de los calculos que solo dependen de k y number_of_hash_functions
		const bool experimento_lsh = true;
		if (experimento_lsh){		
		for (float t : t_values) {
			for (size_t r = 1; r < 7; ++r) {
				//COMPUTATION OF JACCARD & LSH SIMILARITIES
				if(debug) cout << "Candidatos según LSH constantes k " << k << " r " << r << " t " << t << endl;
				time_start = clock();
				set<pair<size_t, size_t> > pairconcidence =  lsh(sm,vf,t,r);
				time_end = clock();
				if(debug_time) cout << "lsh	" << k << "	" << r << "	"  << t << "	" <<  ((float)(time_end - time_start)/CLOCKS_PER_SEC) << "	" << pairconcidence.size() << endl;
				for (auto& kp : pairconcidence){ 
					if(debug) cout << " documentos " << kp.first << " " << kp.second; 
					time_start = clock();
					const double js_sm = jaccard_similarity(kp.first, kp.second, sm);
					time_end = clock();
					if(debug_time) cout << "jacsim_sm	" << k << "	" << r << "	"  << t << "	" << ((float)(time_end - time_start)/CLOCKS_PER_SEC) << endl;
					time_start = clock();
					const double js = jaccard_similarity(cks.at(kp.first), cks.at(kp.second));
					time_end = clock();
					if(debug_time) cout << "jacsim_cks	" << k << "	" << r << "	"  << t << "	" << ((float)(time_end - time_start)/CLOCKS_PER_SEC) << endl;
					if(debug) cout << " jaccard sm " <<  js_sm;
					if(debug) cout << " jaccard directo " << js;
					if(debug) cout << " error relativo " << absolute((js-js_sm)/(js +  0.000001));
					if(debug) cout << endl;
				}
				time_end = clock();
				if(debug) cout << "total elapsed " << ((float)(time_end - time_start))/CLOCKS_PER_SEC << " seconds" << endl;
			  

			  }
		   }
		}
	 }
} 
