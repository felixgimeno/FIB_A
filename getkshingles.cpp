#include <iostream>
#include <string>
#include <set>
using namespace std;

typedef set<string> k_shingles;

/*
 *  str shouldn't have newlines or tabs or multiple consecutive spaces 
 */ 
k_shingles get_k_shingles_from_string(const size_t k, const string str){
	k_shingles ret = set<string> ();
	for (size_t i = 0; i < str.size() ; i+= 1){
		if (i >= k){
			ret.insert(str.substr(i-k,k));
		}
	}
	return ret;
}

double jaccard_similarity(const k_shingles A, const k_shingles B){
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

void print(const k_shingles A){for (const auto& v : A){cout << "\"" << v << "\", ";} cout << endl;}

int main(void){
	const size_t k = 9;
	const string str1 = "minhashing is the best";
	const string str2 = "minhashing is not the best";
	const k_shingles A = get_k_shingles_from_string(k, str1);
	const k_shingles B = get_k_shingles_from_string(k, str2);
	print(A);
	print(B);
	const double j_sim = jaccard_similarity(A, B);
	cout << j_sim << endl;
	return 0;
	
} 
