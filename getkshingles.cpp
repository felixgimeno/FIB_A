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

int main(void){
	return 0;
} 
