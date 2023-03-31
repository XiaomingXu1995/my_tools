#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "CLI11.hpp"
#include "common.hpp"
#include <string>
#include <sys/sysinfo.h>
#include <fstream>
#include <sstream>
#include "robin_hood.h"
//#include <unistd.h>

using namespace std;


int main(int argc, char* argv[]){
	CLI::App app("cal-sse-Jaccard: calculate the Sum of squared error between estimated Jaccard and true Jaccard");
		
	string standard_file = "";
	string rabbit_file = "";
	string mash_file = "";
	string bindash_file = "";
	string dashing2_file = "";
	string output_prefix = "output";

	//int threads = get_nprocs_conf();
	
	auto option_standard_file = app.add_option("--standard", standard_file, "set the standard distance file");
	auto option_rabbit_file = app.add_option("--rabbit", rabbit_file, "set the RabbitKSSD distance file");
	auto option_mash_file = app.add_option("--mash", mash_file, "set the mash distance file");
	auto option_bindash_file = app.add_option("--bindash", bindash_file, "set the bindash distance file");
	auto option_dashing2_file = app.add_option("--dashing2", dashing2_file, "set the dashing2 distance file");
	app.add_option("-o, --output", output_prefix, "set the output prefix file");

	option_standard_file->required();
	
	CLI11_PARSE(app, argc, argv);

	if(!*option_rabbit_file && !*option_mash_file && !*option_bindash_file && !*option_dashing2_file){
		cerr << "Must have an estimated distance file, try --rabbit, --bindash, or --dashing2" << endl;
		return 1;
	}

	//parse the standard_file
	ifstream ifs_standard;
	ifs_standard.open(standard_file);
	if(!ifs_standard){
		cerr << "ERROR: cannot open standard_file: " << standard_file << endl;
		return 1;
	}
	string line;
	robin_hood::unordered_map<string, double> name_jaccard_map;
	while(getline(ifs_standard, line)){
		stringstream ss;
		string query_name, ref_name, common_size;
		double jaccard, mash_distance;
		ss << line;
		ss >> query_name >> ref_name >> common_size >> jaccard >> mash_distance;
		string merged_name = query_name + ref_name;
		name_jaccard_map.insert({merged_name, jaccard});
	}
	ifs_standard.close();

	if(*option_rabbit_file){
		cerr << "=====run the option rabbit file " << endl;
		ifstream ifs_rabbit;
		ifs_rabbit.open(rabbit_file);
		if(!ifs_rabbit){
			cerr << "ERROR: cannot open rabbit_file: " << rabbit_file << endl;
			return 1;
		}
		string output_rabbit = output_prefix + ".rabbit.res";
		ofstream ofs_rabbit(output_rabbit);
		ofs_rabbit << "query_name\tref_name\tstandard_jaccard\testimated_jaccard\terror" << endl;
		double sse_rabbit = 0.0;
		while(getline(ifs_rabbit, line)){
			stringstream ss;
			string query_name, ref_name, common_size;
			double jaccard, mash_distance;
			ss << line;
			ss >> query_name >> ref_name >> common_size >> jaccard >> mash_distance;
			string merged_name = query_name + ref_name;
			if(name_jaccard_map.count(merged_name) == 0){
				cerr << "ERROR: the merged_name: " << merged_name << " is not in the standard name_jaccard_map" << endl;
				return 1;
			}
			double standard_jaccard = name_jaccard_map[merged_name];
			double cur_err = jaccard - standard_jaccard;
			sse_rabbit += cur_err * cur_err;
			ofs_rabbit << query_name << '\t' << ref_name << '\t' << standard_jaccard << '\t' << jaccard << '\t' << cur_err << endl;
		}
		ifs_rabbit.close();
		ofs_rabbit.close();
		cout << "the Sum of Squared Error (SSE) of rabbit_dist is: " << sse_rabbit << endl;
	}
	if(*option_mash_file){
		cerr << "=====run the option mash file " << endl;
		ifstream ifs_mash;
		ifs_mash.open(mash_file);
		if(!ifs_mash){
			cerr << "ERROR: cannot open mash_file: " << mash_file << endl;
			return 1;
		}
		string output_mash = output_prefix + ".mash.res";
		ofstream ofs_mash(output_mash);
		ofs_mash << "query_name\tref_name\tstandard_jaccard\testimated_jaccard\terror" << endl;
		double sse_mash = 0.0;
		while(getline(ifs_mash, line)){
			stringstream ss;
			string query_name, ref_name, jaccard_str;
			double mutation_distance, p_value, jaccard;
			ss << line;
			ss >> ref_name >> query_name >> mutation_distance >> p_value >> jaccard_str;
			string merged_name = query_name + ref_name;
			if(name_jaccard_map.count(merged_name) == 0){
				cerr << "ERROR: the merged_name: " << merged_name << " is not int the standard name_jaccard_map" << endl;
				return 1;
			}
			double standard_jaccard = name_jaccard_map[merged_name];
			int mid_pos = jaccard_str.find_first_of('/');
			int inter_size = stoi(jaccard_str.substr(0, mid_pos));
			int union_size = stoi(jaccard_str.substr(mid_pos+1, jaccard_str.length()-mid_pos-1));
			jaccard = (double)inter_size / union_size;
			double cur_err = jaccard - standard_jaccard;
			sse_mash += cur_err * cur_err;
			ofs_mash << query_name << '\t' << ref_name << '\t' << standard_jaccard << '\t' << jaccard << '\t' << cur_err << endl;
		}
		ifs_mash.close();
		ofs_mash.close();
		cout << "the Sum of Squared Error (SSE) of mash_dist is: " << sse_mash << endl;
	}
	if(*option_bindash_file){
		cerr << "=====run the option bindash file " << endl;
		ifstream ifs_bindash;
		ifs_bindash.open(bindash_file);
		if(!ifs_bindash){
			cerr << "ERROR: cannot open bindash_file: " << bindash_file << endl;
			return 1;
		}
		string output_bindash = output_prefix + ".bindash.res";
		ofstream ofs_bindash(output_bindash);
		ofs_bindash << "query_name\tref_name\tstandard_jaccard\testimated_jaccard\terror" << endl;
		double sse_bindash = 0.0;
		while(getline(ifs_bindash, line)){
			stringstream ss;
			string query_name, ref_name, jaccard_str;
			double mutation_distance, p_value, jaccard;
			ss << line;
			ss >> query_name >> ref_name >> mutation_distance >> p_value >> jaccard_str;
			string merged_name = query_name + ref_name;
			if(name_jaccard_map.count(merged_name) == 0){
				cerr << "ERROR: the merged_name: " << merged_name << " is not int the standard name_jaccard_map" << endl;
				return 1;
			}
			double standard_jaccard = name_jaccard_map[merged_name];
			int mid_pos = jaccard_str.find_first_of('/');
			int inter_size = stoi(jaccard_str.substr(0, mid_pos));
			int union_size = stoi(jaccard_str.substr(mid_pos+1, jaccard_str.length()-mid_pos-1));
			jaccard = (double)inter_size / union_size;
			double cur_err = jaccard - standard_jaccard;
			sse_bindash += cur_err * cur_err;
			ofs_bindash << query_name << '\t' << ref_name << '\t' << standard_jaccard << '\t' << jaccard << '\t' << cur_err << endl;
		}
		ifs_bindash.close();
		ofs_bindash.close();
		cout << "the Sum of Squared Error (SSE) of bindash_dist is: " << sse_bindash << endl;
	}
	if(*option_dashing2_file){
		cerr << "=====run the option dashing2 file " << endl;
		string cmd_get_line = "wc -l " + dashing2_file;
		FILE * fp_cmd = popen(cmd_get_line.c_str(), "r");
		int buffer_len = 1024;
		char * buffer = new char[buffer_len];
		int read_len = fread(buffer, sizeof(char), buffer_len, fp_cmd);
		if(read_len == buffer_len){
			cerr << "the read length is equal buffer length, buffer size not enough" << endl;
			exit(1);
		}
		fclose(fp_cmd);
		string buffer_str(buffer);
		int dashing2_line_number = stoi(buffer_str);
		delete [] buffer;
		#ifdef DEBUG
		cerr << "the dashing2_line_number is: " << dashing2_line_number << endl;
		#endif
		int ref_number = dashing2_line_number - 3;
		ifstream ifs_dashing2;
		ifs_dashing2.open(dashing2_file);
		if(!ifs_dashing2){
			cerr << "ERROR: cannot open dashing2_file: " << dashing2_file << endl;
			return 1;
		}
		string output_dashing2 = output_prefix + ".dashing2.res";
		ofstream ofs_dashing2(output_dashing2);
		ofs_dashing2 << "query_name\tref_name\tstandard_jaccard\testimated_jaccard\terror" << endl;
		double sse_dashing2 = 0.0;
		getline(ifs_dashing2, line);//first line: application description
		getline(ifs_dashing2, line);//second line: application options
		getline(ifs_dashing2, line);//third line: ref_name and query_name list
		stringstream cur_ss;
		cur_ss << line;
		string source, cur_name;
		cur_ss >> source;
		
		int cur_index = 0;
		vector<string> ref_name_arr;
		vector<string> query_name_arr;
		while(cur_ss >> cur_name){
			if(cur_index < ref_number)
				ref_name_arr.push_back(cur_name);
			else
				query_name_arr.push_back(cur_name);
			cur_index++;
		}
		
		int line_index = 0;
		while(getline(ifs_dashing2, line)){
			stringstream ss;
			ss << line;
			string ref_name;
			ss >> ref_name;
			if(ref_name != ref_name_arr[line_index]){
				cerr << "ERROR: mismatched ref_name in dashing2 dist file" << endl;
				exit(1);
			}
			double jaccard;
			size_t cur_query_index = 0;
			while(ss >> jaccard){
				if(cur_query_index >= query_name_arr.size()){
					cerr << "ERROR: mismatched query_distance number and query_name number" << endl;
					exit(1);
				}
				string merged_name = query_name_arr[cur_query_index] + ref_name;
				if(name_jaccard_map.count(merged_name) == 0){
					cerr << "ERROR: the merged_name: " << merged_name << " is not int the standard name_jaccard_map" << endl;
					exit(1);
				}
				double standard_jaccard = name_jaccard_map[merged_name];
				double cur_err = jaccard - standard_jaccard;
				sse_dashing2 += cur_err * cur_err;
				ofs_dashing2 << query_name_arr[cur_query_index] << '\t' << ref_name << '\t' << standard_jaccard << '\t' << jaccard << '\t' << cur_err << endl;
				cur_query_index++;
			}
			line_index++;
			//double mutation_distance, p_value, jaccard;

		}
		ifs_dashing2.close();
		ofs_dashing2.close();
		cout << "the Sum of Squared Error (SSE) of dashing2_dist is: " << sse_dashing2 << endl;
	}

	return 0;
}
