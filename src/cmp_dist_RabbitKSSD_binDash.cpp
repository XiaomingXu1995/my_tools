#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <cassert>
#include <sys/time.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <assert.h>
#include "CLI11.hpp"
#include "robin_hood.h"
#include <cstdio>
#include <sys/sysinfo.h>
#include <omp.h>
#include <queue>
#include "common.hpp"
#include <unistd.h>


using namespace std;

//void get_rabbit_matrix_from_dir(double * rabbit_dist_matrix, string dir_path, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads);
//void get_bindash_matrix(double * bindash_dist_matrix, string bindash_dist_file, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads);
void get_rabbit_matrix_from_dir(float* rabbit_dist_matrix, string dir_path, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads);
void get_bindash_matrix(float* bindash_dist_matrix, string bindash_dist_file, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads);


int main(int argc, char* argv[]){
	CLI::App app{"cmp_dist_RabbitKSSD_binDash: compute the relationship between RabbitKSSD and binDash"};
	string rabbit_dist_file = "bacteria_L3K10.alldist";
	string bindash_dist_file = "bacteria.dist";
	double dist_threshold = 0.05;
	string genome_list_file = "bacteria.list";
	bool is_rabbit_dist_list = false;
	string output_file = "result.out";
	int threads = get_nprocs_conf();
	
	auto option_rabbit_file = app.add_option("--rabbit-file", rabbit_dist_file, "set the distance file of RabbitKSSD");
	auto option_bindash_file = app.add_option("--bindash-file", bindash_dist_file, "set the distance file of binDash");
	auto option_dist_threshold = app.add_option("-D, --threshold", dist_threshold, "set the distance threshold of RabbitKSSD and binDash");
	auto option_genome_list_file = app.add_option("--genome-list", genome_list_file, "set the genome list of RabbitKSSD and binDash, which has been used for sketching and distance computation");
	auto option_output_result = app.add_option("-o, --output", output_file, "the output of distance comparision between RabbitKSSD and bindash");
	app.add_option("-p, --threads", threads, "set the thread number");
	app.add_flag("-l, --rabbit-dist-list", is_rabbit_dist_list, "the input distance file of RabbitKSSD is file list in the directory");

	option_rabbit_file->required();
	option_bindash_file->required();
	option_dist_threshold->required();
	option_genome_list_file->required();
	option_output_result->required();

	//flag_dist_list->required();

	CLI11_PARSE(app, argc, argv);

	double t1 = get_sec();
	ifstream ifs_genome_list, ifs_bindash;
	ifs_genome_list.open(genome_list_file);
	string line;
	robin_hood::unordered_map<string, int> genome_id_map;
	size_t genome_id = 0;
	while(getline(ifs_genome_list, line)){
		genome_id_map.insert({line, genome_id++});
	}
	ifs_genome_list.close();
	cerr << "the number of genome is: " << genome_id_map.size() << endl;
	cerr << "the final genome_id is: " << genome_id << endl;
	size_t genome_number = genome_id_map.size();
	double t2 = get_sec();
	cerr << "=====time of generate genome_id_map is: " << t2 - t1 << endl;

	ifs_bindash.open(bindash_dist_file);
	//double * bindash_dist_matrix = new double[genome_number * genome_number];
	//double * rabbit_dist_matrix = new double[genome_number * genome_number];
	float * bindash_dist_matrix = new float[genome_number * genome_number];
	float * rabbit_dist_matrix = new float[genome_number * genome_number];
	#pragma omp parallel for num_threads(threads)
	for(size_t i = 0; i < genome_number * genome_number; i++){
		bindash_dist_matrix[i] = -1.0;
		rabbit_dist_matrix[i] = -1.0;
	}
	double t3 = get_sec();
	cerr << "=====time of init distance matrix is: " << t3 - t2 << endl;

	get_bindash_matrix(bindash_dist_matrix, bindash_dist_file, genome_id_map, genome_number, threads);

	double t4 = get_sec();
	cerr << "=====time of multithreading update bindash_dist_matrix is: " << t4 - t3 << endl;

	size_t rabbit_line = 0;
	if(is_rabbit_dist_list){
		get_rabbit_matrix_from_dir(rabbit_dist_matrix, rabbit_dist_file, genome_id_map, genome_number, threads);
	}
	else{
		ifstream ifs_dist;
		ifs_dist.open(rabbit_dist_file);
		string cur_line;
		int line_index = 0;
		while(getline(ifs_dist, cur_line)){
			if(line_index == 0){
				line_index++;
				continue;
			}
			stringstream cur_ss;
			string ref_name, query_name, common_union_str;
			double jaccard, dist;
			cur_ss << cur_line;
			cur_ss >> query_name >> ref_name >> common_union_str >> jaccard >> dist;
			if(genome_id_map.count(query_name) == 0 || genome_id_map.count(ref_name) == 0){
				cerr << query_name << '\t' << ref_name << endl;
				exit(1);
			}
			int query_id = genome_id_map[query_name];
			int ref_id = genome_id_map[ref_name];
			rabbit_dist_matrix[ref_id * genome_number + query_id] = dist;
			line_index++;
			rabbit_line++;
		}
		ifs_dist.close();
	}
	double t5 = get_sec();
	cerr << "=====time of updating rabbit_dist_matrix is: " << t5 - t4 << endl;
	cerr << "finished the rabbit dist matrix" << endl;
	cerr << "the line number of rabbit_dist_file: " << rabbit_dist_file << " is: " << rabbit_line << endl;

	FILE * fp_out = fopen(output_file.c_str(), "w+");
	if(!fp_out){
		cerr << "cannot open the output file: " << output_file << endl;
		exit(1);
	}
	//fprintf(fp_out, "query_id\tref_id\tRabbitKSSD_dist\tbindash_dist\terror\terror_rate*100\n");
	FILE ** fp_arr = new FILE*[threads];
	for(int i = 0; i < threads; i++){
		string tmp_file = "tmp." + to_string(i);
		fp_arr[i] = fopen(tmp_file.c_str(), "w");
	}
	#pragma omp parallel for num_threads(threads) schedule(static)
	for(size_t i = 0; i < genome_number * genome_number; i++){
		int tid = omp_get_thread_num();
		if(rabbit_dist_matrix[i] != -1 && bindash_dist_matrix[i] != -1){
			int ref_id = i / genome_number;
			int query_id = i % genome_number;
			double error, error_rate;
			double max_one = std::max(rabbit_dist_matrix[i], bindash_dist_matrix[i]);
			double min_one = std::min(rabbit_dist_matrix[i], bindash_dist_matrix[i]);
			error = max_one - min_one;
			if(min_one == 0.0)
				error_rate = 0.0;
			else
				error_rate = error / min_one * 100;
			//fprintf(fp_out, "%d\t%d\t%lf\t%lf\t%lf\t%lf\n", query_id, ref_id, rabbit_dist_matrix[i], bindash_dist_matrix[i], error, error_rate);
			fprintf(fp_arr[tid], "%d\t%d\t%lf\t%lf\t%lf\t%lf\n", query_id, ref_id, rabbit_dist_matrix[i], bindash_dist_matrix[i], error, error_rate);
		}
	}
	size_t buffer_size = 1LLU << 24;
	char * buffer = new char[buffer_size];
	for(int i = 0; i < threads; i++){
		fclose(fp_arr[i]);
		string tmp_file = "tmp." + to_string(i);
		FILE * cur_fp_in = fopen(tmp_file.c_str(), "r");
		if(!cur_fp_in){
			cerr << "error open the file: " << tmp_file << endl;
			exit(1);
		}
		while(1){
			size_t read_length = fread(buffer, sizeof(char), buffer_size, cur_fp_in);
			fwrite(buffer, sizeof(char), read_length, fp_out);
			if(read_length < buffer_size) break;
		}
		fclose(cur_fp_in);
		remove(tmp_file.c_str());
	}
	fclose(fp_out);
	delete [] buffer;
	delete [] fp_arr;
	double t6 = get_sec();
	cerr << "=====time of comparing the two matrix is: " << t6 - t5 << endl;
	delete [] rabbit_dist_matrix;
	delete [] bindash_dist_matrix;
	cerr << "finished" << endl;
	return 0;
}

//void get_rabbit_matrix_from_dir(double * rabbit_dist_matrix, string dir_path, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads){
void get_rabbit_matrix_from_dir(float* rabbit_dist_matrix, string dir_path, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads){
	string cmd_get_dist_list = "ls " + dir_path;
	FILE * fp_popen = popen(cmd_get_dist_list.c_str(), "r");
	if(!fp_popen){
		cerr << "error popen to execute: " << cmd_get_dist_list << endl;
		exit(1);
	}
	size_t buffer_size = 4 * 1024 * 1024;
	char * buffer = new char [buffer_size];
	string str_dist_file_list("");
	while(1){
		size_t read_length = fread(buffer, sizeof(char), buffer_size, fp_popen);
		str_dist_file_list.append(buffer, read_length);
		if(read_length < buffer_size) break;
	}
	pclose(fp_popen);
	delete [] buffer;
	stringstream ss;
	ss << str_dist_file_list;
	vector<string> dist_file_arr;
	string cur_dist_file;
	while(ss >> cur_dist_file){
		cur_dist_file = dir_path + '/' + cur_dist_file;
		dist_file_arr.push_back(cur_dist_file);
	}
	cerr << "the rabbit dist file number is: " << dist_file_arr.size() << endl;
	#pragma omp parallel for num_threads(threads)
	for(size_t i = 0; i < dist_file_arr.size(); i++){
		FILE * cur_fp = fopen(dist_file_arr[i].c_str(), "r");
		if(!cur_fp){
			cerr << "cannot open the file: " << dist_file_arr[i] << endl;
			exit(1);
		}
		size_t buf_length = 1LLU << 24;
		char * cur_buffer = new char[buf_length];
		size_t cur_buffer_length = buf_length;
		while(1){
			size_t read_length = fread(cur_buffer, sizeof(char), buf_length, cur_fp);
			cur_buffer_length = read_length;
			if(read_length == buf_length){
				size_t start_index = read_length - 1000;//a line length must be less than 1000
				size_t offset = 0;
				for(size_t j = start_index; j < buf_length; j++){
					if(cur_buffer[j] == '\n'){
						offset = -1 * (buf_length - 1 - j);
						fseek(cur_fp, offset, SEEK_CUR);
						break;
					}
				}
				cur_buffer_length = buf_length + offset;
			}
			string cur_str("");
			int element_index = 0;
			int num_element = 5;
			int query_id = 0;
			int ref_id = 0;
			double dist = 1.0;
			for(size_t j = 0; j < cur_buffer_length; j++){
				if(cur_buffer[j] == '\t'){
					if(element_index % num_element == 0){
						if(genome_id_map.count(cur_str) == 0){
							cerr << "query_name is not in the genome_id_map: " << cur_str << endl;
							exit(1);
						}
						query_id = genome_id_map[cur_str];
					}
					else if(element_index % num_element == 1){
						if(genome_id_map.count(cur_str) == 0){
							cerr << "ref_name is not in the genome_id_map: " << cur_str << endl;
							exit(1);
						}
						ref_id = genome_id_map[cur_str];
					}
					cur_str = "";
					element_index++;
				}
				else if(cur_buffer[j] == '\n'){
					dist = stof(cur_str);
					rabbit_dist_matrix[ref_id * genome_number + query_id] = dist;
					cur_str = "";
					element_index++;
				}
				else{
					cur_str +=cur_buffer[j];
				}
			}

			if(read_length < buf_length) break;
		}
		delete [] cur_buffer;
		cerr << "finish the " << dist_file_arr[i] << endl;
		fclose(cur_fp);
	}

}

//void get_bindash_matrix(double * bindash_dist_matrix, string bindash_dist_file, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads){
void get_bindash_matrix(float* bindash_dist_matrix, string bindash_dist_file, robin_hood::unordered_map<string, int>& genome_id_map, size_t genome_number, int threads){
	int buffer_number = 128;
	size_t buffer_size = 1 << 24; //16MB
	char** buffers = new char*[buffer_number];
	queue<char*> empty_queue;
	queue<pair<char*, size_t>> full_queue;
	for(int i = 0; i < buffer_number; i++){
		buffers[i] = new char[buffer_size];
		empty_queue.push(buffers[i]);
	}

	omp_lock_t empty_lock;
	omp_lock_t full_lock;
	omp_init_lock(&empty_lock);
	omp_init_lock(&full_lock);

	bool finished_parse = false;
	int parsed_number = 0;
	int solved_number = 0;

	FILE * fp_in = fopen(bindash_dist_file.c_str(), "r");
	if(!fp_in){
		cerr << "ERROR: get_bindash_matrix(), cannot open the distance file: " << bindash_dist_file << endl;
		exit(1);
	}
	
	#pragma omp parallel num_threads(threads)
	{
		int tid = omp_get_thread_num();
		if(tid == 0)//producer
		{
			while(1){
				omp_set_lock(&empty_lock);
				if(empty_queue.empty()){
					omp_unset_lock(&empty_lock);
					continue;
				}
				char * cur_buffer = empty_queue.front();
				empty_queue.pop();
				omp_unset_lock(&empty_lock);
				size_t read_length = fread(cur_buffer, sizeof(char), buffer_size, fp_in);
				size_t cur_buffer_length = read_length;
				if(read_length == buffer_size){
					size_t start_index = read_length - 1000;//a line length must be less than 1000
					size_t offset = 0;
					for(size_t i = start_index; i < buffer_size; i++){
						if(cur_buffer[i] == '\n'){
							offset = -1 * (buffer_size - 1 - i);
							fseek(fp_in, offset, SEEK_CUR);
							break;
						}
					}
					cur_buffer_length = buffer_size + offset;
				}

				omp_set_lock(&full_lock);
				full_queue.push({cur_buffer, cur_buffer_length});
				parsed_number++;
				//if(parsed_number % 20 == 0) cerr << "the parsed_number is: " << parsed_number << endl;
				omp_unset_lock(&full_lock);
				if(read_length < buffer_size){
					finished_parse = true;
					break;
				}
			}
		}
		else//consumer
		{
			while(1){
				omp_set_lock(&full_lock);
				if(full_queue.empty()){
					omp_unset_lock(&full_lock);
					if(finished_parse && solved_number >= parsed_number) break;
					continue;
				}
				pair<char*, size_t> cur_pair = full_queue.front();
				full_queue.pop();
				omp_unset_lock(&full_lock);

				char* cur_buffer = cur_pair.first;
				size_t cur_length = cur_pair.second;
				string cur_line("");
				string cur_str("");
				int element_index = 0;
				int num_element = 5;
				int query_id = 0;
				int ref_id = 0;
				double dist = 1.0;
				for(size_t i = 0; i < cur_length; i++){
					if(cur_buffer[i] == '\t'){
						if(element_index % num_element == 0){
							if(genome_id_map.count(cur_str) == 0){
								cerr << "query_name is not in the genome_id_map: " << cur_str << endl;
								exit(1);
							}
							query_id = genome_id_map[cur_str];
						}
						else if(element_index % num_element == 1){
							if(genome_id_map.count(cur_str) == 0){
								cerr << "ref_name is not in the genome_id_map: " << cur_str << endl;
								exit(1);
							}
							ref_id = genome_id_map[cur_str];
						}
						else if(element_index % num_element == 2){
							dist = stof(cur_str);
						}
						cur_str = "";
						element_index++;
					}
					else if(cur_buffer[i] == '\n'){
						bindash_dist_matrix[ref_id * genome_number + query_id] = dist;
						cur_str = "";
						element_index++;
					}
					else{
						cur_str +=cur_buffer[i];
					}
				}
				if(cur_line != ""){
					cerr << "the cur_buffer is not end with n, exit" << endl;
					exit(1);
				}
				omp_set_lock(&empty_lock);
				solved_number++;
				empty_queue.push(cur_buffer);
				omp_unset_lock(&empty_lock);
			}
		}
	}
	fclose(fp_in);
	for(int i = 0; i < buffer_number; i++){
		delete [] buffers[i];
	}
	delete [] buffers;

}












