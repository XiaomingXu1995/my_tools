/* Author: Xiaoming Xu
 * Data: 2022/5/12
 * 
 * See the LICENSE.txt file included with this software for licence information.
 */
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


using namespace std;


int main(int argc, char* argv[]){
	CLI::App app{"simulate_mutation_sequence: simulate long genomic DNA sequences with specific genome length and mutation rate"};
	int error_rate = 50;
	int seq_length = 1000000;
	int num_seed_seq = 1;
	int num_each_clust = 10;
	bool separate_file = false;
	string output_prefix = "output";
	int rand_seed = 23;
	auto option_error_rate = app.add_option("-e, --error-rate", error_rate, "set the error_rate* 1000 ");
	auto option_seq_length = app.add_option("-l, --sequence-length", seq_length, "set the length of each sequence");
	auto option_num_seed_seq = app.add_option("-s, --num-seed", num_seed_seq, "set the number of seed sequences");
	auto option_num_each_clust = app.add_option("-c, --num-each-clust", num_each_clust, "set the number of mutated sequence for each cluster");
	auto option_output_prefix = app.add_option("-o, --output", output_prefix, "set the output prefix name of output file (output directory for separate file flag)");
	app.add_flag("--separate-file", separate_file, "each sequence serves as a single file");
	auto option_rand_seed = app.add_option("--rand-seed", rand_seed, "set the seed for random()");

	option_output_prefix->required();
	option_num_each_clust->required();
	option_num_seed_seq->required();
	option_seq_length->required();
	option_error_rate->required();

	CLI11_PARSE(app, argc, argv);

  const char nucs[4] = { 'A','T','G','C'};

	string out_seed_file = output_prefix + "_seed.fna";
	string out_total_file = output_prefix + "_total.fna";
	string out_ground_truth = output_prefix + "_groundTruth";
	string out_dir;
	FILE * fp_ground_truth = fopen(out_ground_truth.c_str(), "w");
	FILE * fp_seed_file;
	FILE * fp_total_file;
	if(separate_file){
		out_dir = output_prefix + "_dir";
		string cmd0 = "mkdir -p " + out_dir;
		int finished = system(cmd0.c_str());
		if(!finished){
			cerr << "finished create the directory: " << out_dir << endl;
		}
	}
	else{
		fp_seed_file = fopen(out_seed_file.c_str(), "w");
		fp_total_file = fopen(out_total_file.c_str(), "w");
	}

	cerr << "the error rate is: " << double(error_rate)/1000 << endl;
	cerr << "the number of clusters is: " << num_seed_seq << endl;
	cerr << "the number of sequences in each cluster is: " << num_each_clust << endl;
	cerr << "the approximate sequence length is: " << seq_length << endl;
	cerr << "the output seed sequences file is: " << out_seed_file << endl;
	cerr << "the output total sequences file is: " << out_total_file << endl;
	cerr << "the groundTruth file is: " << out_ground_truth << endl;


	string key1 = "seqName";
	string key2 = "taxid";
	fprintf(fp_ground_truth, "%s\t%s\n",key1.c_str(), key2.c_str());
	
	for(int i = 0; i < num_seed_seq; i++)
	{
		FILE * cur_fp_seed_file;
		int prefix_length = out_seed_file.find_last_of(".");
		string cur_seed_file_name = out_dir + '/' + out_seed_file.substr(0, prefix_length) + '_' + to_string(i) + ".fna";
		if(separate_file){
			cur_fp_seed_file = fopen(cur_seed_file_name.c_str(), "w");
		}
		struct timeval tv;
		gettimeofday(&tv, NULL);
		if(*option_rand_seed){
			srand(rand_seed);
		}
		else
  		srand(tv.tv_usec);
		string seqName = ">seq_" + to_string(i);
		string groundTruthName = seqName.substr(1);
		if(!separate_file)
			fprintf(fp_ground_truth, "%s\t%d\n", groundTruthName.c_str(), i);
		else
			fprintf(fp_ground_truth, "%s\t%d\n", cur_seed_file_name.c_str(), i);
		string seqComment = "Seed sequence " + to_string(i) + " to generate mutations";
		string infoLine = seqName + '\t' + seqComment + '\n';
		string seedSeq("");
		for(int i = 0; i < seq_length; i++)
		{
			char newC = nucs[random()%4];
			seedSeq += newC;
		}
		//output the seed sequence into outSeedFile
		int infoLineLen = infoLine.length();
		if(!separate_file){
			fwrite(infoLine.c_str(), sizeof(char), infoLineLen, fp_seed_file);
			fwrite(infoLine.c_str(), sizeof(char), infoLineLen, fp_total_file);
		}
		else{
			fwrite(infoLine.c_str(), sizeof(char), infoLineLen, cur_fp_seed_file);
		}
		int seedLen = seedSeq.length();
		string outSeedSeq("");
		for(int k = 0; k < seedLen; k += 80)
		{
			int curLength = std::min(80, seedLen-k);
			string tmpLine = seedSeq.substr(k, curLength);
			outSeedSeq += tmpLine + '\n';
		}
		int outSeedSeqLen = outSeedSeq.length();
		if(!separate_file){
			fwrite(outSeedSeq.c_str(), sizeof(char), outSeedSeqLen, fp_seed_file);
			fwrite(outSeedSeq.c_str(), sizeof(char), outSeedSeqLen, fp_total_file);
		}
		else{
			fwrite(outSeedSeq.c_str(), sizeof(char), outSeedSeqLen, cur_fp_seed_file);
			fclose(cur_fp_seed_file);
		}


		//for generate mutation sequences
		for(int j = 0; j < num_each_clust; j++)
		{
			string mutationName = seqName + "_mutation_" + to_string(j);
			string cur_mutation_file_name;
			FILE * cur_fp_mutation_file;
			if(separate_file){
				int cur_prefix_length = cur_seed_file_name.find_last_of('.');
				cur_mutation_file_name = cur_seed_file_name.substr(0, cur_prefix_length) + "_mutation_" + to_string(j) + ".fna";
				//cur_mutation_file_name = cur_seed_file_name + "_mutation_" + to_string(j);
				cur_fp_mutation_file = fopen(cur_mutation_file_name.c_str(), "w");
			}
			string groundTruthMuName = mutationName.substr(1);
			if(separate_file)
				fprintf(fp_ground_truth, "%s\t%d\n", cur_mutation_file_name.c_str(), i);
			else
				fprintf(fp_ground_truth, "%s\t%d\n", groundTruthMuName.c_str(), i);
			string mutationComment = "mutation sequence " + to_string(j) + " from seedSequence " + to_string(i);
			string mutaInfoLine = mutationName + '\t' + mutationComment + '\n';
			string mutationSeq("");
			for(size_t t = 0; t < seedSeq.length(); t++)
			{
				if(random()%1000 < error_rate){
					int mut = random()%3;
					if(mut == 0)//sub
					{
						while(1){
							char newc = nucs[random()%4];
							if(newc != seedSeq[t]){
								mutationSeq += newc;
								break;
							}
						}
					}
					else if(mut == 1)// ins
					{
						mutationSeq += nucs[random()%4];
						t = t - 1;
					}
					else//del
						continue;
				}//end if mutation
				else// no mutation
					mutationSeq += seedSeq[t];
			}
			int mutaInfoLineLen = mutaInfoLine.length();
			if(!separate_file)
				fwrite(mutaInfoLine.c_str(), sizeof(char), mutaInfoLineLen, fp_total_file);
			else
				fwrite(mutaInfoLine.c_str(), sizeof(char), mutaInfoLineLen, cur_fp_mutation_file);
			int mutationLen = mutationSeq.length();
			string outMutationSeq("");
			for(int k = 0; k < mutationLen; k += 80)
			{
				int curLength = std::min(80, mutationLen-k);
				string tmpLine = mutationSeq.substr(k, curLength);
				outMutationSeq += tmpLine + '\n';
			}
			int outMutationSeqLen = outMutationSeq.length();
			if(!separate_file)
				fwrite(outMutationSeq.c_str(), sizeof(char), outMutationSeqLen, fp_total_file);
			else{
				fwrite(outMutationSeq.c_str(), sizeof(char), outMutationSeqLen, cur_fp_mutation_file);
				fclose(cur_fp_mutation_file);
			}
		}
	}
	if(!separate_file){
		fclose(fp_seed_file);
		fclose(fp_total_file);
	}
	fclose(fp_ground_truth);

	cerr << "finish generate mutation files with multithread " << endl;

  return 0;
}

