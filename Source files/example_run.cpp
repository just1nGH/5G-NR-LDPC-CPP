#include "LDPC.h"
#include <iostream>
#include <random>
#include <chrono> //timing
#include <cstdio> // sprintf
using namespace std;

using chrono::milliseconds;
using std::chrono::duration_cast;
using std::chrono::steady_clock;

vector<double> lin_space(double start, double ed, int num);
void example_run_LDPC();
using namespace std;

int main()
{
	std::cout << "Hello, 5G-NR!" << std::endl;

	example_run_LDPC();

	return 0;
}

void example_run_LDPC()
{
	unsigned nMaxIter = 8; // number of decoders
	int M = 1000;		   // code word length
	double rate = 1.0 / 3; // code rate

	unsigned K = unsigned(ceil(M * rate)); // information length

	// instantiates a POLAR object
	nrLDPC ldpc = nrLDPC(K, rate);

	// random engines
	default_random_engine random_engine;
	bernoulli_distribution bern_dist;
	normal_distribution<double> norm_dist(0, 1);

	// Running parameters
	vector<double> EsN0_dB = lin_space(-4, -2, 9);
	vector<double> N0(EsN0_dB.size(), 0);
	transform(EsN0_dB.begin(), EsN0_dB.end(), N0.begin(), [](const float &x)
			  { return pow(10.0, -x / 10.0); });

	vector<double> ber(N0.size(), 0), bler(N0.size(), 0);
	vector<unsigned> n_bit_errs(N0.size(), 0), n_blk_errs(N0.size(), 0);

	unsigned n_max_blks = 10000;
	vector<bool> fillers(ldpc.getFillerLength(), 0);
	// loop each SNR
	for (unsigned i = 0; i < N0.size(); i++)
	{
		// print progress
		char str[100];
		sprintf(str, "\nNow running EsN0: %.2f dB [%d of %lu]", EsN0_dB[i], i + 1, N0.size());
		cout << str << endl;
		unsigned print_len = 0;

		unsigned n_blks_done = 0;
		clock_t tStart = clock(); // timing

		while ((n_blks_done < n_max_blks) && (n_blk_errs[i] < 100))
		{
			// generate random bit stream
			vector<bool> msg; // = { 1,0,0,1,1,0,1,0,0,0,0,0,1,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,1,0,1,0,1 };
			msg.reserve(K);
			for (unsigned j = 0; j < K; j++)
				msg.push_back(bern_dist(random_engine));

			// add filler bits
			vector<bool> extMsg = msg;
			extMsg.insert(extMsg.end(), fillers.begin(), fillers.end());

			// LDPC encoding
			vector<bool> enc = ldpc.encode(extMsg);
			assert(ldpc.checkSumCodeWord(enc));

			// rate matching
			vector<bool> rm_enc = ldpc.rateMatch(enc, M);

			// BPSK + AWGN
			vector<double> r;
			r.reserve(M);
			for (auto e : rm_enc)
				r.push_back(1 - 2.0 * e + sqrt(N0[i] / 2.0) * norm_dist(random_engine));

			// compute soft bits as LLR
			vector<double> llr;
			llr.reserve(M);
			for (auto e : r)
				llr.push_back(4.0 * e / N0[i]);

			// rate recovery
			vector<double> rr_llr = ldpc.rateRecover(llr);

			// scl decoding
			vector<bool> msg_cap = ldpc.decode(rr_llr, nMaxIter);

			// count errors
			unsigned n_errs = 0;
			for (unsigned j = 0; j < K; j++)
			{
				if (msg[j] != msg_cap[j])
					n_errs++;
			}

			if (n_errs)
			{
				n_bit_errs[i] += n_errs;
				n_blk_errs[i]++;
			}

			n_blks_done += 1;

			ber[i] = n_bit_errs[i] * 1.0 / K / n_blks_done;
			bler[i] = n_blk_errs[i] * 1.0 / n_blks_done;

			// print progress for every 10 blocks
			if (n_blks_done % 10 == 0 || n_blks_done == 1)
			{
				print_len = sprintf(str, "Elapsed time: %.1f seconds, # tx blocks: %d,# error blocks:%d, ber: %.5f, bler %.5f", double(clock() - tStart) / CLOCKS_PER_SEC, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
				cout << std::string(print_len, '\b') << str << flush;
			}
		}

		// print  progress when one SNR is finished
		sprintf(str, "Elapsed time: %.1f seconds, # tx blocks: %d,# error blocks:%d, ber: %.5f, bler %.5f", static_cast<double>(clock() - tStart) / CLOCKS_PER_SEC, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
		cout << std::string(print_len, '\b') << str << flush;
	}

	// print simulation result
	cout << endl;
	cout << "Modulation:" << "BPSK" << endl;
	cout << "[M,R] = [ " << M << "," << rate << "]" << endl;
	cout << "EsN0_dB = [";
	for (auto e : EsN0_dB)
		cout << e << " ";
	cout << "]" << endl;

	cout << "BER = [";
	for (auto e : ber)
		cout << e << " ";
	cout << "]" << endl;

	cout << "BLER = [";
	for (auto e : bler)
		cout << e << " ";
	cout << "]" << endl;
}

vector<double> lin_space(double start, double end, int num)
{
	// catch rarely, throw often
	assert(num >= 2 && "The third parameter must be a positive integer >= 2!");

	int partitions = num - 1;
	vector<double> pts;
	// length of each segment
	double length = (end - start) / partitions;
	// first, not to change
	pts.push_back(start);
	for (int i = 1; i < num - 1; i++)
	{
		pts.push_back(start + i * length);
	}
	// last, not to change
	pts.push_back(end);
	return pts;
}