//---------------------------------------------------------------------------------------------
// 5G NR LDPC encoding/decoding C++ encapsulation&implementation
// supports rate matching and rate recovery
// Dr J Mao 2021 Oct
// juquan.justin.mao@gmail.com
//--------------------------------------------------------------------------------------------

#ifndef NR_LDPC_H
#define NR_LDPC_H

#include <stdint.h> /*uint8_t, uint16_t*/
#include<numeric> /*size_t*/
#include<vector>
#include <assert.h>
#include <iostream>
#include <array>
#include <algorithm> /* sort, transform, rotate*/
#include <functional> /*xor*/

// 3GPP LDPC related tables
extern const uint16_t liftSizeTable[8][8];
extern const uint16_t shiftTableBgn_1[316][10];
extern const uint16_t shiftTableBgn_2[197][10];

//an edge is a connection from a check node to a variable node in a tanner graph
// along with the corresponding shift coefficient
struct edge_t {
	uint16_t cNodeIdx;
	uint16_t vNodeIdx;
	uint16_t nShifts;
};

// a layer includes all the connections(edges) to a check node in a tanner graph
struct layer_t {
	uint16_t edgeStart;
	uint16_t edgeEnd;
};

class nrLDPC {
public:
	std::size_t mKBar;          // info length excluding fillers
	std::size_t mK;             // info length including fillers
	double      mR;             // code rate
	uint8_t     mBGn;           // base graph identifier
	uint16_t    mZc;            // lifting size
	uint8_t     mSetIdx;        // the set that shift coefficients belongs to 0-7
	std::size_t mF;             // number of fillers to make info length is 22*Zc 0r 10*Zc
	std::size_t mN;             // mother code word length
	// edges are sorted by check node idx and variable node idx
	std::vector<edge_t> mEdges;
	// the i-th layer stores all the connections(edges) to the i-th check node in a tanner graph
	std::vector<layer_t> mLayers;

public:
	nrLDPC(std::size_t infoLen, double codeRate);
	// encoding
	std::vector<bool> encode(const std::vector<bool>& msg);
	// decoding
	std::vector<bool> decode(const std::vector<double>& llr, const unsigned nMaxIter);
	// rate matching
	std::vector<bool> rateMatch(const std::vector<bool>& bitsIn, std::size_t nOfBitOut);
	// rate recovery
	std::vector<double> rateRecover(const std::vector<double>& softBitsIn);
	// checkSum a code word
	bool checkSumCodeWord(std::vector<bool>& cw);
	// filler length
	size_t getFillerLength() { return mF; }

public:
	static uint8_t selectBaseGraph(std::size_t KBar, double R);
	static uint16_t selectLiftSize(std::size_t KBar, uint8_t BGn);
	static uint8_t selectShiftSet(uint16_t Zc);
	static std::vector<std::vector<bool>> makeParityCheckMatrix(uint8_t BGn, uint16_t Zc);

private:
	// core decoding module
	std::vector<std::vector<double>> checkNodeOperation(const std::vector<std::vector<double>>& msgIn);
	template<typename T>
	inline std::vector<std::size_t> sort_indexes(const std::vector<T>& v);
	std::vector<std::vector<double>> transposeMat(const std::vector<std::vector<double>>& mat);

	// wrapper functions to enhance readability
	template<typename T>
	inline std::vector<T> circShift(const std::vector<T>& vecIn, const unsigned nShifts);
	inline std::vector<bool> eleWiseXor(const std::vector<bool>& vec1, const std::vector<bool>& vec2);
	inline std::vector<double> eleWisePlus(const std::vector<double>& vec1, const std::vector<double>& vec2);
	inline std::vector<double> eleWiseMinus(const std::vector<double>& vec1, const std::vector<double>& vec2);
};

#endif
