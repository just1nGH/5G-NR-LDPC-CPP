

#include "LDPC.h"
using namespace  std;

nrLDPC::nrLDPC(size_t infoLen, double codeRate )
{
    mKBar = infoLen;
    mR = codeRate;

    // select base graph based on 3GPP 38.212 7.2.2
    mBGn = selectBaseGraph(mKBar,mR);

    // select lifting size
    mZc = selectLiftSize(mKBar,mBGn);

    // select shifting set
    mSetIdx = selectShiftSet(mZc);

    // systematic bits length and parity bits length
    if(mBGn ==1){
        mK = 22*mZc; mN = 68*mZc;
    }else{
        mK = 10*mZc;  mN = 52*mZc;
    }

    // fillers length
    mF = mK- mKBar;

    // build up edges and shifts
    if (mBGn == 1){
        mEdges.reserve(316); mShifts.reserve(316);
        for (unsigned i = 0; i < 316; i++){
            mEdges.push_back(edge_t());
            mEdges[i] = {shiftTableBgn_1[i][0],shiftTableBgn_1[i][1]};
            mShifts.push_back(shiftTableBgn_1[i][mSetIdx+2]);
        }
    }
    else{
        mEdges.reserve(197); mShifts.reserve(197);
        for (unsigned i = 0; i < 197; i++){
            mEdges.push_back(edge_t());
            mEdges[i] = {shiftTableBgn_2[i][0],shiftTableBgn_2[i][1]};
            mShifts.push_back(shiftTableBgn_2[i][mSetIdx+2]);
        }
    }

    // build up layers
    if(mBGn == 1){
        mLayers = {{0,19},{19,38},{38,57},{57,76},{76,79},{79,87},{87,96},{96,103},
                   {103,113},{113,122},{122,129},{129,137},{137,144},{144,150},{150,157},
                   {157,164},{164,170},{170,176},{176,182},{182,188},{188,194},{194,200},
                   {200,205},{205,210},{210,216},{216,221},{221,226},{226,230},{230,235},
                   {235,240},{240,245},{245,250},{250,255},{255,260},{260,265},{265,270},
                   {270,275},{275,279},{279,284},{284,289},{289,293},{293,298},{298,302},
                   {302,307},{307,312},{312,316}};
    }
    else{
        mLayers = {{0,8},{8,18},{18,26},{26,36},{36,40},{40,46},{46,52},
                   {52,58},{58,62},{62,67},{67,72},{72,77},{77,81},{81,86},
                   {86,91},{91,95},{95,100},{100,105},{105,109},{109,113},
                   {113,117},{117,121},{121,124},{124,128},{128,132},{132,135},
                   {135,140},{140,143},{143,147},{147,150},{150,155},{155,158},
                   {158,162},{162,166},{166,170},{170,174},{174,178},{178,181},
                   {181,185},{185,189},{189,193},{193,197}};
    }
}

vector<bool> nrLDPC::encode(const vector<bool>& msg)
{
    size_t  Kb, Cb, totLayers;
    if (mBGn ==1){
        Kb = 22; Cb = 68; totLayers = 46;
    }
    else {
        Kb = 10; Cb = 52; totLayers = 42;
    }

    assert(Kb*mZc == msg.size());

    // initialize encoded bits in nodes(vectors of size Zc)
    // the first Kb nodes correspond to information bits, the rest  parity bits
    vector<vector<bool>> cWord(Cb);
    for (unsigned i = 0; i < Kb; i++){
        cWord[i] = vector<bool>(msg.begin()+i*mZc,msg.begin()+(i+1)*mZc);
    }
    for (unsigned i = Kb; i < Cb; i++){
        cWord[i] = vector<bool>(mZc,0);
    }


    uint16_t vNodeIdx, nShifts, shiftP0;
    vector<bool> tmpWord(mZc,0);

    // solve the first parity node P0
    for (unsigned i = 0; i < 4; i++){
        for(unsigned edgeIdx = mLayers[i].edgeStart; edgeIdx < mLayers[i].edgeEnd; edgeIdx ++ ){
            vNodeIdx = mEdges[edgeIdx].vNodeIdx;
            nShifts = mShifts[edgeIdx];
            // only count information bit nodes
            if(vNodeIdx < Kb){
                tmpWord = cWord[vNodeIdx];
                rotate(tmpWord.rbegin(),tmpWord.rbegin()+ nShifts,tmpWord.rend());
                transform(tmpWord.begin(),tmpWord.end(),cWord[Kb].begin(),cWord[Kb].begin(),bit_xor<>{});
            }
            // find the shift coefficient of P0
            if (vNodeIdx == Kb  && (i == 1 || i == 2)){
                shiftP0 = nShifts;
            }
        }
    }
    // rotate back to get P0
    rotate(cWord[Kb].rbegin(),cWord[Kb].rbegin()+ (mZc-shiftP0),cWord[Kb].rend());

    // solve P1,P2,P3
    for (unsigned i = 0; i < 3; i++){
        for(unsigned edgeIdx = mLayers[i].edgeStart; edgeIdx <  mLayers[i].edgeEnd; edgeIdx ++ ) {
            vNodeIdx = mEdges[edgeIdx].vNodeIdx;
            nShifts = mShifts[edgeIdx];
            if(vNodeIdx < Kb+i){
                tmpWord = cWord[vNodeIdx];
                rotate(tmpWord.rbegin(),tmpWord.rbegin()+ nShifts,tmpWord.rend());
                transform(tmpWord.begin(),tmpWord.end(),cWord[Kb+i+1].begin(),cWord[Kb+i+1].begin(), bit_xor<>{});
            }
        }
    }
    // solve the rest parity node
    for (unsigned i = 4; i < totLayers; i++){
        // not taking the last edge which corresponds the parity position
        for(unsigned edgeIdx = mLayers[i].edgeStart; edgeIdx <  mLayers[i].edgeEnd-1; edgeIdx ++ ) {
            vNodeIdx = mEdges[edgeIdx].vNodeIdx;
            nShifts = mShifts[edgeIdx];
            if(vNodeIdx < Kb+i){
                tmpWord = cWord[vNodeIdx];
                rotate(tmpWord.rbegin(),tmpWord.rbegin()+ nShifts,tmpWord.rend());
                transform(tmpWord.begin(),tmpWord.end(),cWord[Kb+i].begin(), cWord[Kb+i].begin(),bit_xor<>{});
            }
        }
    }
    // flatten codeword and return
    vector<bool> flattenCWord;
    for (auto vec:cWord){
        flattenCWord.insert(flattenCWord.end(), vec.begin(),vec.end());
    }
    return flattenCWord;
}
vector<bool> nrLDPC::decode(const vector<double> &softBitsIn, const unsigned nMaxIter)
{
    //------------------------------------------------------------------------------------------------------
    // [ref] Hocevar, D.E. "A reduced complexity decoder architecture via layered decoding of LDPC codes."
    // In IEEE Workshop on Signal Processing Systems, 2004. SIPS 2004.
    //------------------------------------------------------------------------------------------------------
    assert(softBitsIn.size() == mN);

    // initialize LLR in blocks(nodes), each node with Zc bits
    vector<vector<double>> LLR(mN/mZc);
    for(unsigned i = 0; i < mN/mZc; i++ ){
        LLR[i] = vector<double>(softBitsIn.begin()+i*mZc,softBitsIn.begin()+(i+1)*mZc);
    }

    // find how many parity nodes to use for decoding
    unsigned nMaxLayer;
    if (mBGn == 1)
        // assume tx bits length =  ceil(kBar/R), alternatively can use all layers(slower)
        nMaxLayer = ceil((ceil(mKBar/mR) + mF)/mZc) - 20;
    else{
        nMaxLayer = ceil((ceil(mKBar/mR) + mF)/mZc) - 8;
    }

    // initialize msg from check nodes to vector nodes, each edge correspond a message
    vector<vector<double>> CtoVMsg(mEdges.size());
    for(auto& e : CtoVMsg){
        e = vector<double>(mZc,0);
    }
    // llr updates
    unsigned nLayerEdges,edgeIdx, nShifts, vNodeIdx;
    for (unsigned iIter = 0; iIter < nMaxIter; iIter++){
        for(unsigned iLayer  = 0; iLayer < nMaxLayer; iLayer++ ){
            nLayerEdges = mLayers[iLayer].edgeEnd - mLayers[iLayer].edgeStart;
            // messages from variable nodes to check node
            vector<vector<double>> VtoCMsg(nLayerEdges);
            for(auto& e : VtoCMsg ){
                e = vector<double>(mZc,0);
            }
            for (unsigned iEdge = 0; iEdge < nLayerEdges; iEdge++){

                edgeIdx = mLayers[iLayer].edgeStart+ iEdge;
                vNodeIdx = mEdges[edgeIdx].vNodeIdx;
                nShifts = mShifts[edgeIdx];

                transform(LLR[vNodeIdx].begin(),LLR[vNodeIdx].end(),CtoVMsg[edgeIdx].begin(),LLR[vNodeIdx].begin(),minus<>{});
                VtoCMsg[iEdge] = LLR[vNodeIdx];
                rotate(VtoCMsg[iEdge].rbegin(),VtoCMsg[iEdge].rbegin()+nShifts,VtoCMsg[iEdge].rend());
            }
            //check node operation
            vector<vector<double>> minSumMsgs = checkNodeOperation(VtoCMsg);
            for (unsigned iEdge = 0; iEdge < nLayerEdges; iEdge++){
                edgeIdx = mLayers[iLayer].edgeStart+ iEdge;
                vNodeIdx = mEdges[edgeIdx].vNodeIdx;
                nShifts = mShifts[edgeIdx];

                CtoVMsg[edgeIdx] = minSumMsgs[iEdge];
                rotate(CtoVMsg[edgeIdx].rbegin(),CtoVMsg[edgeIdx].rbegin()+ (mZc-nShifts),CtoVMsg[edgeIdx].rend());
                transform(LLR[vNodeIdx].begin(),LLR[vNodeIdx].end(),CtoVMsg[edgeIdx].begin(),LLR[vNodeIdx].begin(),std::plus<>{});
            }
        }
    }
    // flatten the 2-D vector LLR
    vector<double> flattenLLR;
    for(auto e: LLR){
        flattenLLR.insert(flattenLLR.end(),e.begin(),e.end());
    }
    flattenLLR.erase(flattenLLR.end()-mF,flattenLLR.end());

    // chose information bits
    vector<bool> decBits(mKBar,0);
    transform(flattenLLR.begin(),flattenLLR.end(),decBits.begin(), [](double x){return x <= 0;});
    return decBits;
}
vector<vector<double>> nrLDPC::checkNodeOperation( vector<vector<double>> msgIn)
{
    // ------------------------------------------------------------------------------------------------
    // [ref] Chen, Jinghu, R.M. Tanner, C. Jones, and Yan Li. "Improved min-sum decoding algorithms for
    // irregular LDPC codes." In Proceedings. International Symposium on Information Theory, 2005.
    //-------------------------------------------------------------------------------------------------

    unsigned nNodes = msgIn.size();

    msgIn = reshape(msgIn,mZc,nNodes);
    vector<size_t> sortedIdx(nNodes);
    vector<double> sign(nNodes);
    double min1, min2, parity = 1.0;
    vector<vector<double>> msgOut(mZc);
    for(unsigned i = 0; i < mZc; i++){

        // sort abs(llr)
        sortedIdx = sort_indexes(msgIn[i]);
        transform(msgIn[i].begin(), msgIn[i].end(), sign.begin(), [](double x){return (x>=0)? 1.0: -1.0;});
        for(auto e: sign) parity *= e;

        //minimum and second minimum with offset
        min1 =  (abs(msgIn[i][sortedIdx[0]])> 0.5)? parity*abs(msgIn[i][sortedIdx[0]]): 0;
        min2 =  (abs(msgIn[i][sortedIdx[1]])> 0.5)? parity*abs(msgIn[i][sortedIdx[1]]): 0;

        // assign to output
        msgOut[i]= vector<double>(msgIn[i].size(),min1);
        msgOut[i][sortedIdx[0]] = min2;

        // add signs
        transform(msgOut[i].begin(), msgOut[i].end(), sign.begin(), msgOut[i].begin(), std::multiplies<>{});
    }

    return reshape(msgOut,nNodes,mZc);

}
vector<bool> nrLDPC::rateMatch(const vector<bool> &bitsIn, size_t nOfBitOut)
{
    if (mBGn == 1)
        assert(bitsIn.size() == 68 * mZc);
    else
        assert(bitsIn.size() == 52 * mZc);

    vector<bool> txBufferRing = bitsIn;

    // shortening by removing filling bits
    txBufferRing.erase(txBufferRing.begin()+mK-mF,txBufferRing.begin()+mK);

    // puncturing the first 2*Z
    txBufferRing.erase(txBufferRing.begin(),txBufferRing.begin()+2*mZc);

    // take nOfBitOut bits out of the ring
    vector<bool> bitsOut(nOfBitOut);
    for (unsigned i = 0,j=0; i< nOfBitOut; i++,j++){
        bitsOut[i] = txBufferRing[i%txBufferRing.size()];
    }
    return bitsOut;
}
vector<double> nrLDPC::rateRecover(const vector<double> &softBitsIn)
{
   unsigned rxRingLen = mN - 2*mZc - mF;
   vector<double> rxBufferRing(rxRingLen,0);

   // for received bits longer than the ring
   for(unsigned i=0; i < softBitsIn.size(); i++){
       rxBufferRing[i%rxRingLen] = rxBufferRing[i%rxRingLen] + softBitsIn[i];
   }


   // first 2*Zc with all 0
   vector<double> softBitsOut(2*mZc,0.0);
   // add information soft bits
   softBitsOut.insert(softBitsOut.end(), rxBufferRing.begin(), rxBufferRing.begin()+mKBar);
   //fillers
   vector<double> fillers(mF,numeric_limits<double>::infinity());
   softBitsOut.insert(softBitsOut.end(), fillers.begin(),fillers.end());
   // add parity soft bits
   softBitsOut.insert(softBitsOut.end(), rxBufferRing.begin()+mKBar, rxBufferRing.end());

   return softBitsOut;

}

uint8_t nrLDPC::selectBaseGraph(size_t KBar, double R)
{
    // 3GPP 38.212 7.2.2 LDPC base graph selection
    if( KBar <= 292 || (KBar <= 3824 && R <= 0.67) || R <= 0.25) {
        assert (KBar <= 3840);
        return 2;
    }
    else{
        assert( KBar <= 8448);
        return 1;
    }

}
uint16_t nrLDPC::selectLiftSize(size_t KBar, uint8_t BGn)
{
    // select kb 3GPP 38.212 section 5.2.2
    uint16_t Kb;
    if (BGn == 1)
        Kb = 22;
    else {
        if  (KBar > 640)
            Kb = 10;
        else if (KBar > 560)
            Kb = 9;
        else if (KBar > 192)
            Kb = 8;
        else
            Kb = 6;
    }

    // select lifting size and shifting set index based on 3GPP Table 5.3.2-1
    uint16_t  Zc = 384;
    uint16_t candiZc;
    for (unsigned i = 0; i < 8; i++){
        for(unsigned j = 0; j < 8; j++){
            candiZc = liftSizeTable[i][j];
            if (candiZc*Kb == KBar)
                return candiZc;
            else if(candiZc*Kb > KBar && candiZc < Zc)
                Zc = candiZc;
        }
    }
    return Zc;

}
uint8_t nrLDPC::selectShiftSet(uint16_t Zc)
{
    assert(Zc >=2 && Zc <= 384);
    for (unsigned i = 0; i < 8; i++){
        for(unsigned j = 0; j < 8; j++){
            if (liftSizeTable[i][j] == Zc)
                return i;
        }
    }
    cerr<<" Zc is not valid!";
    return -1;
}
vector<vector<bool>> nrLDPC::makeParityCheckMatrix(uint8_t BGn, const uint16_t Zc)
{

    uint8_t setIdx = selectShiftSet(Zc);
    unsigned numOfBlkRows,numOfBlkCols,numOfShifts;
    if(BGn == 1){
        numOfBlkRows = 46; numOfBlkCols = 68;
        numOfShifts = 316;
    }else{
        numOfBlkRows = 42; numOfBlkCols = 52;
        numOfShifts = 197;
    }

    // initialize parity check matrix (PCM)
    vector<vector<bool>> H(numOfBlkRows);
    for (unsigned i = 0; i <numOfBlkRows*Zc; i++){
        H[i] = vector<bool>(numOfBlkCols*Zc,0);
    }

    unsigned blkRowIdx, blkColIdx, shiftCoeff;
    // a shifted vecOne will be used to fill H
    // vecOne = [0,0,0,1] if Zc = 4;
    vector<bool> vecOne(Zc,0); vecOne.back() = 1;

    for (int i =0; i< numOfShifts; i++){

        if (BGn ==1) {
            blkRowIdx = shiftTableBgn_1[i][0];
            blkColIdx = shiftTableBgn_1[i][1];
            shiftCoeff = shiftTableBgn_1[i][2 + setIdx];
        }
        else{
            blkRowIdx = shiftTableBgn_2[i][0];
            blkColIdx = shiftTableBgn_2[i][1];
            shiftCoeff = shiftTableBgn_2[i][2 + setIdx];
        }

        // shift vector one
        vector<bool> shiftVecOne = vecOne;
        rotate(shiftVecOne.rbegin(),shiftVecOne.rbegin()+ shiftCoeff,shiftVecOne.rend());

        // shift right by 1 for each row in a block
        for (unsigned m =0; m< Zc; m++){
            rotate(shiftVecOne.rbegin(),shiftVecOne.rbegin()+1,shiftVecOne.rend());
            for(unsigned n = 0; n < Zc; n++){
                H[Zc*blkRowIdx+m][Zc*blkColIdx+n] = shiftVecOne[n];
            }
        }
    }

    return H;
}

template<typename T>
inline vector<size_t> nrLDPC::sort_indexes(const vector<T>& v)
{
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return abs(v[i1]) < abs(v[i2]); });

    return idx;
}

//template<typename T>
vector<vector<double>> nrLDPC::reshape(vector<vector<double>>& mat, const unsigned nRows, const unsigned nCols )
{
    vector<double> vec;
    for(auto e: mat)
        vec.insert(vec.end(),e.begin(),e.end());

    assert(vec.size() == nRows * nCols);

    vector<vector<double>> matOut(nRows);
    for(int i=0; i < nRows; i++){
        matOut[i] = vector<double>(vec.begin()+i*nCols, vec.begin()+(i+1)*nCols) ;
    }
    return matOut;
}