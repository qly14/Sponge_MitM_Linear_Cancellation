//
//  main.cpp
//  3r experiment on Ascon
//
//  Created on 2024/12/10.
//

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <cmath>
#include <vector>

using namespace std;

typedef unsigned char UINT8;
typedef unsigned long int UINT32;
typedef unsigned long long int UINT64;


#define random(x) (rand())%x;
#define nrRounds 2
#define nrLanes 5

#define ROL64(a, offset) ((offset != 0) ? ((((UINT64)a) << offset) ^ (((UINT64)a) >> (64-offset))) : a)
#define ROR64(a, offset) ((offset != 0) ? ((((UINT64)a) >> offset) ^ (((UINT64)a) << (64-offset))) : a)
#define TAKE_BIT(x, pos) (((x) >> (pos)) & 0x1)


unsigned int addition_indexm[18]={0,1,4,5,7,8,9,10,11,13,14,15,18,19,20,21,22,23};

unsigned int indexm[14]={2,3,6,12,16,17,26,34, 35, 38, 44, 48, 49, 58};
unsigned int indexb[14]={1,4,11,13,14,17,23,33, 36, 43, 45, 46, 49, 55};
unsigned int indexr[24]={0,5,6,7,15,18,22,25,26,27,28,29, 32, 37, 38, 39, 47, 50, 54, 57, 58, 59, 60, 61};
unsigned int indexrf[24]={5,7,15,18,22,26,27,32, 37, 38, 39, 47, 50, 54, 57, 58, 59, 61};
unsigned int cond0x[6] = {7, 17, 26,39 , 49 , 58};
unsigned int cond1x[18] = {0, 1, 4, 11, 13, 18, 22, 23, 25,32 , 33 , 36 , 43 , 45 , 50 , 54 , 55 , 57 };
unsigned int condplusx[24] = {0, 5, 6, 7, 11, 14, 15, 17, 22, 25, 26, 27,32 , 37 , 38 , 39 , 43 , 46 , 47 , 49 , 54 , 57 , 58 , 59};

void PermutationOnWords(UINT64 *state);
void pL(UINT64 *A);
void pS(UINT64 *A);


void PermutationOnWords(UINT64 state[], int round)
{
    unsigned int i;
    for(i=0; i<round; i++) {
        pS(state);
        pL(state);
    }
}


void pL(UINT64 *A)
{
    A[0] = A[0] ^ ROR64(A[0], 19) ^ ROR64(A[0], 28);
    A[1] = A[1] ^ ROR64(A[1], 61) ^ ROR64(A[1], 39);
    A[2] = A[2] ^ ROR64(A[2], 1) ^ ROR64(A[2], 6);
    A[3] = A[3] ^ ROR64(A[3], 10) ^ ROR64(A[3], 17);
    A[4] = A[4] ^ ROR64(A[4], 7) ^ ROR64(A[4], 41);
}


void pS(UINT64 *A)
{
    unsigned int x;
    UINT64 C[5];

    C[0] = (A[4]&A[1]) ^ (A[3]) ^ (A[2]&A[1]) ^ (A[2]) ^ (A[1]&A[0]) ^ (A[1]) ^ (A[0]);
    C[1] = (A[4]) ^ (A[3]&A[2]) ^ (A[3]&A[1]) ^ (A[3]) ^ (A[2]&A[1]) ^ (A[2]) ^ (A[1]) ^ (A[0]);
    C[2] = (A[4]&A[3]) ^ (A[4]) ^ (A[2]) ^ (A[1]) ^ (0xFFFFFFFFFFFFFFFF);
    C[3] = (A[4]&A[0]) ^ (A[4]) ^ (A[3]&A[0]) ^ (A[3]) ^ (A[2])^ (A[1]) ^ (A[0]);
    C[4] = (A[4]&A[1]) ^ (A[4]) ^ (A[3]) ^ (A[1]&A[0])^ (A[1]);
    for(x=0; x<5; x++)
        A[x] = C[x];
}


void displaystate(UINT64 *state)
{
    unsigned int i;
    for(i=0;i<5;i++)
    {
        printf("%016llx ",(state[i]));
        printf("\n");
    }
    printf("\n");
}

UINT64 generateMP(UINT64 *state)
{
    UINT64 MP=0;
    for (int i=0; i<14; i++) {
        MP = (MP << 1) ^ (TAKE_BIT(state[4], (63-indexm[i])) & TAKE_BIT(state[1], (63-indexm[i]))) ^ TAKE_BIT(state[3], (63-indexm[i])) ^ (TAKE_BIT(state[2], (63-indexm[i])) & TAKE_BIT(state[1], (63-indexm[i])))^
        TAKE_BIT(state[2], (63-indexm[i])) ^ (TAKE_BIT(state[1], (63-indexm[i])) & TAKE_BIT(state[0], (63-indexm[i]))) ^ TAKE_BIT(state[1], (63-indexm[i])) ^ TAKE_BIT(state[0], (63-indexm[i]));
    }
    return MP;
}

UINT64 generateConsume(UINT64 *state)
{
    UINT64 Consume=0;
    //Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-2));
    Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-6));
    //Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-12));
    //Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-16));
    Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-26));

    //Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-34));
    Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-38));
    //Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-44));
    //Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-48));
    Consume = (Consume << 1) ^ TAKE_BIT(state[1], (63-58));

    return Consume;
}

UINT64 generateVerify(UINT64 *state)
{
    UINT64 Verify=0;
    for (int i=0; i<14; i++) {
        Verify = (Verify << 1) ^ (TAKE_BIT(state[4], (63-indexm[i])) & TAKE_BIT(state[1], (63-indexm[i]))) ^ TAKE_BIT(state[3], (63-indexm[i])) ^ (TAKE_BIT(state[2], (63-indexm[i])) & TAKE_BIT(state[1], (63-indexm[i])))^
        TAKE_BIT(state[2], (63-indexm[i])) ^ (TAKE_BIT(state[1], (63-indexm[i])) & TAKE_BIT(state[0], (63-indexm[i]))) ^ TAKE_BIT(state[1], (63-indexm[i])) ^ TAKE_BIT(state[0], (63-indexm[i]));
    }
     for (int i=0; i<18; i++) {
        Verify = (Verify << 1) ^ (TAKE_BIT(state[4], (63-addition_indexm[i])) & TAKE_BIT(state[1], (63-addition_indexm[i]))) ^ TAKE_BIT(state[3], (63-addition_indexm[i])) ^ (TAKE_BIT(state[2], (63-addition_indexm[i])) & TAKE_BIT(state[1], (63-addition_indexm[i])))^
        TAKE_BIT(state[2], (63-addition_indexm[i])) ^ (TAKE_BIT(state[1], (63-addition_indexm[i])) & TAKE_BIT(state[0], (63-addition_indexm[i]))) ^ TAKE_BIT(state[1], (63-addition_indexm[i])) ^ TAKE_BIT(state[0], (63-addition_indexm[i]));
    }
    return Verify;
}

int main(int argc, const char * argv[])
{
    clock_t start_time = clock();
    srand((unsigned)time(NULL));
    
    UINT64 InitialState[5]={0};
    UINT64 TempState[5]={0};
    UINT64 MatchNum = 0;
    UINT64 SucceedNum = 0;
    vector<UINT64> PreImage;

    
    //Init the state with 0
    for(UINT64 i=0;i<5;i++){
        InitialState[i]=0;
    }

    //Set condition gray bits
    for(UINT64 j=0;j<18;j++){
        InitialState[1] |= (UINT64(1)<<(63-cond1x[j]));
        //InitialState[1] |= (UINT64(1)<<(63-cond1x[j]-32));
    }
    for(UINT64 j=0;j<24;j++){
        InitialState[3] |= (UINT64(1)<<(63-condplusx[j]));
        //InitialState[3] |= (UINT64(1)<<(63-condplusx[j]-32));
    }
    //padding
    InitialState[0] |= (UINT64(1)<<(63-62));
    InitialState[0] &= ROL64(~UINT64(1),(63-63));
    
    for(UINT64 Y=0;Y<(UINT64(1)<<2);Y++) {
        
        
        map<UINT64, vector<UINT64>> TableLConsume;
        
        UINT64 Statered = 0;
        UINT64 Fixred_e = 0;
        UINT64 FM3 = 0;
        cout << "Start Build Table Consume!\n";
        //2^{18} for one preimage
        for(UINT64 j=0;j<(UINT64(1)<<18);j++){
            for(UINT64 k=0;k<5;k++){
                TempState[k]=InitialState[k];
            }
            
            //Set blue=0
            for(UINT64 k=0;k<14;k++){
                TempState[0] &= ROL64(~UINT64(1),(63-indexb[k]));
            }
            
            //Traverse the 18 free red bits
            for(UINT64 k=0;k<18;k++) {
                UINT64 temp1=(j>>(17-k))&1;
                if(temp1){
                    TempState[0] |= (UINT64(1)<<(63-indexrf[k]));
                }
                else {
                    TempState[0] &= ROL64(~UINT64(1),(63-indexrf[k]));
                }
            }
            //Compute 6 linear cancellations
            UINT64 temp1=((TempState[0]>>(63-27))&1)^((TempState[0]>>(63-38))&1)^((TempState[0]>>(63-47))&1)^((TempState[0]>>(63-50))&1)^((Y>>(5))&1);
            if(temp1){
                TempState[0] |= (UINT64(1)<<(63-25));
            }
            else {
                TempState[0] &= ROL64(~UINT64(1),(63-25));
            }
            temp1=((TempState[0]>>(63-15))&1)^((Y>>(4))&1);
            if(temp1){
                TempState[0] |= (UINT64(1)<<(63-60));
            }
            else {
                TempState[0] &= ROL64(~UINT64(1),(63-60));
            }
            temp1=((TempState[0]>>(63-61))&1)^((Y>>(3))&1);
            if(temp1){
                TempState[0] |= (UINT64(1)<<(63-0));
            }
            else {
                TempState[0] &= ROL64(~UINT64(1),(63-0));
            }
            temp1=((TempState[0]>>(63-15))&1)^((TempState[0]>>(63-18))&1)^((TempState[0]>>(63-57))&1)^((TempState[0]>>(63-59))&1)^((Y>>(2))&1);
            if(temp1){
                TempState[0] |= (UINT64(1)<<(63-6));
            }
            else {
                TempState[0] &= ROL64(~UINT64(1),(63-6));
            }
            temp1=((TempState[0]>>(63-47))&1)^((Y>>(1))&1);
            if(temp1){
                TempState[0] |= (UINT64(1)<<(63-28));
            }
            else {
                TempState[0] &= ROL64(~UINT64(1),(63-28));
            }
            temp1=((TempState[0]>>(63-32))&1)^((Y>>(0))&1);
            if(temp1){
                TempState[0] |= (UINT64(1)<<(63-29));
            }
            else {
                TempState[0] &= ROL64(~UINT64(1),(63-29));
            }
            
            
            Statered = TempState[0];
            
            PermutationOnWords(TempState, 2);
            UINT64 tmpConsume=generateConsume(TempState);
            UINT64 tmpMatch = generateMP(TempState);
            if(TableLConsume.find(tmpConsume) != TableLConsume.end())
            {
                vector<UINT64> ttmp =TableLConsume[tmpConsume];
                ttmp.push_back(Statered);
                ttmp.push_back(tmpMatch);
                TableLConsume[tmpConsume] = ttmp;
            } else {
                vector<UINT64> ttmp;
                ttmp.push_back(Statered);
                ttmp.push_back(tmpMatch);
                TableLConsume[tmpConsume] = ttmp;
            }
        }
        cout << "Finish Build Table Consume!\n";
        // for CR
        for(UINT64 cr=0;cr<(UINT64(1)<<4);cr++){
            if(TableLConsume.find(cr) != TableLConsume.end()) {
                map<UINT64, vector<UINT64>> TableL1;
                vector<UINT64> CandM = TableLConsume[cr];
                Fixred_e = CandM[0];
                FM3 = CandM[1];
                //cout << "CR: "<< cr <<",  "<< CandM.size() <<"  "<<CandM[0]<<"  "<<  CandM[1]<< endl;
                for (UINT64 k=0; k<CandM.size(); k+=2) {
                    if(TableL1.find(CandM[k+1]) != TableL1.end())
                    {
                        vector<UINT64> ttmp =TableL1[CandM[k+1]];
                        ttmp.push_back(CandM[k]);
                        TableL1[CandM[k+1]] = ttmp;
                    } else {
                        vector<UINT64> ttmp;
                        ttmp.push_back(CandM[k]);
                        TableL1[CandM[k+1]] = ttmp;
                    }
                }
                //Fix Red to e, traverse the blue bits
                for(UINT64 j=0; j<(UINT64(1)<<14); j++){
                    UINT64 Stateblue = 0;
                    for(UINT64 k=0; k<5; k++){
                        TempState[k]=InitialState[k];
                    }
                    //Set red=Fixred_e
                    TempState[0] = Fixred_e;
                    
                    for(UINT64 k=0;k<14;k++)
                    {
                        UINT64 temp1=(j>>(13-k))&1;
                        if(temp1){
                            TempState[0] |= (UINT64(1)<<(63-indexb[k]));
                            Stateblue |= (UINT64(1)<<(63-indexb[k]));
                        }
                        else{
                            TempState[0] &= ROL64(~UINT64(1),(63-indexb[k]));
                            Stateblue &= ROL64(~UINT64(1),(63-indexb[k]));
                        }
                    }
                    
                    PermutationOnWords(TempState, 2);
                    
                    UINT64 Matchfblue = 0;
                    //Compute fb
                    Matchfblue = generateMP(TempState) ^ FM3;
                    
                    if(TableL1.find(Matchfblue) != TableL1.end()) {
                        vector<UINT64> ttmp = TableL1[Matchfblue];
                        MatchNum += ttmp.size();
                        
                        //verify
                        for (UINT64 k=0; k<ttmp.size(); k++) {
                            UINT64 Message = ttmp[k]^Stateblue;
                            TempState[0] = Message;
                            for(int i=1; i<5; i++){
                                TempState[i]=InitialState[i];
                            }
                            
                            PermutationOnWords(TempState, 2);
                            //pS(TempState);
                            UINT64 hash = generateVerify(TempState);
                            if (hash==0) {
                                SucceedNum += 1;
                                PreImage.push_back(Message);
                                //if(outputNum<10){
                                //outputNum+=1;
                                //printf("%016llx \n",(Message));
                                // for(int i=1; i<5; i++){
                                //    printf("%016llx \n",(InitialState[i]));
                                // }
                                //}
                            }
                        }
                        
                    }
                }
            }
        }
    }
    
    
    
    
    cout << "Finish Mitm Attack\n";
    
    clock_t end_time = clock();
    cout << "The time of Mitm Attack is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

    cout << "In total, 2^" << log(double(MatchNum))/log(2.0) << " match are tested!" << endl;
    cout << "In total, 2^" << log(double(SucceedNum))/log(2.0) << " success are found!" << endl;
    for (UINT64 k=0; k<PreImage.size(); k++) {
        printf("%016llx \n",(PreImage[k]));
    }
    return 0;

}



