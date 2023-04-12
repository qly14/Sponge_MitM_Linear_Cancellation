//
//  main.cpp
//  2r experiments on Ascon-XoF
//  
//  Created on 2022/9/3.
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

unsigned int indexm[18]={11,14,16,18,22,27,29,34,35,36,37,38,44,45,46,48,62,63};
unsigned int indexb[17]={0,6,7,8,11,12,13,14,17,18,21,23,24,26,27,30,31};
unsigned int indexr[17]={1,2,3,4,5,9,10,15,16,19,20,22,25,28,29,32,33};
    

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
    unsigned int x, y;
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
    for (int i=0; i<18; i++) {
        MP = (MP << 1) ^ TAKE_BIT(state[0], (63-indexm[i]));
    }
    return MP;
}


int main(int argc, const char * argv[])
{
    srand((unsigned)time(NULL));
    
    UINT64 InitialState[5]={0};
    UINT64 TempState[5]={0};
    //UINT64 FinalState=0;
    UINT64 PreNum=0;

    
    //Init the state with 0
    for(UINT64 i=0;i<5;i++){
        InitialState[i]=0;
    }
    
    
    map<UINT64, vector<UINT64>> TableH;
    for (UINT64 times = 0; times < 32; times ++) {
        //Randomly set gray bits
        for(UINT64 i=34;i<64;i++){
            UINT64 temp=random(2);
            if(temp){
                InitialState[0] |= (UINT64(1)<<(63-i));
            }
            else {
                InitialState[0] &= ROL64(~UINT64(1),(63-i));
            }
        }
        UINT64 Stategray=InitialState[0];
        //Compute fg
        for(UINT64 i=0;i<5;i++){
            TempState[i]=InitialState[i];
        }
        PermutationOnWords(TempState, 1);
        pS(TempState);
        UINT64 Matchfgray = generateMP(TempState);
        //printf("%016llx \n",(Matchfg));
        //displaystate(TempState);
   
        map<UINT64, vector<UINT64>> TableL1;
        //Generate table L1
        UINT64 Statered = 0;
        for(UINT64 j=0;j<(UINT64(1)<<17);j++){
            for(UINT64 k=0;k<5;k++){
                    TempState[k]=InitialState[k];
            }
       
            //Set blue=0
            for(UINT64 k=0;k<17;k++){
                TempState[0] &= ROL64(~UINT64(1),(63-indexb[k]));
            }

            //Traverse the red bits
            for(UINT64 k=0;k<17;k++) {
                UINT64 temp1=(j>>(16-k))&1;
                if(temp1){
                    TempState[0] |= (UINT64(1)<<(63-indexr[k]));
                }
                else {
                    TempState[0] &= ROL64(~UINT64(1),(63-indexr[k]));
                }                   
            }

            Statered = TempState[0];
            //Statered = generateRed(TempState);
            //printf("%016llx \n",(Statered));
            //Compute S^(1)
            PermutationOnWords(TempState, 1);
            pS(TempState);
        
            UINT64 Matchfrg = 0;
            //Compute f'm
            Matchfrg = generateMP(TempState);
            //printf("%016llx \n",(Matchfrg));
            //Store the value for 17 red bits
            if(TableL1.find(Matchfrg) != TableL1.end())
            {
                vector<UINT64> ttmp =TableL1[Matchfrg];
                ttmp.push_back(Statered);
                TableL1[Matchfrg] = ttmp;
            } else {
                vector<UINT64> ttmp;
                ttmp.push_back(Statered);
                TableL1[Matchfrg] = ttmp;
            }        
        }
        //cout << "Build L1!\n";
        
        UINT64 Stateblue = 0;
        //Fix Red to 0, traverse the blue bits
        for(UINT64 j=0; j<(UINT64(1)<<17); j++){
            for(UINT64 k=0; k<5; k++){
                TempState[k]=InitialState[k];
            }
            //Set red=0
            for(UINT64 k=0;k<17;k++){
                TempState[0] &= ROL64(~UINT64(1),(63-indexr[k]));
            }
            for(UINT64 k=0;k<17;k++)
            {
                UINT64 temp1=(j>>(16-k))&1;
                if(temp1){
                    TempState[0] |= (UINT64(1)<<(63-indexb[k]));
                }
                else{
                    TempState[0] &= ROL64(~UINT64(1),(63-indexb[k]));
                }                   
            }
            //Stateblue = generateBlue(TempState);
            Stateblue = TempState[0];
            //printf("%016llx \n",(Stateblue));
            PermutationOnWords(TempState, 1);
            pS(TempState);
        
            UINT64 Matchfblue = 0;
            //Compute f''m
            Matchfblue = generateMP(TempState) ^ Matchfgray;

            if(TableL1.find(Matchfblue) != TableL1.end()) {
                vector<UINT64> ttmp = TableL1[Matchfblue];
            
                for (UINT64 k=0; k<ttmp.size(); k++) {
                    UINT64 Message = ttmp[k]^Stateblue^Stategray;
                    TempState[0] = Message;
                    for(int i=1; i<5; i++)
                        TempState[i]=0;
                 
                    PermutationOnWords(TempState, 1);
                    pS(TempState);  
                    UINT64 hash = generateMP(TempState);
                    if (hash==0) {
                        PreNum += ttmp.size();
                        
                        if(TableH.find(TempState[0]) != TableH.end())
                        {
                            //cout << "Find the collision, !" << endl;
                            printf("Find the collision hash value:%016llx \n",(TempState[0]));
                            vector<UINT64> ttmp1 =TableH[TempState[0]];
                            for (UINT64 k1=0; k1<ttmp1.size(); k1++) {
                                printf("Collision at %016llx ,",(Message));
                                printf("%016llx \n",(ttmp1[k1]));
                            }
                            ttmp1.push_back(Message);
                            TableH[TempState[0]] = ttmp1;
                        } else {
                            vector<UINT64> ttmp1;
                            ttmp1.push_back(Message);
                            TableH[TempState[0]] = ttmp1;
                        }  
                    }   
                }
            }
        }
    }
    //cout << "In total, " << MatchNum << " matches are found!" << endl;
    cout << "In total, 2^" << log(double(PreNum))/log(2.0) << " preimages are found!" << endl;
    return 0;

}


