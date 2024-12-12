#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <tuple>
#include <vector>
#include <map>
#include <cassert>
#include <cmath>

typedef unsigned char tKeccakLane;
typedef unsigned int UINT32;

UINT32 lengthLane = 8;
UINT32 offset[] = {
	0, 1, 6, 4, 3,
	4, 4, 6, 7, 4,
	3, 2, 3, 1, 7,
	1, 5, 7, 5, 0,
	2, 2, 5, 0, 6
};

#define nrLanes 25
#define index(x, y) (((x)%5)+5*((y)%5))
#define GET_BIT(A, pos, length) (((A) >> (length - pos - 1)) & 0x1)
#define ROL(a, offset, length) ((offset != 0) ? ((((tKeccakLane)a) >> offset) ^ (((tKeccakLane)a) << (length-offset))) : a)

using namespace std;

void SetBits(tKeccakLane* A, tuple<UINT32, UINT32, UINT32> pos[], UINT32 value, UINT32 length) {
	for (int i = 0; i < length; i++)
	{
		UINT32 x = get<0>(pos[i]);
		UINT32 y = get<1>(pos[i]);
		UINT32 z = get<2>(pos[i]);
		UINT32 v = ((value >> (length - 1 - i)) & 0x1);
		A[index(x, y)] |= (v << (lengthLane - 1 - z));
	}
}

void getCons(UINT32* cons, tKeccakLane* theta1, tuple <UINT32, UINT32, UINT32> cons_pos3[]) {
	UINT32 x, y, z;
	for (int i = 0; i < 12; i++)
	{
		x = get<0>(cons_pos3[i]);
		y = get<1>(cons_pos3[i]);
		z = get<2>(cons_pos3[i]);
		*cons = ((*cons) << 1) | GET_BIT(theta1[index(x, y)], z, lengthLane);
	}
}

void getMatch(tKeccakLane* A, UINT32* match, tuple <UINT32, UINT32, UINT32> match_pos1[], tuple <UINT32, UINT32, UINT32> match_pos2[], bool verb=0) {
	UINT32 x1, y1, z1, x2, y2, z2;
	for (int i = 0; i < 3; i++)
	{
		x1 = get<0>(match_pos1[i]);
		y1 = get<1>(match_pos1[i]);
		z1 = get<2>(match_pos1[i]);
		x2 = get<0>(match_pos2[i]);
		y2 = get<1>(match_pos2[i]);
		z2 = get<2>(match_pos2[i]);
		if (verb)
			cout << GET_BIT(A[index(x1, y1)], z1, lengthLane) << GET_BIT(A[index(x2, y2)], z2, lengthLane) << endl;
		*match = ((*match) << 1) | ((GET_BIT(A[index(x1, y1)], z1, lengthLane) ^ GET_BIT(A[index(x2, y2)], z2, lengthLane)) & 0x1);
	}
}

void theta(tKeccakLane* A) {
	unsigned int x, y;
	tKeccakLane C[5], D[5];
	for (x = 0; x < 5; x++) {
		C[x] = 0;
		for (y = 0; y < 5; y++)
			C[x] ^= A[index(x, y)];
	}
	//for (int i = 0; i < 5; i++)
	//{
	//	printf("%02x  ", C[i]);
	//}
	//printf("\n");
	for (x = 0; x < 5; x++)
		D[x] = ROL(C[(x + 1) % 5], 1, lengthLane) ^ C[(x + 4) % 5];
	//for (int i = 0; i < 5; i++)
	//{
	//	printf("%02x  ", D[i]);
	//}
	//printf("\n");
	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			A[index(x, y)] ^= D[x];
}
void rho(tKeccakLane* A) {
	unsigned int x, y;
	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			A[index(x, y)] = ROL(A[index(x, y)], offset[index(x, y)], lengthLane);
}
void pi(tKeccakLane* A) {
	unsigned int x, y;
	tKeccakLane tempA[25];
	for (x = 0; x < 5; x++) {
		for (y = 0; y < 5; y++)
			tempA[index(x, y)] = A[index(x, y)];
	}
	for (x = 0; x < 5; x++) {
		for (y = 0; y < 5; y++)
			A[index(0 * x + 1 * y, 2 * x + 3 * y)] = tempA[index(x, y)];
	}
}
void chi(tKeccakLane* A) {
	unsigned int x, y;
	tKeccakLane C[5];
	for (y = 0; y < 5; y++) {
		for (x = 0; x < 5; x++)
			C[x] = A[index(x, y)] ^ ((~A[index(x + 1, y)]) & A[index(x + 2, y)] & 0xff);
		for (x = 0; x < 5; x++)
			A[index(x, y)] = C[x];
	}
}
void Permutation(tKeccakLane* A, int round) {
	for (int r = 0; r < round; r++)
	{
		theta(A);
		rho(A);
		pi(A);
		chi(A);
	}
}

void SolveLinearCons(UINT32 lCons, tuple<UINT32, UINT32, UINT32> red_pos[], vector<UINT32>& redVs) {
	map<tuple<UINT32, UINT32, UINT32>, bool> redV;
	tuple<UINT32, UINT32, UINT32> red_active_pos[15] = {
		{0, 0, 0}, {1, 0, 0}, {3, 0, 0}, {4, 0, 0},
		{3, 0, 1},
		{0, 0, 2}, {2, 0, 2},
		{3, 0, 3},
		{1, 0, 4},
		{2, 0, 5},
		{0, 0, 6}, {2, 0, 6}, {3, 0, 6}, {4, 0, 6},
		{4, 0, 7},
	};
	tuple<UINT32, UINT32, UINT32> red_inactive_pos[8] = {
		{1, 0, 1}, {2, 0, 1}, {4, 0, 1},
		{4, 0, 5},
		{0, 0, 7}, {1, 0, 7}, {2, 0, 7}, {3, 0, 7},
	};
	for (UINT32 i = 0; i < (1 << 15); i++)
	{
		for (UINT32 j = 0; j < 15; j++)
		{
			redV[red_active_pos[j]] = ((i >> (14 - j)) & 0x1);
		}
		redV[{2, 0, 7}] = redV[{0, 0, 0}] ^ GET_BIT(lCons, 0, lengthLane);
		redV[{3, 0, 7}] = redV[{1, 0, 0}] ^ GET_BIT(lCons, 1, lengthLane);
		redV[{0, 0, 7}] = redV[{3, 0, 0}] ^ GET_BIT(lCons, 2, lengthLane);
		redV[{2, 0, 1}] = redV[{0, 0, 2}] ^ GET_BIT(lCons, 3, lengthLane);
		redV[{4, 0, 1}] = redV[{2, 0, 2}] ^ GET_BIT(lCons, 4, lengthLane);
		redV[{1, 0, 1}] = redV[{3, 0, 0}] ^ redV[{0, 0, 2}] ^ GET_BIT(lCons, 3, lengthLane) ^ GET_BIT(lCons, 5, lengthLane);
		redV[{4, 0, 5}] = redV[{3, 0, 6}] ^ redV[{2, 0, 6}] ^ GET_BIT(lCons, 6, lengthLane);
		redV[{1, 0, 7}] = redV[{2, 0, 6}] ^ redV[{3, 0, 0}] ^ GET_BIT(lCons, 2, lengthLane) ^ GET_BIT(lCons, 7, lengthLane);
		UINT32 tmp = 0;
		for (UINT32 j = 0; j < 23; j++)
		{
			tmp = ((tmp << 1) | redV[red_pos[j]]);
		}
		redVs.push_back(tmp);
	}
}

void getTmpM(UINT32* tmpM, UINT32 redV, tuple<UINT32, UINT32, UINT32> red_pos[], tuple <UINT32, UINT32, UINT32> match_pos1[], tuple <UINT32, UINT32, UINT32> match_pos2[]) {
	tKeccakLane initial_red_state[25] = { 0 };
	SetBits(initial_red_state, red_pos, redV, 23);
	Permutation(initial_red_state, 2);
	getMatch(initial_red_state, tmpM, match_pos1, match_pos2);
}

void showState(tKeccakLane* A) {
	for (int j = 0; j < 5; j++)
	{
		for (int i = 0; i < 5; i++)
		{
			printf("%02x, ", A[index(i, j)]);
		}
		printf("\n");
	}
	printf("\n");
}

bool verify_3bit(tKeccakLane* A, bool verb = 0) {
	tKeccakLane tmp[25] = {0};
	for (int i = 0; i < 25; i++)
	{
		tmp[i] = A[i];
	}
	tuple <UINT32, UINT32, UINT32> match_pos1[3] = { {2, 2, 1}, {1, 1, 3}, {0, 0, 7} };
	tuple <UINT32, UINT32, UINT32> match_pos2[3] = { {2, 4, 1}, {1, 3, 3}, {0, 2, 7} };
	Permutation(tmp, 2);
	UINT32 match = 0;
	getMatch(tmp, &match, match_pos1, match_pos2, verb);
	if (match == 0)
		return 1;
	else
		return 0;
}

bool attack_24bit(tKeccakLane* A, bool verb = 0) {
	tuple <UINT32, UINT32, UINT32> pi_pos[24] = {
		{0, 0, 4}, {0, 0, 5}, {0, 0, 6}, {0, 0, 7},
		{2, 1, 0}, {2, 1, 1}, {2, 1, 2}, {2, 1, 3},
		{1, 0, 4}, {1, 0, 5}, {1, 0, 6}, {1, 0, 7},
	    {3, 1, 1}, {3, 1, 2}, {3, 1, 3}, {3, 1, 4},
		{2, 0, 4}, {2, 0, 5}, {2, 0, 6}, {2, 0, 7},
		{4, 1, 4}, {4, 1, 5}, {4, 1, 6}, {4, 1, 7},
	};
	tKeccakLane tmp[25] = { 0 };
	for (int i = 0; i < 25; i++)
	{
		tmp[i] = A[i];
	}
	Permutation(tmp, 2);
	theta(tmp);
	rho(tmp);
	pi(tmp);
	UINT32 x, y, z;
	bool flag = 1;
	for (int i = 0; i < 24; i++)
	{
		x = get<0>(pi_pos[i]);
		y = get<1>(pi_pos[i]);
		z = get<2>(pi_pos[i]);
		if (GET_BIT(tmp[index(x, y)], z, lengthLane) != 0) {
			flag = 0;
			break;
		}
	}
	if (flag) {
		if (verb)
			showState(tmp);
		return 1;
	}
	else
		return 0;
}

void keccak_attack_24bit() {
	tuple <UINT32, UINT32, UINT32> blue_pos[3] = { {4, 0, 3}, {2, 0, 4}, {1, 0, 5} };
	tuple <UINT32, UINT32, UINT32> red_pos[23] = {
		{0, 0, 0}, {1, 0, 0}, {3, 0, 0}, {4, 0, 0},
		{1, 0, 1}, {2, 0, 1}, {3, 0, 1}, {4, 0, 1},
		{0, 0, 2}, {2, 0, 2},
		{3, 0, 3},
		{1, 0, 4},
		{2, 0, 5}, {4, 0, 5},
		{0, 0, 6}, {2, 0, 6}, {3, 0, 6}, {4, 0, 6},
		{0, 0, 7}, {1, 0, 7}, {2, 0, 7}, {3, 0, 7}, {4, 0, 7},
	};
	tuple <UINT32, UINT32, UINT32> nl_cons_pos[12] = {
		{0, 2, 0}, {0, 3, 0},
		{1, 2, 1}, {2, 1, 1}, {4, 3, 1},
		{1, 1, 3},
		{2, 2, 4}, {2, 3, 4},
		{1, 3, 6}, {3, 2, 6},
		{0, 4, 7}, {1, 4, 7}
	};
	tuple <UINT32, UINT32, UINT32> match_pos1[3] = { {2, 2, 1}, {1, 1, 3}, {0, 0, 7} };
	tuple <UINT32, UINT32, UINT32> match_pos2[3] = { {2, 4, 1}, {1, 3, 3}, {0, 2, 7} };
	UINT32 counter = 0;

	// Attack start
	clock_t start_time = clock();
	for (UINT32 lCons = 0; lCons < (UINT32(1) << 8); lCons++) // Linear Constraints
	{
		// Solve the linear constraints
		vector<UINT32> redVs;
		SolveLinearCons(lCons, red_pos, redVs);
		if (counter == 0)
			cout << "After the filter of linear constraints, 2^" << log(double(redVs.size())) / log(2.0) << " solutions are left." << endl;


		// Build Table U
		map<UINT32, vector<UINT32>> TableU; // Hash table U[non-linear cons] -> red cell values, 2^3 under each index
		TableU.clear();
		for (auto redV : redVs) {
			tKeccakLane tmp_red_state[25] = { 0 };
			SetBits(tmp_red_state, red_pos, redV, 23);
			tKeccakLane C[5], D[5];
			UINT32 x, y;
			for (x = 0; x < 5; x++) {
				C[x] = 0;
				for (y = 0; y < 5; y++)
					C[x] ^= tmp_red_state[index(x, y)];
			}
			for (x = 0; x < 5; x++)
				D[x] = ROL(C[(x + 1) % 5], 1, lengthLane) ^ C[(x + 4) % 5];
			for (x = 0; x < 5; x++)
				for (y = 0; y < 5; y++)
					tmp_red_state[index(x, y)] ^= D[x];
			rho(tmp_red_state);
			pi(tmp_red_state);
			chi(tmp_red_state);
			theta(tmp_red_state);
			UINT32 cons = 0;
			getCons(&cons, tmp_red_state, nl_cons_pos);
			TableU[(cons & 0xfff)].push_back(redV);
		}


		for (UINT32 nlCons = 0; nlCons < (1 << 12); nlCons++) // Non-linear Constraints
		{
			UINT32 tmpM = 0;
			if (TableU[nlCons].size() > 0) {
				getTmpM(&tmpM, TableU[nlCons][0], red_pos, match_pos1, match_pos2);
			}
			else
			{
				continue;
			}

			// Build Table L
			map<UINT32, vector<UINT32>> TableL;
			TableL.clear();
			for (auto redV : TableU[nlCons]) // 2^3
			{
				//cout << redV << endl;
				tKeccakLane initial_red_state[25] = { 0 };
				SetBits(initial_red_state, red_pos, redV, 23);
				Permutation(initial_red_state, 2);
				UINT32 match_red = 0;
				getMatch(initial_red_state, &match_red, match_pos1, match_pos2);
				TableL[(match_red) & 0x7].push_back(redV);
			}

			// blue cell computation
			for (UINT32 blueV = 0; blueV < (1 << 3); blueV++) // 2^3
			{
				tKeccakLane initial_blue_state[25] = { 0 };
				SetBits(initial_blue_state, blue_pos, blueV, 3);
				SetBits(initial_blue_state, red_pos, TableU[nlCons][0], 23);
				Permutation(initial_blue_state, 2);
				UINT32 match_blue = 0;
				getMatch(initial_blue_state, &match_blue, match_pos1, match_pos2);
				for (auto redV_match : TableL[match_blue ^ tmpM]) {
					tKeccakLane initial_state[25] = { 0 };
					SetBits(initial_state, blue_pos, blueV, 3);
					SetBits(initial_state, red_pos, redV_match, 23);
					if (attack_24bit(initial_state)) {
						if (counter <= 10) {
							cout << "Solution" << counter << endl;
							showState(initial_state);
							attack_24bit(initial_state, 1);
							
						}
						counter += 1;
					}
				}
			}
		}
	}
	// Attack end
	clock_t end_time = clock();
	cout << "\nThe time of finding all solutions satisfying the 24-bit partial target: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	cout << "\nThe number of all finding solutions: 2^" << log(double(counter)) / log(2.0) << endl;
}

void keccak_attack_3bit() {
	tuple <UINT32, UINT32, UINT32> blue_pos[3] = { {4, 0, 3}, {2, 0, 4}, {1, 0, 5} };
	tuple <UINT32, UINT32, UINT32> red_pos[23] = {
		{0, 0, 0}, {1, 0, 0}, {3, 0, 0}, {4, 0, 0},
		{1, 0, 1}, {2, 0, 1}, {3, 0, 1}, {4, 0, 1},
		{0, 0, 2}, {2, 0, 2},
		{3, 0, 3},
		{1, 0, 4},
		{2, 0, 5}, {4, 0, 5},
		{0, 0, 6}, {2, 0, 6}, {3, 0, 6}, {4, 0, 6},
		{0, 0, 7}, {1, 0, 7}, {2, 0, 7}, {3, 0, 7}, {4, 0, 7},
	};
	tuple <UINT32, UINT32, UINT32> nl_cons_pos[12] = {
		{0, 2, 0}, {0, 3, 0},
		{1, 2, 1}, {2, 1, 1}, {4, 3, 1},
		{1, 1, 3},
		{2, 2, 4}, {2, 3, 4},
		{1, 3, 6}, {3, 2, 6},
		{0, 4, 7}, {1, 4, 7}
	};
	tuple <UINT32, UINT32, UINT32> match_pos1[3] = { {2, 2, 1}, {1, 1, 3}, {0, 0, 7} };
	tuple <UINT32, UINT32, UINT32> match_pos2[3] = { {2, 4, 1}, {1, 3, 3}, {0, 2, 7} };
	UINT32 counter = 0;
	UINT32 error_counter = 0;

	// Attack start
	clock_t start_time = clock();
	for (UINT32 lCons = 0; lCons < (UINT32(1) << 8); lCons++) // Linear Constraints
	{
		// Solve the linear constraints
		vector<UINT32> redVs;
		SolveLinearCons(lCons, red_pos, redVs);
		if (counter == 0)
			cout << "After the filter of linear constraints, 2^" << log(double(redVs.size())) / log(2.0) << " solutions are left." << endl;


		// Build Table U
		map<UINT32, vector<UINT32>> TableU; // Hash table U[non-linear cons] -> red cell values, 2^3 under each index
		TableU.clear();
		for (auto redV : redVs) {
			tKeccakLane tmp_red_state[25] = { 0 };
			SetBits(tmp_red_state, red_pos, redV, 23);
			tKeccakLane C[5], D[5];
			UINT32 x, y;
			for (x = 0; x < 5; x++) {
				C[x] = 0;
				for (y = 0; y < 5; y++)
					C[x] ^= tmp_red_state[index(x, y)];
			}
			for (x = 0; x < 5; x++)
				D[x] = ROL(C[(x + 1) % 5], 1, lengthLane) ^ C[(x + 4) % 5];
			for (x = 0; x < 5; x++)
				for (y = 0; y < 5; y++)
					tmp_red_state[index(x, y)] ^= D[x];
			rho(tmp_red_state);
			pi(tmp_red_state);
			chi(tmp_red_state);
			theta(tmp_red_state);
			UINT32 cons = 0;
			getCons(&cons, tmp_red_state, nl_cons_pos);
			TableU[(cons & 0xfff)].push_back(redV);
		}


		for (UINT32 nlCons = 0; nlCons < (1 << 12); nlCons++) // Non-linear Constraints
		{
			UINT32 tmpM = 0;
			if (TableU[nlCons].size() > 0) {
				getTmpM(&tmpM, TableU[nlCons][0], red_pos, match_pos1, match_pos2);
			}
			else
			{
				continue;
			}
			
			// Build Table L
			map<UINT32, vector<UINT32>> TableL;
			TableL.clear();
			for (auto redV : TableU[nlCons]) // 2^3
			{
				//cout << redV << endl;
				tKeccakLane initial_red_state[25] = { 0 };
				SetBits(initial_red_state, red_pos, redV, 23);
				Permutation(initial_red_state, 2);
				UINT32 match_red = 0;
				getMatch(initial_red_state, &match_red, match_pos1, match_pos2);
				TableL[(match_red) & 0x7].push_back(redV);
			}

			// blue cell computation
			for (UINT32 blueV = 0; blueV < (1 << 3); blueV++) // 2^3
			{
				tKeccakLane initial_blue_state[25] = { 0 };
				SetBits(initial_blue_state, blue_pos, blueV, 3);
				SetBits(initial_blue_state, red_pos, TableU[nlCons][0], 23);
				Permutation(initial_blue_state, 2);
				UINT32 match_blue = 0;
				getMatch(initial_blue_state, &match_blue, match_pos1, match_pos2);
				for (auto redV_match : TableL[match_blue ^ tmpM]) {
					tKeccakLane initial_state[25] = { 0 };
					SetBits(initial_state, blue_pos, blueV, 3);
					SetBits(initial_state, red_pos, redV_match, 23);
					if (verify_3bit(initial_state)) {
						/*if (counter <= 10) {
							verify(initial_state, 1);
							cout << "Solution" << counter << endl;
							showState(initial_state);
						}*/
						counter += 1;
					}
					else
						error_counter += 1;
				}
			}
		}
	}
	// Attack end
	clock_t end_time = clock();
	cout << "\nThe time of finding all solutions satisfying the 3-bit partial target: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	cout << "\nThe number of all finding solutions: 2^" << log(double(counter)) / log(2.0) << " = 2^{8 + 12 + max{3, 3, 3}}" << endl;
	if (error_counter == 0)
		cout << "\nThe number of error solutions: " << error_counter << " , which means each pair passing the matching process is correct!" << endl;
}


int main() {
	keccak_attack_24bit();
	return 0;
}
