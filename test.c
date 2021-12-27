#include "header.h"
#include <stdio.h>

// G, pho, pho1, phoprime 최적화 방안 연구
// 이후에 정리해서 보내기
// 전처리기 사용 고려

const BYTE GIFT_RC[48] = {
	0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3e, 0x3d, 0x3b, 0x37, 0x2f, 0x1e, 0x3c, 0x39, 0x33, 0x27, 0x0e,
	0x1d, 0x3a, 0x35, 0x2b, 0x16, 0x2c, 0x18, 0x30, 0x21, 0x02, 0x05, 0x0b, 0x17, 0x2e, 0x1c, 0x38,
	0x31, 0x23, 0x06, 0x0d, 0x1b, 0x36, 0x2d, 0x1a, 0x34, 0x29, 0x12, 0x24, 0x08, 0x11, 0x22, 0x04
};

//void GIFT_GetRoundConstant() {
//	int initial = 0b000000;
//	for (int i = 1; i < 49; i++) {
//		initial = AFFINE_LFSR(initial);
//		printf("0x%02x, ", initial);
//		if ((i % 16) == 0) printf("\n");
//	}
//}

static void padding(block d, block s, unsigned no_of_bytes) {
    block tmp;
    if (no_of_bytes == 0) {
        memset(tmp, 0, 16);
        tmp[0] = 0x80;
    }
    else if (no_of_bytes < 16) {
        memcpy(tmp, s, no_of_bytes);
        tmp[no_of_bytes] = 0x80;
        memset(tmp + no_of_bytes + 1, 0, 16 - no_of_bytes - 1);
    }
    else {
        memcpy(tmp, s, 16);
    }
    memcpy(d, tmp, 16);
}
/* ------------------------------------------------------------------------- */
static void xor_topbar_block(block d, block s1, half_block s2) {
    unsigned i;
    block tmp;
    XORBLOCK(tmp, s1, s2, 8);
    memcpy(tmp + 8, s1 + 8, 8);
    memcpy(d, tmp, 16);
}
/* ------------------------------------------------------------------------- */
static void double_half_block(half_block d, half_block s) {
    unsigned i;
    half_block tmp;
    /*x^{64} + x^4 + x^3 + x + 1*/
    for (i = 0; i < 7; i++)
        tmp[i] = (s[i] << 1) | (s[i + 1] >> 7);
    tmp[7] = (s[7] << 1) ^ ((s[0] >> 7) * 27);

    memcpy(d, tmp, 8);
}

static void triple_half_block(half_block d, half_block s) {
    unsigned i;
    half_block tmp;
    double_half_block(tmp, s);
    XORBLOCK(d, s, tmp, 8);
}
/* ------------------------------------------------------------------------- */

static void G(block d, block s) {
    unsigned i;
    block tmp;
    /*Y[1],Y[2] -> Y[2],Y[1]<<<1*/
    for (i = 0; i < 8; i++) {
        tmp[i] = s[8 + i];
    }
    for (i = 0; i < 7; i++) {
        tmp[i + 8] = s[i] << 1 | s[i + 1] >> 7;
    }
    tmp[7 + 8] = s[7] << 1 | s[0] >> 7;

    for (i = 0; i < 16; i++)
        d[i] = tmp[i];
}

static void pho1(block d, block Y, block M, int no_of_bytes) {
    block tmpM;
    half_block tmp; // tmp를 16-byte로 선언할 필요 없음.
    // G 함수
    for (int i = 0; i < 7; i++) {
        tmp[i] = Y[i] << 1 | Y[i + 1] >> 7;
    }
    tmp[7] = Y[7] << 1 | Y[0] >> 7;

    memcpy(Y, Y + 8, 8);
    memcpy(Y + 8, tmp, 8);

    padding(tmpM, M, no_of_bytes);
    XORBLOCK(d, Y, tmpM, 16);
}

static void pho(block Y, block M, block X, block C, int no_of_bytes) {
    XORBLOCK(C, Y, M, no_of_bytes);
    pho1(X, Y, M, no_of_bytes);
}

static void phoprime(block Y, block C, block X, block M, int no_of_bytes) {
    XORBLOCK(M, Y, C, no_of_bytes);
    pho1(X, Y, M, no_of_bytes);
}

static uint32_t rowperm(uint32_t S, int B0_pos, int B1_pos, int B2_pos, int B3_pos) {
    uint32_t T = 0;
    int b;
    for (b = 0; b < 8; b++) {
        T |= ((S >> (4 * b + 0)) & 0x1) << (b + 8 * B0_pos);
        T |= ((S >> (4 * b + 1)) & 0x1) << (b + 8 * B1_pos);
        T |= ((S >> (4 * b + 2)) & 0x1) << (b + 8 * B2_pos);
        T |= ((S >> (4 * b + 3)) & 0x1) << (b + 8 * B3_pos);
    }
    return T;
}

void GIFT128_Encrypt(uint8_t* P, const uint8_t* K, uint8_t* C) {
	uint32_t S[4] = { 0, }, T = 0;
	uint16_t W[8] = { 0, }, tmp0, tmp1;

	//		Initialization
	
	// b124 · · · b8  b4 b0	  S0
	// b125 · · · b9  b5 b1	  S1
	// b126 · · · b10 b6 b2 -> S2  =  S
	// b127 · · · b11 b7 b3	  S3
	S[0] = ((uint32_t)P[0] << 24) ^ ((uint32_t)P[1] << 16) ^ ((uint32_t)P[2] << 8) ^ (uint32_t)P[3];
	S[1] = ((uint32_t)P[4] << 24) ^ ((uint32_t)P[5] << 16) ^ ((uint32_t)P[6] << 8) ^ (uint32_t)P[7];
	S[2] = ((uint32_t)P[8] << 24) ^ ((uint32_t)P[9] << 16) ^ ((uint32_t)P[10] << 8) ^ (uint32_t)P[11];
	S[3] = ((uint32_t)P[12] << 24) ^ ((uint32_t)P[13] << 16) ^ ((uint32_t)P[14] << 8) ^ (uint32_t)P[15];

	//	Key State 
	W[0] = ((uint16_t)K[0] << 8) ^ (uint16_t)K[1];
	W[1] = ((uint16_t)K[2] << 8) ^ (uint16_t)K[3];
	W[2] = ((uint16_t)K[4] << 8) ^ (uint16_t)K[5];
	W[3] = ((uint16_t)K[6] << 8) ^ (uint16_t)K[7];
	W[4] = ((uint16_t)K[8] << 8) ^ (uint16_t)K[9];
	W[5] = ((uint16_t)K[10] << 8) ^ (uint16_t)K[11];
	W[6] = ((uint16_t)K[12] << 8) ^ (uint16_t)K[13];
	W[7] = ((uint16_t)K[14] << 8) ^ (uint16_t)K[15];

	// Round 
	for (int r = 0; r < 40; r++) {
		//		SubCells
		S[1] ^= S[0] & S[2];
		S[0] ^= S[1] & S[3];
		S[2] ^= S[0] | S[1];
		S[3] ^= S[2];
		S[1] ^= S[3];
		S[3]  = ~S[3];
		S[2] ^= S[0] & S[1];
		
		S[0] ^= S[3] ^= S[0] ^= S[3];	// Swap S[0] and S[1]

       //uint32_t tmp = S[0]; S[0] = S[1]; S[1] = tmp;

		//		PermBits
		T = 0;
		for (int i = 0; i < 8; i++) {
			T |= ((S[0] >> (4 * i) & 0x1) << i);
			T |= ((S[0] >> (4 * i + 1)) & 0x1) << (i + 8 * 3);
			T |= ((S[0] >> (4 * i + 2)) & 0x1) << (i + 8 * 2);
			T |= ((S[0] >> (4 * i + 3)) & 0x1) << (i + 8 * 1);
		}S[0] = T;

		T = 0;
		for (int i = 0; i < 8; i++) {
			T |= ((S[1] >> (4 * i) & 0x1) << (i + 8 * 1));
			T |= ((S[1] >> (4 * i + 1)) & 0x1) << i;
			T |= ((S[1] >> (4 * i + 2)) & 0x1) << (i + 8 * 3);
			T |= ((S[1] >> (4 * i + 3)) & 0x1) << (i + 8 * 2);
		}S[1] = T;

		T = 0;
		for (int i = 0; i < 8; i++) {
			T |= ((S[2] >> (4 * i) & 0x1) << (i + 8 * 2));
			T |= ((S[2] >> (4 * i + 1)) & 0x1) << (i + 8 * 1);
			T |= ((S[2] >> (4 * i + 2)) & 0x1) << i;
			T |= ((S[2] >> (4 * i + 3)) & 0x1) << (i + 8 * 3);
		}S[2] = T;

		T = 0;
		for (int i = 0; i < 8; i++) {
			T |= ((S[3] >> (4 * i) & 0x1) << (i + 8 * 3));
			T |= ((S[3] >> (4 * i + 1)) & 0x1) << (i + 8 * 2);
			T |= ((S[3] >> (4 * i + 2)) & 0x1) << (i + 8 * 1);
			T |= ((S[3] >> (4 * i + 3)) & 0x1) << i;
		}S[3] = T;

		//		AddRoundKey
		// Extract Round Key from Master Key 
		S[2] ^= ((uint32_t)W[2] << 16) | (uint32_t)W[3]; // U
		S[1] ^= ((uint32_t)W[6] << 16) | (uint32_t)W[7]; // V
		S[3] ^= (0x80000000 ^ GIFT_RC[r]);



		// Update Key State
		tmp0 = ROTR(W[6], 2);	tmp1 = ROTR(W[7], 12);	
		W[6] = W[4];	W[7] = W[5];
		W[4] = W[2];	W[5] = W[3];
		W[2] = W[0];	W[3] = W[1];
		W[0] = tmp0;	W[1] = tmp1;
	}

	//	Copy to Ciphertext
	C[ 0] = S[0] >> 24;	C[ 1] = S[0] >> 16;	C[ 2] = S[0] >> 8;	C[ 3] = S[0];
	C[ 4] = S[1] >> 24;	C[ 5] = S[1] >> 16;	C[ 6] = S[1] >> 8;	C[ 7] = S[1];
	C[ 8] = S[2] >> 24;	C[ 9] = S[2] >> 16;	C[10] = S[2] >> 8;	C[11] = S[2];
	C[12] = S[3] >> 24;	C[13] = S[3] >> 16;	C[14] = S[3] >> 8;	C[15] = S[3];
}

unsigned long long cpucycles() {
	return __rdtsc();
}
/* ------------------------------------------------------------------------- */

static int cofb_crypt(unsigned char* out, unsigned char* k, unsigned char* n,
    unsigned char* a, unsigned alen,
    unsigned char* in, unsigned inlen, int encrypting) {

    unsigned i;
    unsigned emptyA, emptyM;

    if (!encrypting) {
        if (inlen < TAGBYTES) return -1;
        inlen -= TAGBYTES;
    }

    if (alen == 0)
        emptyA = 1;
    else
        emptyA = 0;

    if (inlen == 0)
        emptyM = 1;
    else
        emptyM = 0;

    /*Mask-Gen*/
    block Y, input;
    half_block offset;
    /*nonce is 128-bit*/
    memcpy(input, n, 16);

    GIFT128_Encrypt(input, k, Y); 
    memcpy(offset, Y, 8);

    /*Process AD*/
    /*non-empty A*/
/*full blocks*/
    while (alen > 16) {
        /* X[i] = (A[i] + G(Y[i-1])) + offset */
        pho1(input, Y, a, 16);
        /* offset = 2*offset */
        double_half_block(offset, offset);
        xor_topbar_block(input, input, offset);
        /* Y[i] = E(X[i]) */
        GIFT128_Encrypt(input, k, Y);

        a = a + 16;
        alen -= 16;
    }

    /* last block */
    /* full block: offset = 3*offset */
    /* partial block: offset = 3^2*offset */
    triple_half_block(offset, offset);
    if ((alen % 16 != 0) || (emptyA)) {
        triple_half_block(offset, offset);
    }

    if (emptyM) {
        /* empty M: offset = 3^2*offset */
        triple_half_block(offset, offset);
        triple_half_block(offset, offset);
    }

    /* X[i] = (pad(A[i]) + G(Y[i-1])) + offset */
    pho1(input, Y, a, alen);

    xor_topbar_block(input, input, offset);
    /* Y[a] = E(X[a]) */
    GIFT128_Encrypt(input, k, Y);


    /* Process M */
    /* full blocks */
    while (inlen > 16) {
        double_half_block(offset, offset);
        /* C[i] = Y[i+a-1] + M[i]*/
        /* X[i] = M[i] + G(Y[i+a-1]) + offset */
        if (encrypting) {
            pho(Y, in, input, out, 16);
        }
        else {
            phoprime(Y, in, input, out, 16);
        }

        xor_topbar_block(input, input, offset);
        /* Y[i] = E(X[i+a]) */
        GIFT128_Encrypt(input, k, Y);

        in = in + 16;
        ////////////////////////////
        out = out + 16;
        ////////////////////////////
        inlen -= 16;
    }

    if (!emptyM) {
        /* full block: offset = 3*offset */
        /* empty data / partial block: offset = 3^2*offset */
        triple_half_block(offset, offset);
        if (inlen % 16 != 0) {
            triple_half_block(offset, offset);
        }
        /* last block */
        /* C[m] = Y[m+a-1] + M[m]*/
        /* X[a+m] = M[m] + G(Y[m+a-1]) + offset */
        if (encrypting) {
            pho(Y, in, input, out, inlen);
            out += inlen;
        }
        else {
            phoprime(Y, in, input, out, inlen);
            in += inlen;
        }

        xor_topbar_block(input, input, offset);
        /* T = E(X[m+a]) */
        GIFT128_Encrypt(input, k, Y);
    }

    if (encrypting) {
        memcpy(out, Y, TAGBYTES);
        return 0;
    }
    else
        return (memcmp(in, Y, TAGBYTES) ? -1 : 0);     /* Check for validity */
}

/* ------------------------------------------------------------------------- */

void cofb_encrypt(unsigned char* c, unsigned char* k, unsigned char* n,
    unsigned char* a, unsigned abytes,
    unsigned char* p, unsigned pbytes) 
{
    cofb_crypt(c, k, n, a, abytes, p, pbytes, COFB_ENCRYPT);
}
/* ------------------------------------------------------------------------- */
int cofb_decrypt(unsigned char* p, unsigned char* k, unsigned char* n,
    unsigned char* a, unsigned abytes,
    unsigned char* c, unsigned cbytes) {
    return cofb_crypt(p, k, n, a, abytes, c, cbytes, COFB_DECRYPT);
}

/* ------------------------------------------------------------------------- */

int crypto_aead_encrypt(
    unsigned char* c, unsigned long long* clen,         // 암호문, 암호문 길이
    const unsigned char* m, unsigned long long mlen,    // 메세지, 메세지 길이
    const unsigned char* ad, unsigned long long adlen,  // AD,     AD 길이
    const unsigned char* nsec,                          // 사용 안하는 변수
    const unsigned char* npub,                          // nonce로 들어감
    const unsigned char* k                              // 키
){
    (void)nsec;
    *clen = mlen + TAGBYTES;
    cofb_crypt(c, (unsigned char*)k, (unsigned char*)npub, (unsigned char*)ad,
        adlen, (unsigned char*)m, mlen, COFB_ENCRYPT);
    return 0;
}

int crypto_aead_decrypt(
    unsigned char* m, unsigned long long* mlen,
    unsigned char* nsec,
    const unsigned char* c, unsigned long long clen,
    const unsigned char* ad, unsigned long long adlen,
    const unsigned char* npub,
    const unsigned char* k
)
{
    (void)nsec;
    *mlen = clen - TAGBYTES;
    return cofb_crypt(m, (unsigned char*)k, (unsigned char*)npub,
        (unsigned char*)ad, adlen, (unsigned char*)c, clen, COFB_DECRYPT);
}


/*********************************************************************************************************/
/*******************************************OPEN SOURCE***************************************************/
/*********************************************************************************************************/
void giftb128(uint8_t P[16], const uint8_t K[16], uint8_t C[16]) {
    int round;
    uint32_t S[4], T;
    uint16_t W[8], T6, T7;

    S[0] = ((uint32_t)P[0] << 24) | ((uint32_t)P[1] << 16) | ((uint32_t)P[2] << 8) | (uint32_t)P[3];
    S[1] = ((uint32_t)P[4] << 24) | ((uint32_t)P[5] << 16) | ((uint32_t)P[6] << 8) | (uint32_t)P[7];
    S[2] = ((uint32_t)P[8] << 24) | ((uint32_t)P[9] << 16) | ((uint32_t)P[10] << 8) | (uint32_t)P[11];
    S[3] = ((uint32_t)P[12] << 24) | ((uint32_t)P[13] << 16) | ((uint32_t)P[14] << 8) | (uint32_t)P[15];

    W[0] = ((uint16_t)K[0] << 8) | (uint16_t)K[1];
    W[1] = ((uint16_t)K[2] << 8) | (uint16_t)K[3];
    W[2] = ((uint16_t)K[4] << 8) | (uint16_t)K[5];
    W[3] = ((uint16_t)K[6] << 8) | (uint16_t)K[7];
    W[4] = ((uint16_t)K[8] << 8) | (uint16_t)K[9];
    W[5] = ((uint16_t)K[10] << 8) | (uint16_t)K[11];
    W[6] = ((uint16_t)K[12] << 8) | (uint16_t)K[13];
    W[7] = ((uint16_t)K[14] << 8) | (uint16_t)K[15];

    for (round = 0; round < 40; round++) {
        /*===SubCells===*/
        S[1] ^= S[0] & S[2];
        S[0] ^= S[1] & S[3];
        S[2] ^= S[0] | S[1];
        S[3] ^= S[2];
        S[1] ^= S[3];
        S[3] ^= 0xffffffff;
        S[2] ^= S[0] & S[1];

        T = S[0];
        S[0] = S[3];
        S[3] = T;


        /*===PermBits===*/
        S[0] = rowperm(S[0], 0, 3, 2, 1);
        S[1] = rowperm(S[1], 1, 0, 3, 2);
        S[2] = rowperm(S[2], 2, 1, 0, 3);
        S[3] = rowperm(S[3], 3, 2, 1, 0);

        /*===AddRoundKey===*/
        S[2] ^= ((uint32_t)W[2] << 16) | (uint32_t)W[3];
        S[1] ^= ((uint32_t)W[6] << 16) | (uint32_t)W[7];

        /*Add round constant*/
        S[3] ^= 0x80000000 ^ GIFT_RC[round];

        /*===Key state update===*/
        T6 = (W[6] >> 2) | (W[6] << 14);
        T7 = (W[7] >> 12) | (W[7] << 4);
        W[7] = W[5];
        W[6] = W[4];
        W[5] = W[3];
        W[4] = W[2];
        W[3] = W[1];
        W[2] = W[0];
        W[1] = T7;
        W[0] = T6;
    }

    C[0] = S[0] >> 24;
    C[1] = S[0] >> 16;
    C[2] = S[0] >> 8;
    C[3] = S[0];
    C[4] = S[1] >> 24;
    C[5] = S[1] >> 16;
    C[6] = S[1] >> 8;
    C[7] = S[1];
    C[8] = S[2] >> 24;
    C[9] = S[2] >> 16;
    C[10] = S[2] >> 8;
    C[11] = S[2];
    C[12] = S[3] >> 24;
    C[13] = S[3] >> 16;
    C[14] = S[3] >> 8;
    C[15] = S[3];
}


static int cofb_crypt_open(unsigned char* out, unsigned char* k, unsigned char* n,
    unsigned char* a, unsigned alen,
    unsigned char* in, unsigned inlen, int encrypting) {

    unsigned i;
    unsigned emptyA, emptyM;

    if (!encrypting) {
        if (inlen < TAGBYTES) return -1;
        inlen -= TAGBYTES;
    }

    if (alen == 0)
        emptyA = 1;
    else
        emptyA = 0;

    if (inlen == 0)
        emptyM = 1;
    else
        emptyM = 0;

    /*Mask-Gen*/
    block Y, input;
    half_block offset;
    /*nonce is 128-bit*/
    for (i = 0; i < 16; i++)
        input[i] = n[i];

    giftb128(input, k, Y);
    for (i = 0; i < 8; i++)
        offset[i] = Y[i];


    /*Process AD*/
    /*non-empty A*/
/*full blocks*/
    while (alen > 16) {
        /* X[i] = (A[i] + G(Y[i-1])) + offset */
        pho1(input, Y, a, 16);
        /* offset = 2*offset */
        double_half_block(offset, offset);
        xor_topbar_block(input, input, offset);
        /* Y[i] = E(X[i]) */
        giftb128(input, k, Y);

        a = a + 16;
        alen -= 16;
    }

    /* last block */
    /* full block: offset = 3*offset */
    /* partial block: offset = 3^2*offset */
    triple_half_block(offset, offset);
    if ((alen % 16 != 0) || (emptyA)) {
        triple_half_block(offset, offset);
    }

    if (emptyM) {
        /* empty M: offset = 3^2*offset */
        triple_half_block(offset, offset);
        triple_half_block(offset, offset);
    }

    /* X[i] = (pad(A[i]) + G(Y[i-1])) + offset */
    pho1(input, Y, a, alen);

    xor_topbar_block(input, input, offset);
    /* Y[a] = E(X[a]) */
    giftb128(input, k, Y);


    /* Process M */
    /* full blocks */
    while (inlen > 16) {
        double_half_block(offset, offset);
        /* C[i] = Y[i+a-1] + M[i]*/
        /* X[i] = M[i] + G(Y[i+a-1]) + offset */
        if (encrypting) {
            pho(Y, in, input, out, 16);
        }
        else {
            phoprime(Y, in, input, out, 16);
        }

        xor_topbar_block(input, input, offset);
        /* Y[i] = E(X[i+a]) */
        giftb128(input, k, Y);

        in = in + 16;
        ////////////////////////////
        out = out + 16;
        ////////////////////////////
        inlen -= 16;
    }

    if (!emptyM) {
        /* full block: offset = 3*offset */
        /* empty data / partial block: offset = 3^2*offset */
        triple_half_block(offset, offset);
        if (inlen % 16 != 0) {
            triple_half_block(offset, offset);
        }
        /* last block */
        /* C[m] = Y[m+a-1] + M[m]*/
        /* X[a+m] = M[m] + G(Y[m+a-1]) + offset */
        if (encrypting) {
            pho(Y, in, input, out, inlen);
            out += inlen;
        }
        else {
            phoprime(Y, in, input, out, inlen);
            in += inlen;
        }


        xor_topbar_block(input, input, offset);
        /* T = E(X[m+a]) */
        giftb128(input, k, Y);
    }

    if (encrypting) {
        memcpy(out, Y, TAGBYTES);
        return 0;
    }
    else
        return (memcmp(in, Y, TAGBYTES) ? -1 : 0);     /* Check for validity */
}

int crypto_aead_encrypt_open(
    unsigned char* c, unsigned long long* clen,         // 암호문, 암호문 길이
    const unsigned char* m, unsigned long long mlen,    // 메세지, 메세지 길이
    const unsigned char* ad, unsigned long long adlen,  // AD,     AD 길이
    const unsigned char* nsec,                          // 사용 안하는 변수
    const unsigned char* npub,                          // nonce로 들어감
    const unsigned char* k                              // 키
)
{
    (void)nsec;
    *clen = mlen + TAGBYTES;
    cofb_crypt_open(c, (unsigned char*)k, (unsigned char*)npub, (unsigned char*)ad,
        adlen, (unsigned char*)m, mlen, COFB_ENCRYPT);
    return 0;
}
/*********************************************************************************************************/
void test() {
    /*
        • An encryption key K ∈ {0, 1}128
        • A nonce N ∈ {0, 1}^128
        • Associated data and message A, M ∈ {0, 1}*
        • Ciphertext C ∈ {0, 1}^|M|
        • Tag T ∈ {0, 1}^128
     */
#define LOOPLEN 10000
#define ADLEN 128
#define PLAINLEN 5
#define CIPHERLEN PLAINLEN + TAGBYTES
    const uint8_t K[16] = { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                            0x08 ,0x09 ,0x0A ,0x0B ,0x0C ,0x0D ,0x0E ,0x0F };
    //uint8_t Plaintext[PLAINLEN] = { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08 ,0x09 ,0x0A ,0x0B ,0x0C ,0x0D ,0x0E ,0x0F, };
    uint8_t Plaintext[PLAINLEN] = { 0x00, };
    uint8_t Ciphertext[CIPHERLEN] = { 0x00, };
    const uint8_t AD[ADLEN] = { 0x00, };
    const uint8_t* nsec = NULL;
    uint8_t nonce[16] = { 0x00, };
    uint64_t clen = CIPHERLEN, mlen = PLAINLEN;
#if 0
    unsigned long long start = 0, end = 0;
    start = cpucycles();
    for (int i = 0; i < LOOPLEN; i++) {
        giftb128(Plaintext, K, Ciphertext); // 배포코드
    }
    end = cpucycles();
    printf("%lld\n", (end - start) / LOOPLEN);

    start = cpucycles();
    for (int i = 0; i < LOOPLEN; i++) {
        GIFT128_Encrypt(Plaintext, K, Ciphertext);
    }
    end = cpucycles();
    printf("%lld\n", (end - start) / LOOPLEN);
#endif
    unsigned long long start = 0, end = 0;
    start = cpucycles();
    for (int i = 0; i < LOOPLEN; i++) {
        crypto_aead_encrypt(Ciphertext, &clen,
            Plaintext, mlen,
            AD, 16 * 8,
            nsec,
            nonce, K);
    }
    end = cpucycles();
    printf("mine = %lld\n", (end - start) / LOOPLEN);
    ///////////////////////////////////////////////
    const uint8_t K2[16] = { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                        0x08 ,0x09 ,0x0A ,0x0B ,0x0C ,0x0D ,0x0E ,0x0F };
    //uint8_t Plaintext2[PLAINLEN] = { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08 ,0x09 ,0x0A ,0x0B ,0x0C ,0x0D ,0x0E ,0x0F, };
    uint8_t Plaintext2[PLAINLEN] = { 0x0, };
    uint8_t Ciphertext2[CIPHERLEN] = { 0x00, };
    const uint8_t AD2[ADLEN] = { 0x00, };
    const uint8_t* nsec2 = NULL;
    uint8_t nonce2[16] = { 0x00, };
    uint64_t clen2 = CIPHERLEN, mlen2 = PLAINLEN;

    start = cpucycles();
    for (int i = 0; i < LOOPLEN; i++)
    {
        crypto_aead_encrypt_open(Ciphertext2, &clen2,
            Plaintext2, mlen2,
            AD2, 16 * 8,
            nsec2,
            nonce2, K2);
    }
    end = cpucycles();
    printf("open source  = %lld\n", (end - start) / LOOPLEN);

    printf("\n");
    for (int i = 0; i < CIPHERLEN; i++) {
        printf("%02x ", Ciphertext[i]);
    }
    printf("\n");
    for (int i = 0; i < CIPHERLEN; i++) {
        printf("%02x ", Ciphertext2[i]);
    }
}

int main() {

    test();

	return 0;
}