#pragma once
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <intrin.h>

#define TAGBYTES 16	// GIFT64 : 8 / GIFT128 : 16

typedef unsigned char block[16];
typedef unsigned char half_block[8];
typedef unsigned char BYTE;

#define ROTR(x, n) (((x)>>(n)) | ((x) << (16 - (n))))
#define ROTL(x, n) (((x)<<(n)) | ((x) >> (16 - (n))))
#define AFFINE_LFSR(x) ((((x) << 1) ^ (((x) >> 5) ^ (((x) >> 4) & 0b1) ^ 0b1)) & 0b111111)
#define XORBLOCK(d,s1,s2,n) for(int i = 0; i < n; i++) d[i] = s1[i] ^ s2[i]

#define COFB_ENCRYPT 1
#define COFB_DECRYPT 0


void GIFT128_Encrypt(uint8_t* P, const uint8_t* K, uint8_t* C);