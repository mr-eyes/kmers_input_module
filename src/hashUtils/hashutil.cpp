/*
 * =====================================================================================
 *
 *       Filename:  hashutil.cc
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  04/18/2016 04:49:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
 *   Edited by: Mohamed Abuelanin (mabuelanin@gmail.com) UC Davis
 *
 * =====================================================================================
 */

#include "hashUtils/hashutil.hpp"
#include "Utils/kmer.h"
#include <iostream>
#include <string.h>
#include <unordered_map>

using namespace std;


inline string str_canonical(const string& kmer) {
    auto kmer_rev = kmer;
    std::reverse(kmer_rev.begin(), kmer_rev.end());
    for (size_t j = 0; j < kmer_rev.length(); ++j) {
        if (kmer_rev[j] == 'A') kmer_rev[j] = 'T';
        else if (kmer_rev[j] == 'T') kmer_rev[j] = 'A';
        else if (kmer_rev[j] == 'C') kmer_rev[j] = 'G';
        else if (kmer_rev[j] == 'G') kmer_rev[j] = 'C';
    }
    return kmer < kmer_rev ? kmer : kmer_rev;
}

// MurmurHash3 was written by Austin Appleby, and is placed in the public domain.
void MurmurHash3_x64_128(const void *key, int len, unsigned int seed, void *out) {
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 16;
    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = 0x87c37b91114253d5;
    const uint64_t c2 = 0x4cf5ad432745937f;

    // Body
    const uint64_t *blocks = (const uint64_t *)(data);
    for (int i = 0; i < nblocks; i++) {
        uint64_t k1 = blocks[i * 2 + 0];
        uint64_t k2 = blocks[i * 2 + 1];

        k1 *= c1;
        k1 = (k1 << 31) | (k1 >> (64 - 31));
        k1 *= c2;
        h1 ^= k1;

        h1 = (h1 << 27) | (h1 >> (64 - 27));
        h1 += h2;
        h1 = h1 * 5 + 0x52dce729;

        k2 *= c2;
        k2 = (k2 << 33) | (k2 >> (64 - 33));
        k2 *= c1;
        h2 ^= k2;

        h2 = (h2 << 31) | (h2 >> (64 - 31));
        h2 += h1;
        h2 = h2 * 5 + 0x38495ab5;
    }

    // Tail
    const uint8_t *tail = (const uint8_t *)(data + nblocks * 16);
    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch (len & 15) {
    case 15: k2 ^= ((uint64_t)tail[14]) << 48;
    case 14: k2 ^= ((uint64_t)tail[13]) << 40;
    case 13: k2 ^= ((uint64_t)tail[12]) << 32;
    case 12: k2 ^= ((uint64_t)tail[11]) << 24;
    case 11: k2 ^= ((uint64_t)tail[10]) << 16;
    case 10: k2 ^= ((uint64_t)tail[9]) << 8;
    case 9: k2 ^= ((uint64_t)tail[8]) << 0;
        k2 *= c2;
        k2 = (k2 << 33) | (k2 >> (64 - 33));
        k2 *= c1;
        h2 ^= k2;

    case 8: k1 ^= ((uint64_t)tail[7]) << 56;
    case 7: k1 ^= ((uint64_t)tail[6]) << 48;
    case 6: k1 ^= ((uint64_t)tail[5]) << 40;
    case 5: k1 ^= ((uint64_t)tail[4]) << 32;
    case 4: k1 ^= ((uint64_t)tail[3]) << 24;
    case 3: k1 ^= ((uint64_t)tail[2]) << 16;
    case 2: k1 ^= ((uint64_t)tail[1]) << 8;
    case 1: k1 ^= ((uint64_t)tail[0]) << 0;
        k1 *= c1;
        k1 = (k1 << 31) | (k1 >> (64 - 31));
        k1 *= c2;
        h1 ^= k1;
    };

    // Finalization
    h1 ^= len;
    h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 ^= h1 >> 33;
    h1 *= 0xff51afd7ed558ccd;
    h1 ^= h1 >> 33;
    h1 *= 0xc4ceb9fe1a85ec53;
    h1 ^= h1 >> 33;

    h2 ^= h2 >> 33;
    h2 *= 0xff51afd7ed558ccd;
    h2 ^= h2 >> 33;
    h2 *= 0xc4ceb9fe1a85ec53;
    h2 ^= h2 >> 33;

    h1 += h2;
    h2 += h1;

    ((uint64_t *)out)[0] = h1;
    ((uint64_t *)out)[1] = h2;
}

uint64_t MumurHasher::hash(const string & kmer) {
    string canonical_kmer = str_canonical(kmer);
    const char *c = canonical_kmer.c_str();
    
    uint64_t hash_output[2];
    MurmurHash3_x64_128(canonical_kmer.data(), canonical_kmer.size(), seed, hash_output);

    int64_t hash = static_cast<int64_t>(hash_output[0]);
    if (hash < 0) hash += static_cast<uint64_t>(1) << 64;

    return static_cast<uint64_t>(hash);
  
}

IntegerHasher::IntegerHasher(uint64_t kSize) {
    this->kSize = kSize;
    this->mask = BITMASK(2 * kSize);
}

/*
 *   For any 1<k<=64, let mask=(1<<k)-1. hash_64() is a bijection on [0,1<<k),
 *   which means
 *     hash_64(x, mask)==hash_64(y, mask) if and only if x==y. hash_64i() is
 *     the inversion of
 *       hash_64(): hash_64i(hash_64(x, mask), mask) == hash_64(hash_64i(x,
 *       mask), mask) == x.
 */

// Thomas Wang's integer hash functions. See
// <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
uint64_t IntegerHasher::hash(const string &kmer) {
    uint64_t key = kmer::str_to_canonical_int(kmer);
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

uint64_t IntegerHasher::hash(uint64_t key) {
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

// The inversion of hash_64(). Modified from
// <https://naml.us/blog/tag/invertible>
string IntegerHasher::Ihash(uint64_t key) {
    uint64_t tmp;

    // Invert key = key + (key << 31)
    tmp = (key - (key << 31));
    key = (key - (tmp << 31)) & mask;

    // Invert key = key ^ (key >> 28)
    tmp = key ^ key >> 28;
    key = key ^ tmp >> 28;

    // Invert key *= 21
    key = (key * 14933078535860113213ull) & mask;

    // Invert key = key ^ (key >> 14)
    tmp = key ^ key >> 14;
    tmp = key ^ tmp >> 14;
    tmp = key ^ tmp >> 14;
    key = key ^ tmp >> 14;

    // Invert key *= 265
    key = (key * 15244667743933553977ull) & mask;

    // Invert key = key ^ (key >> 24)
    tmp = key ^ key >> 24;
    key = key ^ tmp >> 24;

    // Invert key = (~key) + (key << 21)
    tmp = ~key;
    tmp = ~(key - (tmp << 21));
    tmp = ~(key - (tmp << 21));
    key = ~(key - (tmp << 21)) & mask;

    return kmer::int_to_str(key, kSize);
}

// QHasher _________________________________-

QHasher::QHasher(uint64_t kSize) {
    this->kSize = kSize;
    this->mask = BITMASK(this->Q);
}

QHasher::QHasher(uint64_t kSize, int Q) {
    this->kSize = kSize;
    this->Q = Q;
    this->key_remainder_bits = (2 * kSize) - Q;
    this->mask = (1 << (Q)) - 1;
}

void QHasher::set_Q(int _Q) {
    this->Q = _Q;
    this->key_remainder_bits = (2 * kSize) - _Q;
}

uint64_t QHasher::merge_Q_R(uint64_t &Qval, uint64_t &R) {
    // cout << "merge_Q_R(" << Qval << "," << R << ") = ";
    // cout << ((Qval << this->key_remainder_bits) | R) << endl;
    return (Qval << this->key_remainder_bits) | R;
}

void QHasher::split_Q_R(uint64_t key, uint64_t &Qval, uint64_t &R) {

    R = key & BITMASK(this->key_remainder_bits);
    Qval = key >> this->key_remainder_bits;
    // cout << "Split_Q_R("<< key <<") = ";
    // cout << "Q: " << _Q << ", R: " << R << endl;
}

uint64_t QHasher::normal_hash(const string &kmer) {
    uint64_t key = kmer::str_to_canonical_int(kmer);
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

uint64_t QHasher::normal_hash(uint64_t key) {
    // cout << "normal_hash(" << key <<")";
    key = (~key + (key << 21U)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24U;
    key = ((key + (key << 3U)) + (key << 8U)) & mask; // key * 265
    key = key ^ key >> 14U;
    key = ((key + (key << 2U)) + (key << 4U)) & mask; // key * 21
    key = key ^ key >> 28U;
    key = (key + (key << 31U)) & mask;
    // cout << " = " << key << endl;
    return key;
}

uint64_t QHasher::normal_Ihash(uint64_t key) {
    // cout << "normal_Ihash (" << key << ") = " << endl;
    uint64_t tmp;

    // Invert key = key + (key << 31)
    tmp = (key - (key << 31U));
    key = (key - (tmp << 31U)) & mask;

    // Invert key = key ^ (key >> 28)
    tmp = key ^ key >> 28U;
    key = key ^ tmp >> 28U;

    // Invert key *= 21
    key = (key * 14933078535860113213ull) & mask;

    // Invert key = key ^ (key >> 14)
    tmp = key ^ key >> 14;
    tmp = key ^ tmp >> 14;
    tmp = key ^ tmp >> 14;
    key = key ^ tmp >> 14;

    // Invert key *= 265
    key = (key * 15244667743933553977ull) & mask;

    // Invert key = key ^ (key >> 24)
    tmp = key ^ key >> 24U;
    key = key ^ tmp >> 24U;

    // Invert key = (~key) + (key << 21)
    tmp = ~key;
    tmp = ~(key - (tmp << 21U));
    tmp = ~(key - (tmp << 21U));
    key = ~(key - (tmp << 21U)) & mask;

    // cout << key << endl;
    return key;
}


uint64_t QHasher::hash(const string &key) {
//	    // cout << "hash()" << endl;
    uint64_t newHash;
    uint64_t Qval;
    uint64_t R;
    uint64_t hashed_Q;
    uint64_t _2bit = kmer::str_to_canonical_int(key);
    // cout << "_2bit: " << _2bit << endl;
    split_Q_R(_2bit, Qval, R);
//    // cout << "splitting| Q: " << _Q << ", R: " << R << endl;
    hashed_Q = normal_hash(Qval);
    newHash = merge_Q_R(hashed_Q, R);
    return newHash;
}

uint64_t QHasher::hash(uint64_t key) {
    uint64_t newHash;
    uint64_t Qval;
    uint64_t R;
    uint64_t hashed_Q;
    split_Q_R(key, Qval, R);
    hashed_Q = normal_hash(Qval);
    newHash = merge_Q_R(hashed_Q, R);
    return newHash;
}

string QHasher::Ihash(uint64_t key) {
    uint64_t _2bit;
    uint64_t Qval;
    uint64_t R;
    uint64_t hashed_Q;
    split_Q_R(key, hashed_Q, R);
    Qval = normal_Ihash(hashed_Q);
    _2bit = merge_Q_R(Qval, R);
    return kmer::int_to_str(_2bit, kSize);
}


// _________ TwoBitsHasher

TwoBitsHasher::TwoBitsHasher(uint64_t kSize) {
    this->kSize = kSize;
}

uint64_t TwoBitsHasher::hash(const string &key) {
    return kmer::str_to_canonical_int(key);
}

uint64_t TwoBitsHasher::hash(uint64_t key) {
    return key;
}

string TwoBitsHasher::Ihash(uint64_t key) {
    return kmer::int_to_str(key, this->kSize);
}


// _________ TwoBitsHasher


uint64_t noncanonical_TwoBitsHasher::hash(const string &key) {
    return kmer::str_to_int(key);
}

// nonCanonical_IntegerHasher

uint64_t noncanonical_IntegerHasher::hash(uint64_t key) {
  key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

// The inversion of hash_64(). Modified from
// <https://naml.us/blog/tag/invertible>
string noncanonical_IntegerHasher::Ihash(uint64_t key) {
  uint64_t tmp;

  // Invert key = key + (key << 31)
  tmp = (key - (key << 31));
  key = (key - (tmp << 31)) & mask;

  // Invert key = key ^ (key >> 28)
  tmp = key ^ key >> 28;
  key = key ^ tmp >> 28;

  // Invert key *= 21
  key = (key * 14933078535860113213ull) & mask;

  // Invert key = key ^ (key >> 14)
  tmp = key ^ key >> 14;
  tmp = key ^ tmp >> 14;
  tmp = key ^ tmp >> 14;
  key = key ^ tmp >> 14;

  // Invert key *= 265
  key = (key * 15244667743933553977ull) & mask;

  // Invert key = key ^ (key >> 24)
  tmp = key ^ key >> 24;
  key = key ^ tmp >> 24;

  // Invert key = (~key) + (key << 21)
  tmp = ~key;
  tmp = ~(key - (tmp << 21));
  tmp = ~(key - (tmp << 21));
  key = ~(key - (tmp << 21)) & mask;

  return kmer::int_to_str(key, kSize);
}

uint64_t noncanonical_IntegerHasher::hash(const string &kmer) {
    uint64_t key = kmer::str_to_int(kmer);
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

// _________ bigKmerHasher

inline string bigKmerRevComplement(string DNAseq){
  reverse(DNAseq.begin(), DNAseq.end());
  for (char & i : DNAseq){
    switch (i){
      case 'A': i = 'T'; break;
      case 'C': i = 'G'; break;
      case 'G': i = 'C'; break;
      case 'T': i = 'A'; break;
    }
  }
  return DNAseq;
}

uint64_t bigKmerHasher::hash(uint64_t key) {
    return Hasher::hash(key);
}

bigKmerHasher::bigKmerHasher(uint64_t kSize) {
    this->kSize = kSize;
}

uint64_t bigKmerHasher::hash(const string &key) {
    return this->hasher(get_canonical_kmer(key));
}

string bigKmerHasher::get_canonical_kmer(const string &kmer) {
    string revComp = bigKmerRevComplement(kmer);
    return (kmer < revComp) ? kmer : revComp;
}
