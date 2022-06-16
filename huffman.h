#include "typedef.h"

typedef struct {
  int n;
  int rootIdx;
  int *left, *right;
  byte *clen;
  u64 *code;
  int *tbl;
  int tbl_width;
} Huffman;

void freeHuffman(Huffman *p);
int decodeFromHuffman(Huffman *h, u64 x);
Huffman *makeHuffmanTree(int n, i64 *freq, int tbl_width,int mode);
i64 Huffman_usedmemory(Huffman *h);
