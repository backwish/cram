#include "typedef.h"
typedef struct{
    int *requiredMinCodeLen;
} CodeAllocator;

CodeAllocator *makeCodeAllocator(int maxDepth);
void freeCodeAllocator(CodeAllocator *self);
u64 codeAllocate(CodeAllocator *self,int size);

typedef struct{
    int ch;
    int len;
} TableElement;
typedef struct{
    int maxDepth;
    TableElement *tbl;    
} AdditionalCodeTable;
AdditionalCodeTable *makeAdditionalCodeTable(int maxDepth);
void freeAdditionalCodeTable(AdditionalCodeTable *self);
void insertToAdditionalCodeTable(AdditionalCodeTable *self,u64 code,int ch,int len);
TableElement decodeFromAdditionalCodeTable(AdditionalCodeTable *self,u64 code);