#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/timeb.h>
#include "cram.h"
#include "mman.h"
#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}
typedef struct timeb mytimestruct;
void mygettime(mytimestruct *t)
{
  ftime(t);
}
double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->time - before->time;
  t += (double)(after->millitm - before->millitm)/1000;
  return t;
}
int correctnessCheck(CRAM* cram,uchar* dest,int n,int b){    
    int s = 0,j = 0;
    uchar *tmpbuf;
    printf("START CORRECTNESS CHECK\n");
    mymalloc(tmpbuf,b);    
    for(s=0;s<n;s+=b){
        cram_read_block(cram,s/b,tmpbuf);
        for(j=0;j<b;j++){
            if(tmpbuf[j]!=dest[s+j]){
                printf("CRAM correctness test failed %d %d %d \n",cram->step*b,s,j);
                free(tmpbuf);
                return 0;
            }
        }
    }
    free(tmpbuf);
    return 1;
}
void printHuffmanSize(CRAM *cram){
    uchar *tmpbuf;
    mymalloc(tmpbuf,1024);
    int i,ret;
    ret = 0;
    for(i=0;i<cram->nb;i++){
        ret+=cram_read_block(cram,i,tmpbuf);
    }
    printf("HUFFMAN SIZE bpc : %1.5f\n",(double)ret/cram->n);
    free(tmpbuf);
}
int main(int argc, char *argv[]){
    CRAM *cram;
    i64 s,t,i,j,n, n2, ss,ret;
    uchar *buf;
    MMAP *map;
    uchar *text2;
    int b, bu, mode,changeMultiple,testMode,testResult;    
    mytimestruct before,after;    
    double retBpc,answer;
    double outputOurs[4][4] = {{3.03,3.53,2.86,2.75},{6.24,6.24,2.98,2.75},{2.46,3.75,2.59,2.41},{6.12,6.12,2.68,2.41}};
    double outputSada[3] = {9.59,2.48,2.25};    
    b = 1024;    
    testMode = 0;bu = 4;mode = 2;changeMultiple = 2;        
    if (argc>=4){
        testMode = atoi(argv[3]);
        bu = atoi(argv[4]);
        mode = atoi(argv[5]);
        changeMultiple = atoi(argv[6]);
    }
    answer = mode==0 ? outputSada[bu/2] : outputOurs[(mode-1)*2+changeMultiple/4][bu==4 ? 3 : bu];

    printf("CRAM TEST START bu : %d, mode : %d, changeMultiple : %d\n",bu,mode,changeMultiple);    

    cram = cram_initialize(argv[1], b, bu, mode,changeMultiple);
    n = cram->n;    

    map = mymmap(argv[2]);
    text2 = map->addr;
    n2 = map->len;
    if (n2 > n) n2 = n;
    mymalloc(buf, 1024);

    
    //printHuffmanSize(cram);
    mygettime(&before);
        
    for(s=0;s<n2;s+=b){
        if(testMode==0 && s % ((n2/b/10)*b) == 0){
            fprintf(stderr, "%ld \r", s / ((n2/b/100)*b));  fflush(stderr);
            printf("%d %1.2f\n", s / ((n2/b/100)*b), (double)cram_usedmemory(cram)*8/n);
        }
        t = s/b;
        //t = s+b-1;  if (t >= n) t = n-1;
        cram_write(cram, t, &text2[s]);
        cram_read_block(cram,t,buf);
        for(j=0;j<b;j++){
            if(buf[j]!=text2[s+j]){
                printf("CRAM correctness test failed %d %d %d \n",cram->step*b,s,j);
                free(buf);
                return 0;
            }
        }
    }
    mygettime(&after);
    //printHuffmanSize(cram);
    ret = cram_usedmemory(cram);
    retBpc = (double)ret*8/n;
    if(answer - retBpc < -0.03 || answer - retBpc > 0.03){
        printf("CRAM INVALID IMPLEMENTATION\n");
        return 0;
    }
    if(testMode==0){
        fprintf(stderr, "rewrite %f sec b = %d\n", mylaptime(&before,&after), b);
        fprintf(stderr, "used memory %ld bytes (%1.2f bpc)\n", ret, retBpc);
        printf("%d %1.2f\n", 100, retBpc);        
    }    
    testResult = correctnessCheck(cram,text2,n2,b);
    if(testResult==1){
        printf("CRAM correctness test SUCCESS\n");
    }
    
}