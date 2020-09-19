#ifndef RS_H
#define RS_H

#include <stdint.h>

#define mm 8 /* RS code over GF(2**8)  */
//#define nn  255     /* nn=2**mm -1   length of codeword */
//#define tt  55      /* number of errors that can be corrected */
//#define kk  (nn-(2*tt))     /* kk = nn-2*tt  */



struct rs_s {


    uint8_t nn;
    uint8_t tt;
    uint8_t kk;
    uint8_t msg_size;
    short alpha_to[256];
    short index_of[256];
    short gg[256];
    short recd[256];
};
typedef struct rs_s rs_t;

struct packet {

    // uint8_t data[kk];
    // uint8_t ECF[nn-kk];
    uint8_t data[255];
    uint8_t ECF[255];
};
typedef struct packet packet_t;

void init_rs(rs_t* rs, uint8_t lnn, uint8_t ltt);
void encode_rs(packet_t* packet, rs_t* rs);
void decode_rs(packet_t* packet, rs_t* rs);
#endif // RS_H
