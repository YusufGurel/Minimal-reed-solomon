#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * ----------FLEXIBLE TELECOMMAND DECODER (FTCD)----------
 *
 *  Authors:    Wouter Doppenberg & David Rijlaarsdam
 *  Date:       07-2019
 *  Course:     AE4S15 Space Embedded Systems - TU Delft
 *  Instructor: Alessandra Menicucci
 *
 * -------------------------------------------------------
 */

// CUSTOM HEADERS
#include "rs.h"
#include <string.h>

packet_t pack;
rs_t rst;
int main()
{

    init_rs(&rst, 255, 55);

    memset((uint8_t*)&pack, 0, sizeof(packet_t));



    for(int i = 0; i < 145; i++) {
        pack.data[i] = i % 255;
    }

    encode_rs(&pack, &rst);

    printf("ORIGINAL MESSAGE\n");
    for(int i = 0; i < 255; i++)
        printf("%d,", pack.data[i]);
    printf("\n");

    for(int i = 0; i < 55; i++)
        pack.data[10 + i] = 0;

    printf("CORRUPTED MESSAGE\n");
    for(int i = 0; i < 255; i++)
        printf("%d,", pack.data[i]);
    printf("\n");
    // write_to_file(MSG_SIZE, transfer.packs[0].data, "CORRUPTED_MESSAGE.csv");

    /*
     * And here.
     */

    decode_rs(&pack, &rst);

    printf("DECODED MESSAGE\n");
    for(int i = 0; i < 255; i++)
        printf("%d,", pack.data[i]);
    printf("\n");

    return 0;
}