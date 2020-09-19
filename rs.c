#include "rs.h"

short pp[9] = { 1, 0, 1, 1, 1, 0, 0, 0, 1 };
// short alpha_to[nn+1], index_of[nn+1], gg[nn-kk+1], recd[nn];

void generate_gf(rs_t* rs)
{
    register int i, mask;

    mask = 1;
    rs->alpha_to[mm] = 0;
    for(i = 0; i < mm; i++) {
        rs->alpha_to[i] = mask;
        rs->index_of[rs->alpha_to[i]] = i;
        if(pp[i] != 0)
            rs->alpha_to[mm] ^= mask;
        mask <<= 1;
    }
    rs->index_of[rs->alpha_to[mm]] = mm;
    mask >>= 1;
    for(i = mm + 1; i < rs->nn; i++) {
        if(rs->alpha_to[i - 1] >= mask)
            rs->alpha_to[i] = rs->alpha_to[mm] ^ ((rs->alpha_to[i - 1] ^ mask) << 1);
        else
            rs->alpha_to[i] = rs->alpha_to[i - 1] << 1;
        rs->index_of[rs->alpha_to[i]] = i;
    }
    rs->index_of[0] = -1;
}

void gen_poly(rs_t* rs)
{
    register int i, j;

    rs->gg[0] = 2;
    rs->gg[1] = 1;
    for(i = 2; i <= rs->nn - rs->kk; i++) {
        rs->gg[i] = 1;
        for(j = i - 1; j > 0; j--)
            if(rs->gg[j] != 0)
                rs->gg[j] = rs->gg[j - 1] ^ rs->alpha_to[(rs->index_of[rs->gg[j]] + i) % rs->nn];
            else
                rs->gg[j] = rs->gg[j - 1];
        rs->gg[0] = rs->alpha_to[(rs->index_of[rs->gg[0]] + i) % rs->nn];
    }

    for(i = 0; i <= rs->nn - rs->kk; i++)
        rs->gg[i] = rs->index_of[rs->gg[i]];
}

void init_rs(rs_t* rs, uint8_t lnn, uint8_t ltt)
{
    rs->nn = lnn;
    rs->tt = ltt;
    rs->kk = rs->nn - 2 * rs->tt;
    
    generate_gf(rs);
    gen_poly(rs);
}


void encode_rs(packet_t* packet, rs_t* rs)
{
    register int i, j;
    int feedback;

    for(i = 0; i < rs->nn - rs->kk; i++)
        packet->ECF[i] = 0;
    for(i = rs->kk - 1; i >= 0; i--) {
        feedback = rs->index_of[packet->data[i] ^ packet->ECF[rs->nn - rs->kk - 1]];
        if(feedback != -1) {
            for(j = rs->nn - rs->kk - 1; j > 0; j--) {
                if(rs->gg[j] != -1)
                    packet->ECF[j] = packet->ECF[j - 1] ^ rs->alpha_to[(rs->gg[j] + feedback) % rs->nn];
                else
                    packet->ECF[j] = packet->ECF[j - 1];
            }

            packet->ECF[0] = rs->alpha_to[(rs->gg[0] + feedback) % rs->nn];
        } else {
            for(j = rs->nn - rs->kk - 1; j > 0; j--)
                packet->ECF[j] = packet->ECF[j - 1];
            packet->ECF[0] = 0;
        }
    }
}
void decode_rs(packet_t* packet, rs_t* rs)
{

    register int i, j, u, q;

    for(i = 0; i < rs->nn - rs->kk; i++)
        rs->recd[i] = packet->ECF[i];
    for(i = 0; i < rs->kk; i++)
        rs->recd[i + rs->nn - rs->kk] = packet->data[i];

    for(i = 0; i < rs->nn; i++)
        rs->recd[i] = rs->index_of[rs->recd[i]];

    int elp[rs->nn - rs->kk + 2][rs->nn - rs->kk], d[rs->nn - rs->kk + 2], l[rs->nn - rs->kk + 2],
        u_lu[rs->nn - rs->kk + 2], s[rs->nn - rs->kk + 1];
    int count = 0, syn_error = 0, root[rs->tt], loc[rs->tt], z[rs->tt + 1], err[rs->nn], reg[rs->tt + 1];

    for(i = 1; i <= rs->nn - rs->kk; i++) {
        s[i] = 0;
        for(j = 0; j < rs->nn; j++)
            if(rs->recd[j] != -1)
                s[i] ^= rs->alpha_to[(rs->recd[j] + i * j) % rs->nn];

        if(s[i] != 0)
            syn_error = 1;
        s[i] = rs->index_of[s[i]];
    };

    if(syn_error) {

        d[0] = 0;
        d[1] = s[1];
        elp[0][0] = 0;
        elp[1][0] = 1;
        for(i = 1; i < rs->nn - rs->kk; i++) {
            elp[0][i] = -1;
            elp[1][i] = 0;
        }
        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;
        u = 0;

        do {
            u++;
            if(d[u] == -1) {
                l[u + 1] = l[u];
                for(i = 0; i <= l[u]; i++) {
                    elp[u + 1][i] = elp[u][i];
                    elp[u][i] = rs->index_of[elp[u][i]];
                }
            } else

            {
                q = u - 1;
                while((d[q] == -1) && (q > 0))
                    q--;

                if(q > 0) {
                    j = q;
                    do {
                        j--;
                        if((d[j] != -1) && (u_lu[q] < u_lu[j]))
                            q = j;
                    } while(j > 0);
                };

                if(l[u] > l[q] + u - q)
                    l[u + 1] = l[u];
                else
                    l[u + 1] = l[q] + u - q;

                for(i = 0; i < rs->nn - rs->kk; i++)
                    elp[u + 1][i] = 0;
                for(i = 0; i <= l[q]; i++)
                    if(elp[q][i] != -1)
                        elp[u + 1][i + u - q] = rs->alpha_to[(d[u] + rs->nn - d[q] + elp[q][i]) % rs->nn];
                for(i = 0; i <= l[u]; i++) {
                    elp[u + 1][i] ^= elp[u][i];
                    elp[u][i] = rs->index_of[elp[u][i]];
                }
            }
            u_lu[u + 1] = u - l[u + 1];

            if(u < rs->nn - rs->kk) {
                if(s[u + 1] != -1)
                    d[u + 1] = rs->alpha_to[s[u + 1]];
                else
                    d[u + 1] = 0;
                for(i = 1; i <= l[u + 1]; i++)
                    if((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
                        d[u + 1] ^= rs->alpha_to[(s[u + 1 - i] + rs->index_of[elp[u + 1][i]]) % rs->nn];
                d[u + 1] = rs->index_of[d[u + 1]];
            }
        } while((u < rs->nn - rs->kk) && (l[u + 1] <= rs->tt));

        u++;
        if(l[u] <= rs->tt) {

            for(i = 0; i <= l[u]; i++)
                elp[u][i] = rs->index_of[elp[u][i]];

            for(i = 1; i <= l[u]; i++)
                reg[i] = elp[u][i];
            count = 0;
            for(i = 1; i <= rs->nn; i++) {
                q = 1;
                for(j = 1; j <= l[u]; j++)
                    if(reg[j] != -1) {
                        reg[j] = (reg[j] + j) % rs->nn;
                        q ^= rs->alpha_to[reg[j]];
                    };
                if(!q) {
                    root[count] = i;
                    loc[count] = rs->nn - i;
                    count++;
                };
            };
            if(count == l[u]) {

                for(i = 1; i <= l[u]; i++) {
                    if((s[i] != -1) && (elp[u][i] != -1))
                        z[i] = rs->alpha_to[s[i]] ^ rs->alpha_to[elp[u][i]];
                    else if((s[i] != -1) && (elp[u][i] == -1))
                        z[i] = rs->alpha_to[s[i]];
                    else if((s[i] == -1) && (elp[u][i] != -1))
                        z[i] = rs->alpha_to[elp[u][i]];
                    else
                        z[i] = 0;
                    for(j = 1; j < i; j++)
                        if((s[j] != -1) && (elp[u][i - j] != -1))
                            z[i] ^= rs->alpha_to[(elp[u][i - j] + s[j]) % rs->nn];
                    z[i] = rs->index_of[z[i]];
                };

                for(i = 0; i < rs->nn; i++) {
                    err[i] = 0;
                    if(rs->recd[i] != -1)
                        rs->recd[i] = rs->alpha_to[rs->recd[i]];
                    else
                        rs->recd[i] = 0;
                }
                for(i = 0; i < l[u]; i++) {
                    err[loc[i]] = 1;
                    for(j = 1; j <= l[u]; j++)
                        if(z[j] != -1)
                            err[loc[i]] ^= rs->alpha_to[(z[j] + j * root[i]) % rs->nn];
                    if(err[loc[i]] != 0) {
                        err[loc[i]] = rs->index_of[err[loc[i]]];
                        q = 0;
                        for(j = 0; j < l[u]; j++)
                            if(j != i)
                                q += rs->index_of[1 ^ rs->alpha_to[(loc[j] + root[i]) % rs->nn]];
                        q = q % rs->nn;
                        err[loc[i]] = rs->alpha_to[(err[loc[i]] - q + rs->nn) % rs->nn];
                        rs->recd[loc[i]] ^= err[loc[i]];
                    }
                }
            } else
                for(i = 0; i < rs->nn; i++)
                    if(rs->recd[i] != -1)
                        rs->recd[i] = rs->alpha_to[rs->recd[i]];
                    else
                        rs->recd[i] = 0;
        } else
            for(i = 0; i < rs->nn; i++)
                if(rs->recd[i] != -1)
                    rs->recd[i] = rs->alpha_to[rs->recd[i]];
                else
                    rs->recd[i] = 0;
    } else
        for(i = 0; i < rs->nn; i++)
            if(rs->recd[i] != -1)
                rs->recd[i] = rs->alpha_to[rs->recd[i]];
            else
                rs->recd[i] = 0;

    for(i = 0; i < rs->nn - rs->kk; i++)
        packet->ECF[i] = rs->recd[i];
    for(i = 0; i < rs->kk; i++)
        packet->data[i] = rs->recd[i + rs->nn - rs->kk];
}


