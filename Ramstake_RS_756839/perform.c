#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <stdint.h>

#include "ramstake.h"
#include "csprng.h"

uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

int main( int argc, char ** argv )
{
    unsigned long randomness;
    csprng rng;
    int i;
    int num_trials, trial_index;
    int num_successes;
    int trial_success;

    unsigned char * rng_seed;

    ramstake_public_key pk;
    ramstake_secret_key sk;
    ramstake_ciphertext c;

    unsigned char * sk_seed;
    unsigned char * c_seed;
    unsigned char * pk_bytes;
    unsigned char * sk_bytes;
    unsigned char * c_bytes;
    unsigned char * key1;
    unsigned char * key2;

    double time;

    struct timespec clock_total_start;
    struct timespec clock_total_stop;
    struct timespec clock_keygen_start;
    struct timespec clock_keygen_stop;
    struct timespec clock_encaps_start;
    struct timespec clock_encaps_stop;
    struct timespec clock_decaps_start;
    struct timespec clock_decaps_stop;

    uint64_t cycles_total_start;
    uint64_t cycles_total_stop;
    uint64_t cycles_keygen_start;
    uint64_t cycles_keygen_stop;
    uint64_t cycles_encaps_start;
    uint64_t cycles_encaps_stop;
    uint64_t cycles_decaps_start;
    uint64_t cycles_decaps_stop;

    if( argc != 3 || strlen(argv[2]) % 2 != 0 )
    {
        printf("usage: ./test [num trials, eg 13] [random seed, eg d13d13deadbeef]\n");
        printf("(And take note that time is measured in seconds so for a meaningful granularity set the number of trials to at least several hundred.)\n");
        return 0;
    }

    /* grab randomness */
    csprng_init(&rng);

    rng_seed = malloc(strlen(argv[2])/2);
    for( i = 0 ; i < strlen(argv[2])/2 ; ++i )
    {
        sscanf(argv[2] + 2*i, "%2hhx", &rng_seed[i]);
    }
    csprng_seed(&rng, strlen(argv[2])/2, rng_seed);
    free(rng_seed);

    randomness = csprng_generate_ulong(&rng);

    printf("randomness: %lu\n", randomness);

    /* grab trial number */
    num_trials = atoi(argv[1]);
    printf("num trials: %i\n", num_trials);

    /* allocate memory */
    sk_seed = malloc(RAMSTAKE_SEED_LENGTH * num_trials);
    csprng_generate(&rng, RAMSTAKE_SEED_LENGTH * num_trials, sk_seed);
    c_seed = malloc(RAMSTAKE_SEED_LENGTH * num_trials);
    csprng_generate(&rng, RAMSTAKE_SEED_LENGTH * num_trials, c_seed);
    pk_bytes = malloc(RAMSTAKE_PUBLIC_KEY_LENGTH * num_trials);
    sk_bytes = malloc(RAMSTAKE_SECRET_KEY_LENGTH * num_trials);
    c_bytes = malloc(RAMSTAKE_CIPHERTEXT_LENGTH * num_trials);
    key1 = malloc(RAMSTAKE_KEY_LENGTH * num_trials);
    key2 = malloc(RAMSTAKE_KEY_LENGTH * num_trials);

    /* run trials */
    num_successes = 0;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_total_start); cycles_total_start = rdtsc();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_keygen_start); cycles_keygen_start = rdtsc();
    for( trial_index = 0 ; trial_index < num_trials ; ++trial_index )
    {
        ramstake_public_key_init(&pk);
        ramstake_secret_key_init(&sk);

        ramstake_keygen(&sk, &pk, sk_seed + trial_index*RAMSTAKE_SEED_LENGTH, 0);

        ramstake_export_public_key(pk_bytes + trial_index*RAMSTAKE_PUBLIC_KEY_LENGTH, pk);
        ramstake_export_secret_key(sk_bytes + trial_index*RAMSTAKE_SECRET_KEY_LENGTH, sk);

        ramstake_public_key_destroy(pk);
        ramstake_secret_key_destroy(sk);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_keygen_stop); cycles_keygen_stop = rdtsc();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_encaps_start); cycles_encaps_start = rdtsc();
    for( trial_index = 0 ; trial_index < num_trials ; ++trial_index )
    {
        ramstake_public_key_init(&pk);
        ramstake_ciphertext_init(&c);

        ramstake_import_public_key(&pk, pk_bytes + trial_index*RAMSTAKE_PUBLIC_KEY_LENGTH);

        ramstake_encaps(&c, key1 + trial_index*RAMSTAKE_KEY_LENGTH, pk, c_seed + trial_index*RAMSTAKE_SEED_LENGTH, 0);

        ramstake_export_ciphertext(c_bytes + trial_index*RAMSTAKE_CIPHERTEXT_LENGTH, c);

        ramstake_public_key_destroy(pk);
        ramstake_ciphertext_destroy(c);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_encaps_stop); cycles_encaps_stop = rdtsc();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_decaps_start); cycles_decaps_start = rdtsc();
    for( trial_index = 0 ; trial_index < num_trials ; ++trial_index )
    {
        ramstake_secret_key_init(&sk);
        ramstake_ciphertext_init(&c);

        ramstake_import_ciphertext(&c, c_bytes + trial_index*RAMSTAKE_CIPHERTEXT_LENGTH);
        ramstake_import_secret_key(&sk, sk_bytes + trial_index*RAMSTAKE_SECRET_KEY_LENGTH);

        trial_success = ramstake_decaps(key2 + trial_index*RAMSTAKE_KEY_LENGTH, c, sk, 0);

        if( trial_success < 0 )
        {
            printf("failure found at trial index %i.\n", trial_index);
            printf("sk seed: ");
            for( i = 0 ; i < RAMSTAKE_SEED_LENGTH ; ++i )
                printf("%02x", sk_seed[trial_index*RAMSTAKE_SEED_LENGTH + i]);
            printf("\n");
            printf("c seed: ");
            for( i = 0 ; i < RAMSTAKE_SEED_LENGTH ; ++i )
                printf("%02x", c_seed[trial_index*RAMSTAKE_SEED_LENGTH + i]);
            printf("\n");
            break;
        }

        num_successes += (0 == strncmp((const char*)key1 + trial_index*RAMSTAKE_KEY_LENGTH, (const char*)key2 + trial_index*RAMSTAKE_KEY_LENGTH, RAMSTAKE_KEY_LENGTH));

        ramstake_secret_key_destroy(sk);
        ramstake_ciphertext_destroy(c);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_total_stop); cycles_total_stop = rdtsc();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &clock_decaps_stop); cycles_decaps_stop = rdtsc();

    /* free memory */
    free(sk_seed);
    free(c_seed);
    free(sk_bytes);
    free(pk_bytes);
    free(c_bytes);
    free(key1);
    free(key2);

    /* report on results */
    printf("Success rate: %i/%i.\n", num_successes, num_trials);

    time = 1000.0*(clock_total_stop.tv_sec - clock_total_start.tv_sec);
    printf("Total time: %f ms / %i trials = %f ms\n", time, num_trials, time/num_trials);
    printf("Total cycles: %lu / %i trials = %f\n", cycles_total_stop - cycles_total_start, num_trials, 1.0*(cycles_total_stop - cycles_total_start)/num_trials);

    time = 1000.0*(clock_keygen_stop.tv_sec - clock_keygen_start.tv_sec);
    printf("Keygen time: %f ms / %i trials = %f ms\n", time, num_trials, time/num_trials);
    printf("Keygen cycles: %lu / %i trials = %f\n", cycles_keygen_stop - cycles_keygen_start, num_trials, 1.0*(cycles_keygen_stop - cycles_keygen_start)/num_trials);

    time = 1000.0*(clock_encaps_stop.tv_sec - clock_encaps_start.tv_sec);
    printf("Encaps time: %f ms / %i trials = %f ms\n", time, num_trials, time/num_trials);
    printf("Encaps cycles: %lu / %i trials = %f\n", cycles_encaps_stop - cycles_encaps_start, num_trials, 1.0*(cycles_encaps_stop - cycles_encaps_start)/num_trials);

    time = 1000.0*(clock_decaps_stop.tv_sec - clock_decaps_start.tv_sec);
    printf("Decaps time: %f ms / %i trials = %f ms\n", time, num_trials, time/num_trials);
    printf("Decaps cycles: %lu / %i trials = %f\n", cycles_decaps_stop - cycles_decaps_start, num_trials, 1.0*(cycles_decaps_stop - cycles_decaps_start)/num_trials);

    return 0;
}

