#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <cmath>
#include <vector>
#include <chrono>

/* =====  Falcon-Lazy includes  ===== */
extern "C" {
#include "falcon-lazy/falcon.h"
}

/* =====  Dilithium-pqm4 includes  ===== */
extern "C" {
#include "dilithium-pqm4/api.h"
}

/* =====  Ed25519 (original tweetnacl fork) ===== */
extern "C" {
#include "ed25519/src/ed25519.h"
}

#if defined(__x86_64__) || defined(_M_X64)
#  define HAS_RDTSC 1
#  include <x86intrin.h>
static inline uint64_t rdtsc() {
    unsigned aux;
    return __rdtscp(&aux);          // serialising read
}
#else
#  define HAS_RDTSC 0
static inline uint64_t rdtsc() { return 0; }
#endif

/* ------------------------------------------------------------------ */
struct Stats {
    unsigned long long min_cycle = ~0ULL, max_cycle = 0;
    long double sum_cycle = 0, sumsq_cycle = 0;
    long double sum_us = 0;
    size_t n = 0;

    void add(uint64_t cycles, long double us) {
        min_cycle = std::min(min_cycle, cycles);
        max_cycle = std::max(max_cycle, cycles);
        sum_cycle += cycles;
        sumsq_cycle += (long double)cycles * cycles;
        sum_us += us;
        ++n;
    }
    void report(const char *label) const {
        long double avg_c = sum_cycle / n;
        long double var_c = (sumsq_cycle / n) - (avg_c * avg_c);
        printf("%s:\n", label);
        printf("    runs           : %zu\n", n);
        printf("    avg cycles     : %.0Lf\n", avg_c);
        printf("    min / max      : %llu / %llu\n", min_cycle, max_cycle);
        printf("    stdev cycles   : %.1Lf\n", std::sqrt(var_c));
        printf("    avg time (us)  : %.1Lf\n\n", sum_us / n);
    }
};

static inline uint64_t clock_now_us() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return ts.tv_sec * 1000000ULL + ts.tv_nsec / 1000ULL;
}

/* ------------------------------------------------------------------ */
int main() {
    constexpr size_t BENCH_ROUNDS = 1000;

    /* -------------  Falcon-Lazy ------------- */
    {
        unsigned logn = 9;                // 512-byte parameters
        shake256_context rng;
        char seed[16] = {0};
        shake256_init_prng_from_seed(&rng, seed, sizeof(seed));

        size_t pk_len = FALCON_PUBKEY_SIZE(logn);
        size_t sk_len = FALCON_PRIVKEY_SIZE(logn);
        size_t sig_len_ct = FALCON_SIG_CT_SIZE(logn);

        std::vector<uint8_t> pk(pk_len), sk(sk_len), sig(sig_len_ct);
        std::vector<uint8_t> tmp_keygen(FALCON_TMPSIZE_KEYGEN(logn));
        std::vector<uint8_t> tmp_sign(FALCON_TMPSIZE_SIGNDYN(logn));
        std::vector<uint8_t> tmp_verify(FALCON_TMPSIZE_VERIFY(logn));

        Stats st_keygen, st_sign, st_verify;
        const char msg[] = "data1";
        for (size_t i = 0; i < BENCH_ROUNDS; ++i) {
            /* Key-gen */
            uint64_t t0_us = clock_now_us();
            uint64_t c0 = rdtsc();
            falcon_keygen_make(&rng, logn, sk.data(), sk_len,
                               pk.data(), pk_len,
                               tmp_keygen.data(), tmp_keygen.size());
            uint64_t t1_us = clock_now_us();
            uint64_t c1 = rdtsc();
            st_keygen.add(c1 - c0, t1_us - t0_us);

            /* Sign (dynamic lazy) */
            t0_us = clock_now_us();
            c0 = rdtsc();
            falcon_sign_dyn_lazy(&rng,
                                 sig.data(), &sig_len_ct, FALCON_SIG_CT,
                                 pk.data(), pk_len,
                                 sk.data(), sk_len,
                                 msg, sizeof(msg) - 1,
                                 tmp_sign.data(), tmp_sign.size());
            t1_us = clock_now_us();
            c1 = rdtsc();
            st_sign.add(c1 - c0, t1_us - t0_us);

            /* Verify */
            t0_us = clock_now_us();
            c0 = rdtsc();
            falcon_verify(sig.data(), sig_len_ct, FALCON_SIG_CT,
                          pk.data(), pk_len,
                          msg, sizeof(msg) - 1,
                          tmp_verify.data(), tmp_verify.size());
            t1_us = clock_now_us();
            c1 = rdtsc();
            st_verify.add(c1 - c0, t1_us - t0_us);
        }
        printf("========== Falcon-Lazy ==========\n");
        st_keygen.report("Key generation");
        st_sign.report("Sign (dyn-lazy)");
        st_verify.report("Verify");
    }

    /* -------------  Dilithium (config.h selects mode) ------------- */
    {
        Stats st_keygen, st_sign, st_verify;
        constexpr size_t MLEN = 32;
        uint8_t m[MLEN] = {0}, sm[MLEN + CRYPTO_BYTES], m2[MLEN + CRYPTO_BYTES];
        size_t smlen, mlen2;
        uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];

        for (size_t i = 0; i < BENCH_ROUNDS; ++i) {
            /* keypair */
            uint64_t t0_us = clock_now_us();
            uint64_t c0 = rdtsc();
            crypto_sign_keypair(pk, sk);
            uint64_t t1_us = clock_now_us();
            uint64_t c1 = rdtsc();
            st_keygen.add(c1 - c0, t1_us - t0_us);

            /* sign */
            t0_us = clock_now_us();
            c0 = rdtsc();
            crypto_sign(sm, &smlen, m, MLEN, sk);
            t1_us = clock_now_us();
            c1 = rdtsc();
            st_sign.add(c1 - c0, t1_us - t0_us);

            /* verify (open) */
            t0_us = clock_now_us();
            c0 = rdtsc();
            crypto_sign_open(m2, &mlen2, sm, smlen, pk);
            t1_us = clock_now_us();
            c1 = rdtsc();
            st_verify.add(c1 - c0, t1_us - t0_us);
        }
        printf("========== Dilithium (mode %d) ==========\n", DILITHIUM_MODE);
        st_keygen.report("Key generation");
        st_sign.report("Sign");
        st_verify.report("Open/verify");
    }

    /* -------------  Ed25519 ------------- */
    {
        Stats st_keygen, st_sign, st_verify;
        constexpr size_t MLEN = 32;
        unsigned char pk[32], sk[64], seed[32], sig[64], msg[MLEN] = {0};

        for (size_t i = 0; i < BENCH_ROUNDS; ++i) {
            ed25519_create_seed(seed);
            /* keygen */
            uint64_t t0_us = clock_now_us();
            uint64_t c0 = rdtsc();
            ed25519_create_keypair(pk, sk, seed);
            uint64_t t1_us = clock_now_us();
            uint64_t c1 = rdtsc();
            st_keygen.add(c1 - c0, t1_us - t0_us);

            /* sign */
            t0_us = clock_now_us();
            c0 = rdtsc();
            ed25519_sign(sig, msg, MLEN, pk, sk);
            t1_us = clock_now_us();
            c1 = rdtsc();
            st_sign.add(c1 - c0, t1_us - t0_us);

            /* verify */
            t0_us = clock_now_us();
            c0 = rdtsc();
            ed25519_verify(sig, msg, MLEN, pk);
            t1_us = clock_now_us();
            c1 = rdtsc();
            st_verify.add(c1 - c0, t1_us - t0_us);
        }
        printf("========== Ed25519 ==========\n");
        st_keygen.report("Key generation");
        st_sign.report("Sign");
        st_verify.report("Verify");
    }

    return 0;
}
