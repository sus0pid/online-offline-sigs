#include "benchmark/benchmark.h"
#include "testlib.h"

static void lazy_offline(benchmark::State& state) {
    // Perform setup here
    const uint64_t logn = 9;
    const uint64_t n = 1 << logn;
    inner_shake256_context rng;
    inner_shake256_init(&rng);
    falcon_key_t key = keygen(logn, &rng);
    vec_modQ hq = to_vec_modQ(key.h);
    std::vector<int16_t> sig(n);
// use a random hash of message
    std::vector<uint16_t> hm(n);
    for (uint64_t i=0; i<n; ++i) {
        hm[i]=rand()%F_Q;
    }
    std::vector<int8_t> sample1(n);
    std::vector<int8_t> sample2(n);
    std::vector<uint16_t> sample_target(n);
    std::vector<fpr> f_FFT(n);
    std::vector<fpr> g_FFT(n);
    std::vector<fpr> F_FFT(n);
    std::vector<fpr> G_FFT(n);

    for (auto _ : state) {
        sign_dyn_lazy_offline(&rng, key.f.data(), key.g.data(), key.F.data(), key.G.data(), key.h.data(), logn,
                              sample1.data(), sample2.data(), sample_target.data(), f_FFT.data(), g_FFT.data(),
                              F_FFT.data(), G_FFT.data());
    }
}

static void lazy_online(benchmark::State& state) {
    // Perform setup here
    const uint64_t logn = 9;
    const uint64_t n = 1 << logn;
    inner_shake256_context rng;
    inner_shake256_init(&rng);
    falcon_key_t key = keygen(logn, &rng);
    vec_modQ hq = to_vec_modQ(key.h);
    std::vector<int16_t> sig(n);
    // use a random hash of message
    std::vector<uint16_t> hm(n);
    for (uint64_t i=0; i<n; ++i) {
        hm[i]=rand()%F_Q;
    }
    std::vector<int8_t> sample1(n);
    std::vector<int8_t> sample2(n);
    std::vector<uint16_t> sample_target(n);
    std::vector<fpr> f_FFT(n);
    std::vector<fpr> g_FFT(n);
    std::vector<fpr> F_FFT(n);
    std::vector<fpr> G_FFT(n);

    for (uint64_t i=0; i<1; ++i) {
        sign_dyn_lazy_offline(&rng, key.f.data(), key.g.data(), key.F.data(), key.G.data(), key.h.data(), logn,
                              sample1.data(), sample2.data(), sample_target.data(), f_FFT.data(), g_FFT.data(),
                              F_FFT.data(), G_FFT.data());
    }
    std::vector<uint16_t> orig_sample_target = sample_target;
    for (auto _ : state) {
        sign_dyn_lazy_online(sample1.data(), sample2.data(), sample_target.data(), sig.data(),
                             f_FFT.data(), g_FFT.data(), F_FFT.data(), G_FFT.data(), hm.data(), logn,
                             nullptr);
    }
}

static void falcon_orig(benchmark::State& state) {
    const uint64_t logn = 9;
    const uint64_t n = 1 << logn;
    inner_shake256_context rng;
    inner_shake256_init(&rng);
    falcon_key_t key = keygen(logn, &rng);
    vec_modQ hq = to_vec_modQ(key.h);
    std::vector<int16_t> sig(n);
    std::vector<uint16_t> hm(n);
// use a random hash of message
    for (uint64_t i = 0; i < n; ++i) {
        hm[i] = rand() % F_Q;
    }
    uint8_t *tmp = (uint8_t *) aligned_alloc(64, 1024 * 1024);
    for (auto _ : state) {
        falcon_inner_sign_dyn(
                sig.data(),
                &rng,
                key.f.data(), key.g.data(),
                key.F.data(), key.G.data(),
                hm.data(),
                logn, tmp);
    }
    free(tmp);
}


// Register the function as a benchmark
BENCHMARK(lazy_offline);
BENCHMARK(lazy_online);
BENCHMARK(falcon_orig);
