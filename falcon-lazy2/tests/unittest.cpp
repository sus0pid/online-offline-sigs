#include <gtest/gtest.h>
#include <vector>
#include <random>

#define restrict
extern "C" {
#include "../inner.h"
}

#define EXPORT extern "C"

#define REQUIRE_DRAMATICALLY(req_contition, error_msg)                                                        \
  do {                                                                                                        \
    if (!(req_contition)) {                                                                                   \
      std::cerr << "REQUIREMENT FAILED at " << __FILE__ << ":" << __LINE__ << ": " << error_msg << std::endl; \
      abort();                                                                                                \
    }                                                                                                         \
  } while (0)

#define F_Q 12289

int64_t posmod(int64_t a, int64_t q) {
    int64_t t = a%q;
    return t<0?t+q:t;
}
int64_t centermod(int64_t a, int64_t q) {
    int64_t t = a%q;
    if (t>q/2) return t-q;
    if (t<-q/2) return t+q;
    return t;
}

struct modQ {
    int64_t v;
    modQ(int64_t x=0): v(posmod(x,F_Q)) {}
};

modQ operator+(modQ a, modQ b) { return posmod(a.v+b.v, F_Q); }
modQ operator-(modQ a, modQ b) { return posmod(a.v-b.v, F_Q); }
modQ operator*(modQ a, modQ b) { return posmod(a.v*b.v, F_Q); }
modQ& operator+=(modQ& a, modQ b) { a.v = posmod(a.v+b.v, F_Q); return a; }
modQ& operator-=(modQ& a, modQ b) { a.v = posmod(a.v-b.v, F_Q); return a; }
modQ& operator*=(modQ& a, modQ b) { a.v = posmod(a.v*b.v, F_Q); return a; }
bool operator==(modQ& a, modQ b) { return a.v == b.v; }
std::ostream& operator<<(std::ostream& out, modQ x) { return out << x.v; }

void convert(modQ* res, const int8_t* x, uint64_t n) {
    for (uint64_t k=0; k<n; ++k) {
        res[k]=modQ(x[k]);
    }
}
void convert(modQ* res, const int16_t* x, uint64_t n) {
    for (uint64_t k=0; k<n; ++k) {
        res[k]=modQ(x[k]);
    }
}
void convert(modQ* res, const uint16_t* x, uint64_t n) {
    for (uint64_t k=0; k<n; ++k) {
        res[k]=modQ(x[k]);
    }
}
void starproduct(modQ* res, const modQ* x, const modQ* y, uint64_t n) {
    for (uint64_t k=0; k<n; ++k) {
        modQ rk = 0;
        for (uint64_t i=0; i<=k; ++i) {
            rk += x[k-i] * y[i];
        }
        for (uint64_t i=k+1; i<n; ++i) {
            rk -= x[n+k-i] * y[i];
        }
        res[k]=rk;
    }
}

typedef std::vector<modQ> vec_modQ;
vec_modQ to_vec_modQ(const std::vector<int8_t>& v) {
    vec_modQ res(v.size());
    convert(res.data(), v.data(), v.size());
    return res;
}
vec_modQ to_vec_modQ(const std::vector<uint16_t>& v) {
    vec_modQ res(v.size());
    convert(res.data(), v.data(), v.size());
    return res;
}
vec_modQ to_vec_modQ(const std::vector<int16_t>& v) {
    vec_modQ res(v.size());
    convert(res.data(), v.data(), v.size());
    return res;
}
vec_modQ starproduct(const vec_modQ& x, const vec_modQ& y) {
    uint64_t n = x.size();
    REQUIRE_DRAMATICALLY(y.size()==n, "wrong size");
    vec_modQ res(n);
    starproduct(res.data(), x.data(), y.data(), n);
    return res;
}
vec_modQ operator-(const vec_modQ& x, const vec_modQ& y) {
    uint64_t n = x.size();
    REQUIRE_DRAMATICALLY(y.size()==n, "wrong size");
    vec_modQ res(n);
    for (uint16_t i=0; i<n; ++i) {
        res[i] = x[i] - y[i];
    }
    return res;
}
vec_modQ operator+(const vec_modQ& x, const vec_modQ& y) {
    uint64_t n = x.size();
    REQUIRE_DRAMATICALLY(y.size()==n, "wrong size");
    vec_modQ res(n);
    for (uint16_t i=0; i<n; ++i) {
        res[i] = x[i] + y[i];
    }
    return res;
}
bool operator==(const vec_modQ& x, const vec_modQ& y) {
    uint64_t n = x.size();
    REQUIRE_DRAMATICALLY(y.size()==n, "wrong size");
    for (uint64_t i=0; i<n; ++i) {
        if (x[i].v != y[i].v) return false;
    }
    return true;
}
std::default_random_engine& randgen() {
    static std::default_random_engine eng;
    return eng;
}
uint64_t random_u64() {
    static std::uniform_int_distribution<uint64_t> dist;
    return dist(randgen());
}
uint64_t random_i64() {
    static std::uniform_int_distribution<int64_t> dist;
    return dist(randgen());
}

double print_statistics(const std::vector<double>& vec) {
    uint64_t n = vec.size();
    double minv = 1./0.;
    double maxv = -1./0.;
    double sumv = 0;
    double sumsq = 0;
    for (uint64_t i=0; i<n; ++i) {
        double v = vec[i];
        if (v > maxv) maxv = v;
        if (v < minv) minv = v;
        sumv += v;
        sumsq += v * v;
    }
    double mean = sumv/n;
    double var = sumsq/n - mean*mean;
    std::cout
            << "min: " << minv << "\n"
            << "max: " << maxv << "\n"
            << "mean: " << mean << "\n"
            << "stdv: " << sqrt(var) << "\n"
            << "norm: " << sqrt(sumsq) << std::endl;
    return sqrt(sumsq);
}

struct falcon_key_t {
    std::vector<int8_t> f;
    std::vector<int8_t> g;
    std::vector<int8_t> F;
    std::vector<int8_t> G;
    std::vector<uint16_t> h;
};

falcon_key_t keygen(uint64_t logn, inner_shake256_context* rng) {
    falcon_key_t res;
    const uint64_t n = 1<<logn;
    uint8_t* tmp = (uint8_t*) aligned_alloc(64, 1024*1024);
    res.f.resize(n);
    res.g.resize(n);
    res.F.resize(n);
    res.G.resize(n);
    res.h.resize(n);
    Zf(keygen)(
            rng,
            res.f.data(),res.g.data(),
            res.F.data(),res.G.data(),
            res.h.data(),logn, tmp);
    free(tmp);
    return res;
}

TEST(falcon, keygen) {
    const uint64_t logn = 9;
    //const uint64_t n = 1 << logn;
    inner_shake256_context rng;
    inner_shake256_init(&rng);
    falcon_key_t key = keygen(logn, &rng);
    // verify that fh=g mod q
    // verify that Fh=G mod q
    vec_modQ fq = to_vec_modQ(key.f);
    vec_modQ gq = to_vec_modQ(key.g);
    vec_modQ Fq = to_vec_modQ(key.F);
    vec_modQ Gq = to_vec_modQ(key.G);
    vec_modQ hq = to_vec_modQ(key.h);
    ASSERT_EQ(starproduct(fq, hq), gq);
    ASSERT_EQ(starproduct(Fq, hq), Gq);
}

TEST(falcon, original_sig) {
    const uint64_t logn = 9;
    const uint64_t n = 1 << logn;
    inner_shake256_context rng;
    inner_shake256_init(&rng);
    falcon_key_t key = keygen(logn, &rng);
    vec_modQ hq = to_vec_modQ(key.h);
    std::vector<int16_t> sig(n);
    std::vector<uint16_t> hm(n);
    // use a random hash of message
    for (uint64_t i=0; i<n; ++i) {
        hm[i]=rand()%F_Q;
    }
    uint8_t* tmp = (uint8_t*) aligned_alloc(64, 1024*1024);
    falcon_inner_sign_dyn(
            sig.data(),
            &rng,
            key.f.data(), key.g.data(),
            key.F.data(), key.G.data(),
            hm.data(),
            logn, tmp);
    free(tmp);
    // compute the full uncompressed signature
    vec_modQ sigq = to_vec_modQ(sig);
    vec_modQ hsigq = starproduct(hq, sigq);
    std::vector<double> full_sig(2*n);
    for (uint64_t i=0; i<n; i++) {
        full_sig[i] = centermod(hm[i]-hsigq[i].v, F_Q);
        full_sig[i+n] = sig[i];
    }
    // print statistics about this vector and check the norm
    double norm = print_statistics(full_sig);
    ASSERT_LE(norm, 6000);
}

TEST(falcon, lazy_sig) {
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
    uint8_t* tmp = (uint8_t*) aligned_alloc(64, 1024*1024);
    falcon_inner_sign_dyn_lazy(
            sig.data(),
            &rng,
            key.f.data(), key.g.data(),
            key.F.data(), key.G.data(),
            key.h.data(),
            hm.data(),
            logn, tmp);
    free(tmp);
    // compute the full uncompressed signature
    vec_modQ sigq = to_vec_modQ(sig);
    vec_modQ hsigq = starproduct(hq, sigq);
    std::vector<double> full_sig(2*n);
    for (uint64_t i=0; i<n; i++) {
        full_sig[i] = centermod(hm[i]-hsigq[i].v, F_Q);
        full_sig[i+n] = sig[i];
    }
    // print statistics about this vector and check the norm
    double norm = print_statistics(full_sig);
    ASSERT_LE(norm, 6000);
}

EXPORT void sample_gaussian(int8_t *res,
                     sampler_context* spc,
                     fpr isigma,
                     unsigned logn);

TEST(falcon, sample_gaussian) {
    for (const uint64_t logn: {9,10}) {
        const uint64_t n = 1 << logn;
        inner_shake256_context rng;
        inner_shake256_init(&rng);
        double sigma = 2.; // between 1 and 2
        fpr isigma = FPR(1./sigma);
        sampler_context sc;
        Zf(prng_init)(&sc.p, &rng);
        sc.sigma_min = fpr_sigma_min[logn];
        std::vector<int8_t> res(n);
        sample_gaussian(res.data(), &sc, isigma, logn);
        std::vector<double> st(n);
        for (uint64_t i=0; i<n; ++i) st[i]=res[i];
        double norm = print_statistics(st);
        ASSERT_LE(norm, 3*n);
    }
}

/** x0 - h.x1 */
EXPORT void compute_target(
        uint16_t* res,
        uint16_t* h,
        int8_t* x0, int8_t* x1, unsigned logn);

TEST(falcon, compute_target) {
    for (const uint64_t logn: {9,10}) {
        const uint64_t n = 1 << logn;
        std::vector<uint16_t> actual(n);
        std::vector<uint16_t> h(n);
        std::vector<int8_t> x0(n);
        std::vector<int8_t> x1(n);
        for (uint16_t i=0; i<n; ++i) {
            h[i] = posmod(random_u64(), F_Q);
            x0[i] = centermod(random_u64(), 65);
            x1[i] = centermod(random_u64(), 65);
        }
        vec_modQ hq = to_vec_modQ(h);
        vec_modQ x0q = to_vec_modQ(x0);
        vec_modQ x1q = to_vec_modQ(x1);
        vec_modQ expect = x0q - starproduct(x1q , hq);
        compute_target(actual.data(), h.data(), x0.data(), x1.data(), logn);
        ASSERT_EQ(to_vec_modQ(actual), expect);
    }
}
