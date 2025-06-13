#include <NTT.hpp>

using namespace mcl;
using std::vector;

static int log2(int N)
{
    int log2_N = 0; int temp = 1;
    while (temp < N)
    {
        temp *= 2; log2_N++;
    }
    return log2_N;
}

static Fr pow_P(const Fr& w, size_t a) {
    Fr p;
    Fr::pow(p, w, a);
    return p;
}

void NTT_initialize()
{
    Fr x;
    //order: 2^16
    x.setStr("3953746291237847556699976569425934098192757447134953284900320948457849631132", 10);
    g = x;
    inv_g = 1 / g;
}

vector<Fr> iter_NTT(const vector<Fr>& f)
{
    size_t N = f.size();
    if (N == 1) return f;
    vector<Fr> arr(N); int* arr1 = new int[N];
    size_t T = N / 2; size_t H = N; int temp = 1;
    int log2_N = log2(N);
    arr1[0] = 0;
    for (int i = 0; i < log2_N; i++)
    {
        for (size_t j = 0; j < N; j += H)
            arr1[j + T] = arr1[j] + temp;
        T /= 2; H /= 2; temp *= 2;
    }
    for (int i = 0; i < N; i++) arr[i] = f[arr1[i]];
    T = 2; H = 1;
    Fr* root = new Fr[log2_N];
    root[log2_N - 1] = pow_P(g, order / N);
    for (int i = log2_N - 2; i >= 0; i--)
        root[i] = (root[i + 1] * root[i + 1]);
    Fr* temp_arr = new Fr[N];
    for (int i = 0; i < log2_N; i++)
    {
        Fr TthRoot = root[i];
        for (size_t j = 0; j < N; j += T)
        {
            Fr t = 1;
            for (int k = 0; k < H; k++)
            {
                temp_arr[k] = arr[j + k] + arr[j + k + H] * t;
                temp_arr[k + H] = arr[j + k] - arr[j + k + H] * t;
                t *= TthRoot;
            }
            for (int k = 0; k < T; k++) arr[j + k] = temp_arr[k];
        }
        T *= 2; H *= 2;
    }
    delete[] arr1; delete[] temp_arr; delete[] root;
    return arr;
}

vector<Fr> inv_NTT(const vector<Fr>& f)
{
    size_t N = f.size();
    vector<Fr> arr(N); vector<Fr> res(N);
    for (int i = 0; i < N; i++) arr[i] = f[i];
    int log2_N = log2(N);
    size_t T = N; size_t H = N / 2;
    Fr w = pow_P(inv_g, order / N);
    for (int i = 0; i < log2_N; i++) {
        for (size_t j = 0; j < N; j += T)
        {
            Fr a = 1;
            for (int k = 0; k < H; k++)
            {
                res[j + k] = (arr[j + k] + arr[j + k + H]) / 2;
                res[j + k + H] = (arr[j + k] - arr[j + k + H]) / 2 * a;
                a *= w;
            }
        }
        w *= w; T /= 2; H /= 2;
        arr.swap(res);
    }
    size_t* temp = new size_t[N];
    temp[0] = 0; int t = 1;
    for (H = N / 2; H > 0; H /= 2) {
        for (int j = 0; j < t; j++) temp[j + t] = temp[j] + H;
        t *= 2;
    }
    for (int i = 0; i < N; i++) res[i] = arr[temp[i]];
    delete[] temp;
    return res;
}

vector<Fr> poly_mult(const vector<Fr>& f, const vector<Fr>& h)
{
    size_t df = f.size() - 1; size_t dh = h.size() - 1;
    size_t n = df + dh + 1;
    int N = 1;
    while (N < n) N *= 2;
    vector<Fr> fN(N); vector<Fr> hN(N);
    for (size_t i = 0; i < df + 1; i++) fN[i] = f[i];
    for (size_t i = df + 1; i < N; i++) fN[i] = 0;
    for (size_t i = 0; i < dh + 1; i++) hN[i] = h[i];
    for (size_t i = dh + 1; i < N; i++) hN[i] = 0;
    fN = iter_NTT(fN); hN = iter_NTT(hN);
    for (int i = 0; i < N; i++) fN[i] = fN[i] * hN[i];
    vector<Fr> tres = inv_NTT(fN);
    vector<Fr> res(tres.begin(), tres.begin() + n);
    return res;
}