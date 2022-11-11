#include <cmath>
#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream> 
#include <iostream>
#include "murmur3.h"
#include <sstream>
#include <filesystem>
#include <windows.h>
using namespace std;

class vHLL {
public:
    vector<int> ms_;
    int b_;
    int m_;
    int N_;
    uint32_t SEED = 313;

    vHLL(int N, int b) {
        b_ = b;
        m_ = (1 << b);
        N_ = N;
        ms_.resize(N);
    }

    void insert(const char* src, const char* dst) {
        uint32_t f = 0;
        uint32_t e = 0;
        MurmurHash3_x86_32(src, 32, SEED, &f);
        MurmurHash3_x86_32(dst, 32, SEED, &e);
        int p = e & ((1 << b_) - 1);
        int q = nlz5(e >> b_) + 1;
        uint32_t index = 0;
        string tp = to_string(f ^ p);
        const char* key = tp.c_str();
        MurmurHash3_x86_32(key, 32, SEED, &index);
        index %= N_;
        ms_[index] = max(ms_[index], q);
    }

    double estimate(const char* src, double n) const {
        uint32_t f = 0;
        MurmurHash3_x86_32(src, 32, SEED, &f);
        double sum = 0;
        double zero = 0;
        for (int i = 0; i < m_; i++) {
            uint32_t index = 0;
            string tp = to_string(f ^ i);
            const char* key = tp.c_str();
            MurmurHash3_x86_32(key, 32, SEED, &index);
            index %= N_;
            sum += 1.0 / (1 << ms_[index]);
            if (ms_[index] == 0) zero++;
        }
        double ns = 1.0 * alpha(m_) * m_ * m_ / sum;
        if (ns < 2.5 * m_ && zero != 0) {
            ns = -m_*log2(zero / m_);
        }
        double nf = (1.0 * N_ * m_ / (N_ - m_)) * (ns / m_ - n / N_);
        return round(nf + 0.5);
    }

    double getTotal() {
        double sum = 0;
        for (int i = 0; i < N_; i++) sum += 1.0/(1 << ms_[i]);
        return 1.0 * alpha(N_) * N_ * N_ / sum;
    }

    double alpha(int t) const {
        switch (t) {
        case 16:
            return 0.673;
        case 32:
            return 0.697;
        case 64:
            return 0.709;
        default:
            return 0.7213 / (1 + 1.079 / t);
        }
    }

    uint32_t reverseBits(uint32_t n)
    {
        uint32_t j = 0;
        for (int i = 0; i < 32; i++)
        {
            j <<= 1;      //j���ƣ��ұ߲�0
            j |= (n & 1);  //��n�����λȡ�����
            n >>= 1;      //n����
        }
        return j;
    }

    // These from the book "Hacker's Delight"
    uint32_t pop(uint32_t x) {
        x = x - ((x >> 1) & 0x55555555);
        x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
        x = (x + (x >> 4)) & 0x0F0F0F0F;
        x = x + (x << 8);
        x = x + (x << 16);
        return x >> 24;
    }

    uint32_t nlz5(uint32_t x) {
        x = reverseBits(x);
        x = x | (x >> 1);
        x = x | (x >> 2);
        x = x | (x >> 4);
        x = x | (x >> 8);
        x = x | (x >> 16);
        return pop(~x);
    }
};

void driver(int N, int b, int P) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train";
    string res_data = "E://DeskTop//res";
    const char* src;
    const char* dst;
    string buf, t1, t2;
    char idx = '1';
    double runt = 0;
    for (auto& p : filesystem::directory_iterator(pkt_data)) {
        vHLL hll = vHLL(N, b);
        inf.open(p.path().string(), ios::in);
        ouf.open(res_data + "//size" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0;
        LARGE_INTEGER tt1, tt2, tc;
        QueryPerformanceFrequency(&tc);
        QueryPerformanceCounter(&tt1);
        while (getline(inf, buf)) {
            if (sample == 100) sample = 0;
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t1 = buf.substr(0, i);
                    t2 = buf.substr(i + 1);
                    break;
                }
            }
            src = t1.c_str();
            dst = t2.c_str();
            flow[src].insert(dst);
            if (sample++ < P) hll.insert(src, dst);
        }
        QueryPerformanceCounter(&tt2);
        runt += (double)(tt2.QuadPart - tt1.QuadPart) / (double)tc.QuadPart;
        cout << flow.size() << endl;
        double n = hll.getTotal();
        double em = 0;
        size_t am = 0;
        for (auto& cur : flow) {
            string f = cur.first;
            unordered_set<string> ss = cur.second;
            double esi = hll.estimate(f.c_str(), n);
            am = max(am, ss.size());
            em = max(em, esi);
            ouf << ss.size() << " " << esi << endl;
        }
        cout << am << " " << em << endl;
        inf.close();
        ouf.close();
    }
    cout << "����ʱ�䣺" << runt / 5 << endl;
}

int main() {
    vHLL hll = vHLL(1677722, 5);
    const char *ip = "10.23.44.1";
    const char *d = "192.111.12.10";
    hll.insert(ip, d);
    double n = hll.getTotal();
    cout << n << endl;
    cout << hll.estimate(ip, n) << endl;
    //driver(1677722, 5, 10);
    /*for (int i = 0; i < 100; i++) {
        uint32_t res;
        MurmurHash3_x86_32(ip, 32, 313, &res);
        cout << res << endl;
    }*/
}