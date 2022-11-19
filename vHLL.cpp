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
    uint32_t SEED = 131;

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

    vector<double> estimate(const char* src, double n) const {
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
        double ave = m_ / sum;
        //ave = 0.008661963939488782 * ave * ave + 1.2152514269450094 * ave - 0.21991023804331547;
        double ns = 1.0 * alpha(m_) * m_ * ave;
        double nf = (1.0 * N_ * m_ / (N_ - m_)) * (ns / m_ - n / N_);
        vector<double> ret = {m_ / sum, round(nf+0.5)};
        return ret;
    }

    double getTotal() {
        double sum = 0, zero = 0;
        for (int i = 0; i < N_; i++) {
            sum += 1.0 / (1 << ms_[i]);
        }
        double ave = N_ / sum;
        //ave = 0.008661963939488782 * ave * ave + 1.2152514269450094 * ave -0.21991023804331547;
        double ns = 1.0 * alpha(N_) * N_ * ave;
        return ns;
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
            j <<= 1;      
            j |= (n & 1);  
            n >>= 1;     
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

void driver(int N, int b, int P, string res_data) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train";
    const char* src;
    const char* dst;
    string buf, t1, t2;
    char idx = '1';
    double runt = 0;
    for (auto& p : filesystem::directory_iterator(pkt_data)) {
        vHLL hll = vHLL(N, b);
        inf.open(p.path().string(), ios::in);
        //ouf.open(res_data + "//size" + idx + ".txt"); idx++;
        ouf.open(res_data + "//"+to_string(P) + "ret" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0;
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
        cout << flow.size() << endl;
        double n = hll.getTotal();
        double em = 0;
        cout << n << endl;
        size_t am = 0;
        for (auto& cur : flow) {
            string f = cur.first;
            unordered_set<string> ss = cur.second;
            vector<double> esi = hll.estimate(f.c_str(), n);
            double ave = esi[0], nf = esi[1];
            am = max(am, ss.size());
            em = max(em, nf);
            //ouf << ss.size() << " " << nf << endl;
            ouf << ss.size() << " " << nf << " " << ave << " " << (double)P/100 << endl;
        }
        cout << am << " " << em << endl;
        inf.close();
        ouf.close();
    }
    cout << "Completed!"  << endl;
}

void Throughout(int N, int b, int P) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train";
    const char* src;
    const char* dst;
    string buf, t1, t2;
    double runt = 0;
    int sample = 0, flowNum = 0;
    for (auto& p : filesystem::directory_iterator(pkt_data)) {
        vHLL hll = vHLL(N, b);
        inf.open(p.path().string(), ios::in);
        while (getline(inf, buf)) {
            flowNum++;
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
            auto start = std::chrono::system_clock::now();
            if (sample++ < P) hll.insert(src, dst);
            auto end = std::chrono::system_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            runt += (double)(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
        }
        inf.close();
    }
    cout << "pkt_num:" << flowNum << " "
        "sampling_rate:" << 1.0 * P / 100 << " "
        "vHLL_insert_time:" << runt << endl;
}

int main() {
    /*vHLL hll = vHLL(1677722, 5);
    const char *ip = "10.23.44.1";
    const char *d = "192.111.12.10";
    hll.insert(ip, d);
    double n = hll.getTotal();
    cout << n << endl;
    cout << hll.estimate(ip, n)[1] << endl;
    auto start = std::chrono::system_clock::now();
    if (n == 0) n++;
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    cout << (double)(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;*/
    int p[5] = {10, 30, 50, 80, 100};
    string r1 = "E://DeskTop//res";
    string r2 = "E://DeskTop//res//base";
    for (int i = 3; i < 4; i++) {
        //Throughout(1677722, 5, P);
        driver(1677722, 5, p[i], r1);
    }
    /*for (int i = 0; i < 100; i++) {
        uint32_t res;
        MurmurHash3_x86_32(ip, 32, 313, &res);
        cout << res << endl;
    }*/
}