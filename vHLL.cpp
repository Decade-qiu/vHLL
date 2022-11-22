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
    int b_ = 0;
    int m_ = 0;
    int N_ = 0;
    int LEN = 32;
    uint32_t SEED = 0;

    vHLL(int N, int b) {
        b_ = b;
        m_ = (1 << b);
        N_ = N;
        ms_.resize(N);
    }

    void insert(const uint32_t* src, const uint32_t* dst) {
        uint32_t f = 0;
        uint32_t e = 0;
        MurmurHash3_x86_32(src, LEN, SEED, &f);
        MurmurHash3_x86_32(dst, LEN, SEED, &e);
        int p = e & ((1 << b_) - 1);
        int q = nlz5(e >> b_) + 1;
        uint32_t index = 0;
        string tp = to_string(f ^ p);
        MurmurHash3_x86_32((const uint32_t*)tp.c_str(), LEN, SEED, &index);
        index %= N_;
        ms_[index] = max(ms_[index], q);
    }

    double estimate(const uint32_t* src, double n) const {
        uint32_t f = 0;
        MurmurHash3_x86_32(src, LEN, SEED, &f);
        double sum = 0;
        double zero = 0;
        for (int i = 0; i < m_; i++) {
            uint32_t index = 0;
            string tp = to_string(f ^ i);
            MurmurHash3_x86_32((const uint32_t*)tp.c_str(), LEN, SEED, &index);
            index %= N_;
            sum += 1.0 / (1 << ms_[index]);
            if (ms_[index] == 0) zero++;
        }
        double ave = m_ / sum;
        //ave = 0.008661963939488782 * ave * ave + 1.2152514269450094 * ave - 0.21991023804331547;
        double ns = 1.0 * alpha(m_) * m_ * ave;
        double nf = (1.0 * N_ * m_ / (N_ - m_)) * (ns / m_ - n / N_);
        return nf;
    }

    vector<double> getParm(const uint32_t* src, double n) const {
        uint32_t f = 0;
        MurmurHash3_x86_32(src, LEN, SEED, &f);
        double sum = 0;
        double zero = 0;
        for (int i = 0; i < m_; i++) {
            uint32_t index = 0;
            string tp = to_string(f ^ i);
            MurmurHash3_x86_32((const uint32_t*)tp.c_str(), LEN, SEED, &index);
            index %= N_;
            sum += 1.0 / (1 << ms_[index]);
            if (ms_[index] == 0) zero++;
        }
        double ave = m_ / sum;
        //ave = 0.008661963939488782 * ave * ave + 1.2152514269450094 * ave - 0.21991023804331547;
        double ns = 1.0 * alpha(m_) * m_ * ave;
        double nf = (1.0 * N_ * m_ / (N_ - m_)) * (ns / m_ - n / N_);
        vector<double> ret = { 1.0*m_ / sum, round(nf + 0.5) };
        return ret;
    }

    double getTotal() {
        double sum = 0, zero = 0;
        for (int i = 0; i < N_; i++) {
            sum += 1.0 / pow(2, ms_[i]);
        }
        double ave = (double)N_ / sum;
        //ave = -0.0001697495894330087 * ave * ave + 0.03709623856264059 * ave * ave + 1.1906057074760585 * ave - 0.16254656144756013;
        //ave = -2.3771430362972123e-05*ave*ave + 0.008661963939488782 * ave * ave + 1.2152514269450094 * ave -0.21991023804331547;
        //ave = -1.5077544524883974e-06 * ave * ave + 0.0007378796456244439 * ave * ave + 1.128860364312081 * ave -0.1535521480349578;
        double ns = 1.0 * alpha(N_) * (double)N_ * ave;
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

void pkt_sampling_driver(int N, int b, int P, string res_data) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train//";
    string txt[] = {"00.txt","01.txt","02.txt","03.txt","04.txt"};
    string buf = "1",t1 = "1", t2 = "1";
    char idx = '1';
    double runt = 0;
    for (auto& p : txt) {
        vHLL hll = vHLL(N, b);
        inf.open(pkt_data+p, ios::in);
        //ouf.open(res_data + "//size" + idx + ".txt"); idx++;
        ouf.open(res_data + "//"+to_string(P) + "ret" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0, count = 0;
        while (getline(inf, buf)) {
            if (sample == 100) sample = 0;
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t1 = buf.substr(0, i);
                    t2 = buf.substr(i + 1);
                    break;
                }
            }
            flow[t1].insert(t2);
            if (sample++ < P) hll.insert((const uint32_t*)t1.c_str(), (const uint32_t*)t2.c_str());
        }
        cout << flow.size() << endl;
        double n = hll.getTotal();
        double em = 0;
        cout << n << endl;
        size_t am = 0;
        for (auto& cur : flow) {
            string f = cur.first;
            uint32_t ss = cur.second.size();
            vector<double> esi = hll.getParm((const uint32_t*)f.c_str(), n);
            double ave = esi[0], nf = esi[1];
            am = max(am, ss);
            em = max(em, nf);
            //ouf << ss.size() << " " << nf << endl;
            ouf << ss << " " << nf << " " << ave << " " << (double)P/100 << endl;
        }
        cout << am << " " << em << endl;
        inf.close();
        ouf.close();
    }
    cout << "Completed!"  << endl;
}

void element_sampling_driver(int N, int b, int P, string res_data) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train//";
    string txt[] = { "00.txt","01.txt","02.txt","03.txt","04.txt" };
    string buf = "1", t1 = "1", t2 = "1";
    char idx = '1';
    double runt = 0;
    for (auto& p : txt) {
        vHLL hll = vHLL(N, b);
        inf.open(pkt_data + p, ios::in);
        //ouf.open(res_data + "//size" + idx + ".txt"); idx++;
        ouf.open(res_data + "//" + to_string(P) + "ret" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0, count = 0;
        while (getline(inf, buf)) {
            if (sample == 100) sample = 0;
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t1 = buf.substr(0, i);
                    t2 = buf.substr(i + 1);
                    break;
                }
            }
            flow[t1].insert(t2);
            uint32_t e = 0;
            MurmurHash3_x86_32((const uint32_t*)t1.c_str(), 32, 0, &e);
            if (100.0*e <= 1.0*UINT32_MAX*P) hll.insert((const uint32_t*)t1.c_str(), (const uint32_t*)t2.c_str());
        }
        cout << flow.size() << endl;
        double n = hll.getTotal();
        double em = 0;
        cout << n << endl;
        size_t am = 0;
        for (auto& cur : flow) {
            string f = cur.first;
            uint32_t ss = cur.second.size();
            vector<double> esi = hll.getParm((const uint32_t*)f.c_str(), n);
            double ave = esi[0], nf = esi[1];
            am = max(am, ss);
            em = max(em, nf);
            //ouf << ss.size() << " " << nf << endl;
            ouf << ss << " " << nf << " " << ave << " " << (double)P / 100 << endl;
        }
        cout << am << " " << em << endl;
        inf.close();
        ouf.close();
    }
    cout << "Completed!" << endl;
}

void pkt_sampling_Throughout(int N, int b, int P) {
    ifstream inf;
    string pkt_data = "E://DeskTop//train//";
    string txt[] = { "00.txt","01.txt","02.txt","03.txt","04.txt" };
    string buf, t1, t2;
    double runt = 0;
    int sample = 0, flowNum = 0, cnt = 0;
    for (auto& p : txt) {
        vector<vector<string>> pkt;
        inf.open(pkt_data+p, ios::in);
        while (getline(inf, buf)) {
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t1 = buf.substr(0, i);
                    t2 = buf.substr(i + 1);
                    break;
                }
            }
            vector<string> tp = {t1, t2};
            pkt.push_back(tp);
        }
        flowNum += pkt.size();
        vHLL hll = vHLL(N, b);
        clock_t time1 = clock();
        for (int i = 0; i < pkt.size();i++, sample++) {
            if (sample == 100) sample = 0;
            if (sample < P) {
                hll.insert((const uint32_t*)pkt[i][0].c_str(), (const uint32_t*)pkt[i][1].c_str());
            }
        }
        clock_t time2 = clock();
        runt += (double)(time2 - time1) / CLOCKS_PER_SEC;
        inf.close();
        if (cnt == 1) break;
        ++cnt;
    }
    /*cout << "pkt:" << flowNum << " "
        "P:" << 1.0 * P / 100 << " "
        "thoughout:" << (flowNum / 1000000.0) / runt << "mpps " 
        "ave_insert_time"<< runt*1000000000.0/flowNum << "ns" << endl*/
    cout << (flowNum / 1000000.0) / runt << ", ";
    //cout << flowNum << " " << runt << endl;
}

void element_sampling_Throughout(int N, int b, int P) {
    ifstream inf;
    string pkt_data = "E://DeskTop//train//";
    string txt[] = { "00.txt","01.txt","02.txt","03.txt","04.txt" };
    string buf, t1, t2;
    double runt = 0;
    int sample = 0, flowNum = 0, cnt = 0, ss = 0, ma = 0;
    for (auto& p : txt) {
        vector<vector<string>> pkt;
        inf.open(pkt_data + p, ios::in);
        while (getline(inf, buf)) {
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t1 = buf.substr(0, i);
                    t2 = buf.substr(i + 1);
                    break;
                }
            }
            vector<string> tp = { t1, t2};
            pkt.push_back(tp);
        }
        flowNum += pkt.size();
        vHLL hll = vHLL(N, b);
        clock_t time1 = clock();
        for (int i = 0; i < pkt.size(); i++) {
            uint32_t e = 0;
            MurmurHash3_x86_32((const uint32_t*)pkt[i][1].c_str(), 32, 0, &e);
            if (e < pow(2, 32)*P/100) {
                hll.insert((const uint32_t*)pkt[i][0].c_str(), (const uint32_t*)pkt[i][1].c_str());
            }
        }
        clock_t time2 = clock();
        runt += (double)(time2 - time1) / CLOCKS_PER_SEC;
        inf.close();
        if (cnt == 1) break;
        ++cnt;
    }
    /*cout << "pkt:" << flowNum << " "
        "P:" << 1.0 * P / 100 << " "
        "thoughout:" << (flowNum / 1000000.0) / runt << "mpps "
        "ave_insert_time"<< runt*1000000000.0/flowNum << "ns" << endl*/
    cout << (flowNum / 1000000.0) / runt << ", ";
}

int main() {
    int p[5] = {10, 30, 50, 80, 100};
    string pkt = "e://desktop//res";
    string element = "e://desktop//res//element";
    string r = "e://desktop//res//base";
    for (int i = 0; i < 4; i++) {
        //pkt_sampling_Throughout(1677722, 5, p[i]);
        //element_sampling_Throughout(1677722, 5, p[i]);
        //pkt_sampling_driver(1677722, 5, p[i], pkt);
        element_sampling_driver(1677722, 5, p[i], element);
    }
    /*vHLL hll = vHLL(1677722, 5);
    const char *ip = "10.23.44.1";
    const char *d = "192.111.12.10";
    hll.insert(ip, d);
    double n = hll.getTotal();
    cout << n << endl;
    cout << hll.estimate(ip, n)[1] << endl;*/
}
