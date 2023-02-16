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
    vector<int>seed_gen;
    int b_ = 0;
    int m_ = 0;
    int N_ = 0;
    int LEN = 13;
    uint32_t SEED = 0;


    vHLL(int N, int b) {
        b_ = b;
        m_ = (1 << b);
        N_ = N;
        ms_.resize(N);
        seed_gen.resize(0);
        srand((int)time(0));
        unordered_set<int> seed_set;
        while (seed_set.size() != m_) {
            int cur = rand();
            if (seed_set.count(cur) == 0) {
                seed_gen.push_back(cur);
                seed_set.insert(cur);
            }
        }
        /*for (int i = 0; i < m_; i++) cout << seed_gen[i] << " ";
        cout << endl;*/
    }

    void insert(const uint32_t* src, const uint32_t* dst) {
        uint32_t f = 0;
        uint32_t e = 0;
        MurmurHash3_x86_32(src, LEN, SEED, &f);
        MurmurHash3_x86_32(dst, LEN, SEED, &e);
        int p = e >> (32 - b_);
        int q = e - (p << (32 - b_));
        int leftmost_index = 0;
        while (q > 0) {
            leftmost_index++;
            q >>= 1;
        }
        leftmost_index = 32 - b_ - leftmost_index + 1;
        uint32_t index = 0;
        string tp = to_string(f ^ seed_gen[p]);
        MurmurHash3_x86_32((const uint32_t*)tp.c_str(), LEN, SEED, &index);
        index %= N_;
        ms_[index] = max(ms_[index], leftmost_index);
    }

    vector<double> getParm(const uint32_t* src, double n) const {
        uint32_t f = 0;
        MurmurHash3_x86_32(src, LEN, SEED, &f);
        double sum = 0;
        double zero = 0;
        vector<double> ret;
        for (int i = 0; i < m_; i++) {
            uint32_t index = 0;
            string tp = to_string(f ^ seed_gen[i]);
            MurmurHash3_x86_32((const uint32_t*)tp.c_str(), LEN, SEED, &index);
            index %= N_;
            ret.push_back(ms_[index]);
            sum += 1.0 / (1 << ms_[index]);
            if (ms_[index] == 0) zero++;
        }
        double ave = m_ / sum;
        double ns = 1.0 * alpha(m_) * m_ * ave;
        if (ns < 2.5 * m_ && zero != 0) {
            ns = -m_ * log(zero/m_);
        }
        double nf = (1.0 * N_ * m_ / (N_ - m_)) * (ns / m_ - n / N_);
        double nn = round(nf + 0.5);
        if (nn < 1) nn = 1;
        ret.push_back(nn);
        return ret;
    }

    double getTotal() {
        double sum = 0, zero = 0;
        for (int i = 0; i < N_; i++) {
            sum += 1.0 / (1 << ms_[i]);
            if (ms_[i] == 0) zero++;
        }
        double ave = (double)N_ / sum;
        double ns = 1.0 * alpha(N_) * (double)N_ * ave;
        if (ns < 2.5 * N_ && zero != 0) {
            ns = -N_ * log(zero/ N_);
        }
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
};

void pkt_sampling_driver(int N, int b, int P, string res_data) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train//";
    string train[] = {"00.txt","01.txt","02.txt","03.txt","04.txt"};
    string test[] = { "05.txt","06.txt"};
    string buf = "1",t1 = "1", t2 = "1";
    char idx = '1';
    for (auto& p : test) {
        vHLL hll = vHLL(N, b);
        inf.open(pkt_data+p, ios::in);
        ouf.open(res_data + "//"+to_string(P) + "ret" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0;
        while (getline(inf, buf)) {
            if (sample == 10) sample = 0;
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t2 = buf.substr(0, i);
                    t1 = buf.substr(i + 1);
                    break;
                }
            }
            flow[t1].insert(t2);
            if (sample++ < P/10) hll.insert((const uint32_t*)t1.c_str(), (const uint32_t*)t2.c_str());
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
            double nf = esi[32], ave = esi[33];
            am = max(am, ss);
            em = max(em, nf);
            ouf << ss << " " << nf << " " << ave << " ";
            for (int i = 0; i < hll.m_; i++) ouf << esi[i] << " ";
            ouf << endl;
        }
        cout << am << " " << em << endl;
        inf.close();
        ouf.close();
        break;
    }
    cout << "Completed!"  << endl;
}

void element_sampling_driver(int N, int b, int P, string res_data) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train//";
    string txt[] = { "06.txt","01.txt","02.txt","03.txt","04.txt" };
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
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t2 = buf.substr(0, i);
                    t1 = buf.substr(i + 1);
                    break;
                }
            }
            flow[t1].insert(t2);
            count++;
            uint32_t e = 0;
            MurmurHash3_x86_32((const uint32_t*)t2.c_str(), 32, 1331, &e);
            if (100.0 * e <= 1.0 * UINT32_MAX * P) {
                sample++;
                hll.insert((const uint32_t*)t1.c_str(), (const uint32_t*)t2.c_str());
            }
        }
        cout << flow.size() << " " << count << " " << sample << endl;
        double n = hll.getTotal();
        double em = 0;
        cout << n << endl;
        size_t am = 0;
        for (auto& cur : flow) {
            string f = cur.first;
            uint32_t ss = cur.second.size();
            vector<double> esi = hll.getParm((const uint32_t*)f.c_str(), n);
            double nf = esi[(1<<b)] * 100 / P;
            am = max(am, ss);
            em = max(em, nf);
            ouf << ss << " " << nf <<  endl;
        }
        cout << am << " " << em << endl;
        inf.close();
        ouf.close();
        break;
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
                    t2 = buf.substr(0, i);
                    t1 = buf.substr(i + 1);
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
        if (cnt == 0) break;
        ++cnt;
    }
    /*cout << "pkt:" << flowNum << " "
        "P:" << 1.0 * P / 100 << " "
        "thoughout:" << (flowNum / 1000000.0) / runt << "mpps " 
        "ave_insert_time"<< runt*1000000000.0/flowNum << "ns" << endl*/
    cout << (flowNum / 1000000.0) / runt << ", ";
    if (P == 100) cout << endl;
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
                    t2 = buf.substr(0, i);
                    t1 = buf.substr(i + 1);
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
            hll.insert((const uint32_t*)pkt[i][0].c_str(), (const uint32_t*)pkt[i][1].c_str());
        }
        clock_t time2 = clock();
        runt += (double)(time2 - time1) / CLOCKS_PER_SEC;
        inf.close();
        if (cnt == 0) break;
        ++cnt;
    }
    /*cout << "pkt:" << flowNum << " "
        "P:" << 1.0 * P / 100 << " "
        "thoughout:" << (flowNum / 1000000.0) / runt << "mpps "
        "ave_insert_time"<< runt*1000000000.0/flowNum << "ns" << endl*/
    cout << (flowNum / 1000000.0) / runt << ", ";
    if (P == 100) cout << endl;
}

void vHLL_driver(int N, int b, string res_data) {
    ifstream inf;
    ofstream ouf;
    string pkt_data = "E://DeskTop//train//";
    string test[] = { "00.txt"};
    string buf = "1", t1 = "1", t2 = "1";
    char idx = '1';
    for (auto& p : test) {
        vHLL hll = vHLL(N, b);
        inf.open(pkt_data + p, ios::in);
        ouf.open(res_data + "//" + "size" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0;
        while (getline(inf, buf)) {
            if (sample == 10) sample = 0;
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t2 = buf.substr(0, i);
                    t1 = buf.substr(i + 1);
                    break;
                }
            }
            flow[t1].insert(t2);
            hll.insert((const uint32_t*)t1.c_str(), (const uint32_t*)t2.c_str());
        }
        double n = hll.getTotal();
        for (auto& cur : flow) {
            string f = cur.first;
            uint32_t ss = cur.second.size();
            vector<double> esi = hll.getParm((const uint32_t*)f.c_str(), n);
            double nf = esi[(1<<b)];
            ouf << ss << " " << nf << " " << 0 << " ";
            for (int i = 0; i < hll.m_; i++) ouf << esi[i] << " ";
            ouf << endl;
        }
        inf.close();
        ouf.close();
        break;
    }
    cout << "Completed!" << endl;
}

void count_pkt_nums() {
    string train[] = { "00.txt","01.txt","02.txt","03.txt","04.txt" };
    string var[] = { "05.txt" };
    string test[] = { "06.txt"};
    ifstream inf;
    string pkt_data = "E://DeskTop//train//";
    string buf, t1, t2;
    int a = 0, b = 0, c = 0;
    for (auto& p : train) {
        inf.open(pkt_data + p, ios::in);
        while (getline(inf, buf)) a++;
        inf.close();
    }
    for (auto& p : var) {
        inf.open(pkt_data + p, ios::in);
        while (getline(inf, buf)) b++;
        inf.close();
    }
    for (auto& p : test) {
        inf.open(pkt_data + p, ios::in);
        while (getline(inf, buf)) c++;
        inf.close();
    }
    cout << a << " " << b << " " << c << endl;
}

int main() {
    int p[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    string pkt = "e://desktop//res";
    string element = "e://desktop//res//element";
    string KB256 = "e://desktop//res//KB256";
    string test = "e://desktop//res//test";
    //vHLL_driver(1677722, 9, pkt);
    for (int i = 0; i < 10; i++) {
        //pkt_sampling_Throughout(1677722, 9, p[i]);
        //pkt_sampling_driver(1677722, 5, p[i], pkt);
        //pkt_sampling_driver(419430, 5, p[i], KB256);
        //pkt_sampling_driver(1677722, 5, p[i], test);
        //pkt_sampling_driver(419430, 5, p[i], test);
    }
    for (int i = 3; i < 7; i++) {
        //element_sampling_Throughout(1677722, 9, p[i]);
        element_sampling_driver(1677722, 9, p[i], element);
        //element_sampling_driver(419430, 5, p[i], element);
    }
}
