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
    int LEN = 32;
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
        //string tp = to_string(f ^ p);
        MurmurHash3_x86_32((const uint32_t*)tp.c_str(), LEN, SEED, &index);
        index %= N_;
        ms_[index] = max(ms_[index], leftmost_index);
    }

    vector<double> getParm(const uint32_t* src, double n, int P) const {
        uint32_t f = 0;
        MurmurHash3_x86_32(src, LEN, SEED, &f);
        double sum = 0;
        double zero = 0;
        vector<double> ret;
        for (int i = 0; i < m_; i++) {
            uint32_t index = 0;
            string tp = to_string(f ^ seed_gen[i]);
            //string tp = to_string(f ^ i);
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
        if (P == 10) nf = +2.115340123799178e-17 * nf * nf * nf * nf * nf * nf * nf - 1.1256140728656618e-13 * nf * nf * nf * nf * nf * nf + 2.249795999138142e-10 * nf * nf * nf * nf * nf - 2.047892577515146e-07 * nf * nf * nf * nf + 7.935351082867359e-05 * nf * nf * nf - 0.008644351623126996 * nf * nf + 2.1973808873316045 * nf - 1.9577326577386387;
        else if (P == 20) nf = -1.1476009300998414e-18 * nf * nf * nf * nf * nf * nf * nf + 6.175660565393717e-15 * nf * nf * nf * nf * nf * nf - 1.1908710546093112e-11 * nf * nf * nf * nf * nf + 9.766204510632295e-09 * nf * nf * nf * nf - 3.6076510787704976e-06 * nf * nf * nf + 0.0011484450317673339 * nf * nf + 1.4221517162156627 * nf - 1.110003955048971;
        else if (P == 30) nf = -7.724100960823974e-20 * nf * nf * nf * nf * nf * nf * nf + 1.185213998384578e-15 * nf * nf * nf * nf * nf * nf - 4.991272942117753e-12 * nf * nf * nf * nf * nf + 8.549907906449684e-09 * nf * nf * nf * nf - 6.259473229940573e-06 * nf * nf * nf + 0.001894217467541572 * nf * nf + 1.1790722109507126 * nf - 1.0081630791990335;
        else if (P == 40) nf = +2.0517887861068614e-19 * nf * nf * nf * nf * nf * nf * nf - 1.4303112386843354e-15 * nf * nf * nf * nf * nf * nf + 3.325516034268529e-12 * nf * nf * nf * nf * nf - 2.7335177248263074e-09 * nf * nf * nf * nf + 1.6708135263850997e-07 * nf * nf * nf + 0.0005565316477764249 * nf * nf + 1.1225656488248998 * nf - 0.9589160649479493;
        else if (P == 50) nf = +9.057785861550707e-20 * nf * nf * nf * nf * nf * nf * nf - 7.560850022590245e-16 * nf * nf * nf * nf * nf * nf + 2.2295213708430094e-12 * nf * nf * nf * nf * nf - 2.732149572031114e-09 * nf * nf * nf * nf + 1.2325767096087433e-06 * nf * nf * nf - 6.009712185341652e-06 * nf * nf + 1.100770820309592 * nf - 1.116471531980017;
        else if (P == 60) nf = +3.720794992102691e-20 * nf * nf * nf * nf * nf * nf * nf - 3.1593786810697714e-16 * nf * nf * nf * nf * nf * nf + 9.255739100641945e-13 * nf * nf * nf * nf * nf - 1.0274828557120334e-09 * nf * nf * nf * nf + 2.2165877246052306e-07 * nf * nf * nf + 0.000244699661669841 * nf * nf + 1.026698180972445 * nf - 1.0239748196286018;
        else if (P == 70) nf = +9.43900546218361e-21 * nf * nf * nf * nf * nf * nf * nf - 8.941649490294502e-17 * nf * nf * nf * nf * nf * nf + 2.8174050149643037e-13 * nf * nf * nf * nf * nf - 2.909395298749796e-10 * nf * nf * nf * nf - 7.442696628945217e-08 * nf * nf * nf + 0.0002693612702674616 * nf * nf + 0.9846310308653047 * nf - 0.9941755252978641;
        else if (P == 80) nf = -2.730632379269233e-21 * nf * nf * nf * nf * nf * nf * nf + 3.875466779168853e-17 * nf * nf * nf * nf * nf * nf - 2.121694745242347e-13 * nf * nf * nf * nf * nf + 5.555881209755576e-10 * nf * nf * nf * nf - 7.038562424257159e-07 * nf * nf * nf + 0.00043317313202772677 * nf * nf + 0.9449526498354037 * nf - 0.9818122938855343;
        else if (P == 90) nf = +1.6779656926474064e-23 * nf * nf * nf * nf * nf * nf * nf + 2.8789688363077657e-19 * nf * nf * nf * nf * nf * nf - 9.978030410293229e-15 * nf * nf * nf * nf * nf + 6.701778064373672e-11 * nf * nf * nf * nf - 1.7817855535423915e-07 * nf * nf * nf + 0.00022876103285848727 * nf * nf + 0.9371138992890357 * nf - 1.115580832199844;
        else if (P == 100) nf = +4.807349035426421e-22 * nf * nf * nf * nf * nf * nf * nf - 7.94224238851793e-18 * nf * nf * nf * nf * nf * nf + 4.5269073699623554e-14 * nf * nf * nf * nf * nf - 1.0235531433111993e-10 * nf * nf * nf * nf + 4.796199857745417e-08 * nf * nf * nf + 0.00012496979534414095 * nf * nf + 0.9268351772522895 * nf - 1.1353611591862454;
        double nn = round(nf + 0.5);
        if (nn < 1) nn = 1;
        ret.push_back(nn);
        ret.push_back(ave);
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

    vector<double> Elem_getParm(const uint32_t* src, double n, int P) const {
        uint32_t f = 0;
        MurmurHash3_x86_32(src, LEN, SEED, &f);
        double sum = 0;
        double zero = 0;
        vector<double> ret;
        for (int i = 0; i < m_; i++) {
            uint32_t index = 0;
            string tp = to_string(f ^ seed_gen[i]);
            //string tp = to_string(f ^ i);
            MurmurHash3_x86_32((const uint32_t*)tp.c_str(), LEN, SEED, &index);
            index %= N_;
            ret.push_back(ms_[index]);
            sum += 1.0 / (1 << ms_[index]);
            if (ms_[index] == 0) zero++;
        }
        double ave = m_ / sum;
        double ns = 1.0 * alpha(m_) * m_ * ave;
        double nf = (1.0 * N_ * m_ / (N_ - m_)) * (ns / m_ - n / N_);
        double nn = round(nf + 0.5);
        ret.push_back(nn);
        ret.push_back(ave);
        return ret;
    }

    double Elem_getTotal(int P) {
        double sum = 0, zero = 0;
        for (int i = 0; i < N_; i++) {
            sum += 1.0 / (1 << ms_[i]);
            if (ms_[i] == 0) zero++;
        }
        double ave = (double)N_ / sum;
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
        //ouf.open(res_data + "//size" + idx + ".txt"); idx++;
        ouf.open(res_data + "//"+to_string(P) + "ret" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0, count = 0;
        while (getline(inf, buf)) {
            if (sample == 100) sample = 0;
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t2 = buf.substr(0, i);
                    t1 = buf.substr(i + 1);
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
            vector<double> esi = hll.getParm((const uint32_t*)f.c_str(), n, P);
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
    for (auto& p : txt) {
        vHLL hll = vHLL(N, b);
        inf.open(pkt_data + p, ios::in);
        //ouf.open(res_data + "//size" + idx + ".txt"); idx++;
        ouf.open(res_data + "//" + to_string(P) + "ret" + idx + ".txt"); idx++;
        unordered_map<string, unordered_set<string>> flow;
        int sample = 0, count = 0, sum = 0;
        while (getline(inf, buf)) {
            if (sample == 100) sample = 0;
            for (int i = 0; i < buf.size(); i++) {
                if (buf[i] == ' ') {
                    t2 = buf.substr(0, i);
                    t1 = buf.substr(i + 1);
                    break;
                }
            }
            flow[t1].insert(t2);
            uint32_t e = 0;
            MurmurHash3_x86_32((const uint32_t*)t2.c_str(), 32, 0, &e);
            sum++;
            if (e <= pow(2, 32) * P / 100) {
                count++;
                hll.insert((const uint32_t*)t1.c_str(), (const uint32_t*)t2.c_str());
            }
        }
        cout << "P: " << 1.0*P/100 << " flow: " << flow.size() << endl;
        cout << "pkt: " << sum << " smp_pkt: " << count << endl;
        double n = hll.Elem_getTotal(P);
        double em = 0;
        cout << n << endl;
        size_t am = 0;
        for (auto& cur : flow) {
            string f = cur.first;
            uint32_t ss = cur.second.size();
            vector<double> esi = hll.Elem_getParm((const uint32_t*)f.c_str(), n, 0);
            double nf = esi[32]*100/P;
            am = max(am, ss);
            em = max(em, nf);
            ouf << ss << " " << nf << endl;
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
            uint32_t e = 0;
            MurmurHash3_x86_32((const uint32_t*)pkt[i][1].c_str(), 32, 0, &e);
            if (e <= pow(2, 32)*P/100) {
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
    int p[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    string pkt = "e://desktop//res";
    string element = "e://desktop//res//element";
    string r = "e://desktop//res//base";
    for (int i = 9; i >= 0; i--) {
        //pkt_sampling_Throughout(1677722, 5, p[i]);
        pkt_sampling_driver(1677722, 5, p[i], pkt);
    }
    for (int i = 9; i >= 0; i--) {
        //element_sampling_Throughout(1677722, 5, p[i]);
        //element_sampling_driver(1677722, 5, p[i], element);
    }
    /*int b_ = 5;
    int e = 234546345;
    cout << e << endl;
    int p = e >> (32 - b_);
    int q = e - (p << (32 - b_));
    int leftmost_index = 0;
    while (q > 0) {
        leftmost_index++;
        q >>= 1;
    }
    leftmost_index = 32 - b_ - leftmost_index + 1;
    cout << p << " " << leftmost_index << endl;*/
    /*vHLL hll = vHLL(1677722, 5);
    const char *ip = "10.23.44.1";
    const char *d = "192.111.12.10";
    hll.insert(ip, d);
    double n = hll.getTotal();
    cout << n << endl;
    cout << hll.estimate(ip, n)[1] << endl;*/
}
