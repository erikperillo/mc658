/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "knapsack.h"
#include <map>
#include <functional>

static bool debug = true;

#define str(x) #x
#define pr(x) ({if(debug) {cout << str(x)": " << x << endl;}})
#define prv(x) ({if(debug) {cout << str(x)": "; print_vec(x);}})

///Preencher aqui para facilitar a correcao.
// Nome1: Erik de Godoy Perillo
// RA1: 135582
// Nome2:
// RA2:
///
// Bactracking function:
///
bool bt(int n, int d, int B, vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int t){
	return false;
}

#include <iostream>
using namespace std;

template <class T>
void print_vec(vector<T> vec, bool nl=true)
{
    cout << "[";
    for(int i=0; i<(int(vec.size()))-1; i++)
        cout << vec[i] << ", ";
    if(vec.size() > 0)
        cout << vec[vec.size()-1];
    cout << "]";
    if(nl)
        cout << endl;
}


struct client
{
    client(int p, int w, int c, int d): p(p), w(w), c(c) {
        rel_p = p/float(w+d);
    }
    int p;
    int w;
    int c;
    float rel_p;
    void print()
    {
        cout << "client p: " << p << " | w: " << w << " | c: " << c
            << " | rel_p: " << rel_p << endl;
    }
};

bool gt(const client& a, const client& b)
{
    return a.rel_p > b.rel_p;
}
#include <algorithm>

struct compar
{
    compar(const vector<int>& p, const vector<int>& w, const vector<int>& c,
        int d): p(p), w(w), c(c), d(d) {;}

    bool operator()(size_t i, size_t j) const
    {
        return (p[i]/float(w[i])) > (p[j]/float(w[j]));
    }

    const vector<int>& p;
    const vector<int>& w;
    const vector<int>& c;
    int d;
};

vector<int> seq(int start, int end)
{
    vector<int> sequence;
    for(int i=start; i<=end; i++)
        sequence.push_back(i);
    return sequence;
}

template <class T>
void swap(T* a, T* b)
{
    T aux = *a;
    *a = *b;
    *b = aux;
}

template<class T>
vector<T> map_to_indexes(const vector<T>& vec, const vector<int>& indexes)
{
    vector<T> mapped;
    for(int i=0; i<(int)vec.size(); i++)
        mapped.push_back(vec[indexes[i]]);
    return mapped;
}

int value(const vector<int>& sol, const vector<int>& p, int k)
{
    int val = 0;
    for(int i=0; i<=k; i++)
        val += sol[i]?p[i]:0;
    return val;
}

int get_weight(int k, int d,
    const vector<int>& sol,
    const vector<int>& w, const vector<int>& c)
{
    int no_div_w = 0;
    set<int> classes;

    for(int i=0; i<=k; i++)
        if(sol[i])
        {
            classes.insert(c[i]);
            no_div_w += w[i];
        }

    return no_div_w + max(0, ((int)(classes.size())-1))*d;
}

/*template <class T>
T min(T a, T b)
{
    return (a > b)?a:b;
}*/

static int bnb_counter = 0;

set<int> get_classes_set(const vector<int>& c)
{
    set<int> classes;
    for(int i=0; i<(int)c.size(); i++)
        classes.insert(c[i]);
    return classes;
}

bool promising2(int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int best_val, int weight, int val)
{
    double best_rel_val = 0;
    for(int i=k; i<n; i++)
    {
        best_rel_val = max(best_rel_val, p[i]/(double)w[i]);
    }
    int add_val = (int)((B - weight)*best_rel_val);

    return (val + add_val) > best_val;
}

set<int> classes_in_sol(int k, const vector<int>& sol, const vector<int>& c)
{
    set<int> clases;

    for(int i=0; i<k; i++)
        if(sol[i])
            clases.insert(c[i]);

    return clases;
}

bool promising3(int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int best_val, int weight, int val)
{
    int rem_val = 0;
    int rem_w = 0;
    set<int> cl = classes_in_sol(n, sol, c);
    vector<float> rel_val(n-k, 0);
    for(int i=0; i<n-k; i++)
        if(cl.find(c[i]) == cl.end())
            rel_val[i] = p[i]/(float)(w[i]+d);
        else
            rel_val[i] = p[i]/(float)w[i];
    vector<bool> used(n-k, false);

    int a = 0;
    while(true)
    {
        bool all_used = true;
        for(int i=0; i<n-k; i++)
            if(!used[i])
            {
                all_used = false;
                break;
            }
        if(all_used)
            break;

        float max_rel_val = -1;
        int max_id = -1;

        for(int i=0; i<n-k; i++)
            if(!used[i] && rel_val[i] > max_rel_val)
            {
                max_rel_val = rel_val[i];
                max_id = i;
            }

        used[max_id] = true;

        cout << "a = " << a++ << endl;
        bool is_new = false;
        if(cl.find(c[max_id]) == cl.end())
        {
            is_new = true;
            if(weight + rem_w + w[max_id] + d > B)
                continue;
            cl.insert(c[max_id]);
            for(int j=0; j<n-k; j++)
                if(c[j] == c[max_id])
                    rel_val[j] = p[j]/(float)w[j];
        }
        else if(weight + rem_w + w[max_id] > B)
            continue;

        rem_w += w[max_id] + (is_new?d:0);
        rem_val += p[max_id];
    }

    return (val + rem_val) > best_val;
}

inline bool promising(int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int best_val, int weight, int val)
{
    int rem_val = 0;
    int rem_weight = 0;
    int i;
    /*set<int> cl;
    for(int i=0; i<n; i++)
        if(sol[i] == 1)
            cl.insert(c[i]);*/

    for(i=k; i<n; i++)
    {
        int wi = w[i];

        if(weight + rem_weight + wi <= B)
        {
            rem_val += p[i];
            rem_weight += wi;
            //cl.insert(c[i]);
        }
    }

    return (val + rem_val) > best_val;
}

inline bool promising(int k, int n, int d, int B,
    const vector<int> &p, const vector<int> &w,
    int& weight, int& val,
    vector<int>& c_count, const vector<int>& c_hist, int& n_c,
    vector<int>& other_c, int& best_val)
{
    int used_c_rem = 0;
    int n_others = 0;
    for(int i=0; i<(int)c_count.size(); i++)
    {
        if(c_count[i] > 0)
            used_c_rem += (c_hist[i] - c_count[i]);
        else if(c_hist[i] > 0)
        {
            other_c[n_others++] = c_hist[i];
            for(int j=n_others-1; j>0 && other_c[j] > other_c[j-1]; j--)
                swap(&other_c[j], &other_c[j-1]);
        }
    }

    int rem_val = 0;
    int rem_weight = 0;
    int j=0;
    for(; (j+k)<n && j<used_c_rem; j++)
        if(weight + rem_weight + w[k+j] <= B)
        {
            rem_val += p[k+j];
            rem_weight += w[k+j];
        }
    for(int m=0; m<n_others; m++)
    {
        rem_weight += (m == 0 && n_c == 0)?0:d;
        //rem_weight += (it == other_c.begin() && n_c == 0)?0:d;
        if(weight + rem_weight >= B)
            break;
        for(; (j+k)<n && j<other_c[m]; j++)
            if(weight + rem_weight + w[k+j] <= B)
            {
                rem_val += p[k+j];
                rem_weight += w[k+j];
            }
    }

    return val + rem_val > best_val;
}

void _bnb(
    int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int& weight, int& val,
    vector<int>& c_count, const vector<int>& c_hist, int& n_c,
    vector<int>& other_c,
    vector<int>& best, int& best_val, int t)
{
    bnb_counter++;

    if(k == n)
        return;

    if(weight > B)
        return;

    if(val > best_val)
    {
        best_val = val;
        best = sol;
        cout << "BEST VAL = " << best_val << " ON LVL = " << k
            << " (weight = " << weight << ", get_weight() = "
            << get_weight(k, d, sol, w, c) << ")" << endl;
    }

    if(!promising(k, n, d, B, p, w, weight, val,
        c_count, c_hist, n_c, other_c, best_val))
        return;

    //including node
    sol[k] = 1;
    weight += w[k] + ((c_count[c[k]] == 0 && n_c > 0)?d:0);
    n_c += (c_count[c[k]] == 0);
    c_count[c[k]] += 1;
    val += p[k];
    _bnb(k+1, n, d, B,
        p, w, c,
        sol, weight, val, c_count, c_hist, n_c, other_c,
        best, best_val, t);

    //excluding node
    sol[k] = 0;
    c_count[c[k]] -= 1;
    n_c -= (c_count[c[k]] == 0);
    weight -= w[k] + ((c_count[c[k]] == 0 && n_c > 0)?d:0);
    val -= p[k];
    _bnb(k+1, n, d, B,
        p, w, c,
        sol, weight, val, c_count, c_hist, n_c, other_c,
        best, best_val, t);
}

///
// Branch and Bound function
///
//
//
bool bnb(
    int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int t)
{
    bnb_counter = 0;
    vector<int> sort_map = seq(0, n-1);
    vector<int> best(n, 0);
    compar comparator(p, w, c, d);
    int best_val = 0;
    int weight = 0;
    int val = 0;
    int max_c = 0;
    int n_c = 0;
    for(int i=0; i<n; i++)
        if(c[i] > max_c)
            max_c = c[i];
    vector<int> c_count(max_c+2, 0);
    vector<int> c_hist(max_c+2, 0);
    vector<int> c_other(max_c+2, 0);
    for(int i=0; i<n; i++)
        c_hist[c[i]]++;

    sort(sort_map.begin(), sort_map.end(), comparator);
    p = map_to_indexes(p, sort_map);
    w = map_to_indexes(w, sort_map);
    c = map_to_indexes(c, sort_map);

    _bnb(0, n, d, B,
        p, w, c,
        sol, weight, val, c_count, c_hist, n_c, c_other,
        best, best_val, t);

    sol = best;

    cout << "counter: " << bnb_counter << " n: " << n << ", d: " << d
        << ", B: " << B << endl;
    cout << "weight: " << get_weight(n-1, d, best, w, c)
        << ", val: " << value(sol, p, n-1) << endl;

    return true;
}


template <class T>
T maxx(T a, T b)
{
    return (a > b)?a:b;
}

static int bf_counter = 0;

//brute force
void _bf(int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, vector<int>& best, int t)
{
    bf_counter++;

    if(k == n)
        return;

    int weight = get_weight(n-1, d, sol, w, c);
    int val = value(sol, p, n-1);
    int best_val = value(best, p, n-1);

    if(weight <= B && val > best_val)
        best = sol;

    for(int pick=0; pick<=1; pick++)
    {
        sol[k] = pick;
        //prv(sol);
        _bf(k+1, n, d, B, p, w, c, sol, best, t);
    }
}

bool bf(int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int t)
{
    bf_counter = 0;
    vector<int> best(n, 0);

    _bf(0, n, d, B, p, w, c, sol, best, t);

    sol = best;

    cout << "counter: " << bf_counter << " n: " << n << ", d: " << d
        << ", B: " << B << endl;
    int weight = get_weight(n-1, d, best, w, c);
    int val = value(best, p, n-1);
    cout << "weight: " << weight << ", val: " << val << endl;

    return true;
}
