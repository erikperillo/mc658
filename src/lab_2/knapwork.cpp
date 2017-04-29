/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "knapsack.h"
#include <map>

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

bool promising(int k, int n, int d, int B,
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

#include <set>
void _bnb(
    int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, vector<int>& best, int& best_val, int t)
{
    bnb_counter++;

    if(k == n)
        return;

    int weight = get_weight(n-1, d, sol, w, c);
    int val = value(sol, p, n-1);

    if(weight > B)
        return;

    if(val > best_val)
    {
        best_val = val;
        best = sol;
        cout << "BEST VAL = " << best_val << " ON LVL = " << k << endl;
    }

    //int best_rem_val = (int)((B - weight)*(p[k]/float(w[k]+d)));
    //if(val + best_rem_val <= best_val)
     //   return;
    if(!promising(k, n, d, B, p, w, c, sol, best_val, weight, val))
        return;

    for(int pick=1; pick>=0; pick--)
    {
        sol[k] = pick;
        //prv(sol);
        _bnb(k+1, n, d, B, p, w, c, sol, best, best_val, t);
    }
}
int get_n_classes(const vector<int>& c)
{
    set<int> classes;
    for(int i=0; i<(int)c.size(); i++)
        classes.insert(c[i]);
    return (int)classes.size();
}


void _bnb2(
    int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, vector<int>& best, int& best_val, int t)
{
    //bnb_counter++;

    //_bnb(k, n, d, B, p, w, c, sol, best, best_val, t);
    /*_bnb(0, n/16, d, B, p, w, c, sol, best, best_val, t);
    for(int i=0; i<n; i++)
        sol[i] = 0;
    _bnb(0, n/8, d, B, p, w, c, sol, best, best_val, t);
    for(int i=0; i<n; i++)
        sol[i] = 0;*/
    /*_bnb(0, n/10, d, B, p, w, c, sol, best, best_val, t);
    for(int i=0; i<n; i++)
        sol[i] = 0;*/
    /*int n_classes = get_n_classes(c);
    int step = n_classes;
    for(int depth=1; depth<=32; depth+=1)
    {
        cout << "on depth = " << depth << endl;
        _bnb(0, depth, d, B, p, w, c, sol, best, best_val, t);
        for(int i=0; i<n; i++)
            sol[i] = 0;
    }
    cout << "depth = " << 64 << endl;
    _bnb(0, 64, d, B, p, w, c, sol, best, best_val, t);
        for(int i=0; i<n; i++)
            sol[i] = 0;
    cout << "depth = " << 256 << endl;
    _bnb(0, 256, d, B, p, w, c, sol, best, best_val, t);
        for(int i=0; i<n; i++)
            sol[i] = 0;
    cout << "FINAL BOSS" << endl;*/
    _bnb(0, n, d, B, p, w, c, sol, best, best_val, t);
}

vector<int> best_classes(int n,
    const vector<int>& c, const vector<int>& p, const vector<int>& w)
{
    set<int> classes = get_classes_set(c);
    map<int, double> rel_values;
    vector<double> best_val;
    vector<int> best_c;

    for(set<int>::iterator it = classes.begin(); it != classes.end(); ++it)
    {
        int clase = *it;
        int tot_val = 0;
        int tot_w = 0;
        for(int i=0; i<n; i++)
            if(c[i] == clase)
            {
                tot_val += p[i];
                tot_w += w[i];
            }
        double val = tot_val/(double)tot_w;

        best_val.push_back(val);
        best_c.push_back(clase);
        for(int i=(int)best_c.size()-1; i>0 && best_val[i] > best_val[i-1]; i--)
        {
            swap(&best_val[i], &best_val[i-1]);
            swap(&best_c[i], &best_c[i-1]);
        }
    }

    for(int i=0; i<(int)best_c.size(); i++)
        cout << "clase " << best_c[i] << ": " << best_val[i] << endl;

    return best_c;
}

vector<int> clase_index2(int n, vector<int>& c, vector<int>& best_clases)
{
    vector<int> new_indexes(n, 0);

    int pos = 0;
    for(int i=0; i<(int)best_clases.size(); i++)
    {
        int clase = best_clases[i];
        for(int j=0; j<n; j++)
        {
            if(c[j] == clase)
                new_indexes[pos++] = j;
        }
    }

    return new_indexes;
}

vector<int> clase_index(int n, const vector<int>& c)
{
    vector<bool> used(n, false);
    vector<int> new_index(n, 0);
    int n_classes = get_n_classes(c);
    set<int> classes;
    for(int i=0; i<n; i++)
        classes.insert(c[i]);

    for(int k=0; k<n; k++)
    {
        bool found = false;
        set<int>::iterator it = classes.begin();
        advance(it, k%n_classes);

        while(!found)
        {
            int clase = *it;
            for(int i=0; i<n; i++)
                if(!used[i] && c[i] == clase)
                {
                    found = true;
                    used[i] = true;
                    new_index[k] = i;
                    break;
                }
            advance(it, 1);
            if(it == classes.end())
                it = classes.begin();
            //cout << "clase = " << clase << endl;
        }
    }

    return new_index;
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
    int n_classes = get_n_classes(c);

    sort(sort_map.begin(), sort_map.end(), comparator);
    prv(p);
    prv(w);
    prv(c);
    p = map_to_indexes(p, sort_map);
    w = map_to_indexes(w, sort_map);
    c = map_to_indexes(c, sort_map);
    cout << "1" << endl;
    prv(p);
    prv(w);
    prv(c);
    //vector<int> best_clases = best_classes(n, c, p, w);
    /*vector<int> what = clase_index(n, c);
    p = map_to_indexes(p, what);
    w = map_to_indexes(w, what);
    c = map_to_indexes(c, what);
    cout << "2" << endl;
    prv(p);
    prv(w);
    prv(c);*/

    //_bnb(0, n, d, B, p, w, c, sol, best, best_val, t);
    _bnb2(0, n, d, B, p, w, c, sol, best, best_val, t);

    sol = best;

    cout << "counter: " << bnb_counter << " n: " << n << ", d: " << d
        << ", B: " << B << endl;
    int weight = get_weight(n-1, d, best, w, c);
    int val = value(best, p, n-1);
    cout << "weight: " << weight << ", val: " << val << endl;

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
