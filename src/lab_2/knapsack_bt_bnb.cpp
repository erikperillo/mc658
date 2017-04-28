/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "knapsack.h"

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

bool promising(int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int best_val, int weight, int val)
{
    int rem_val = 0;
    int rem_weight = 0;
    int i;
    for(i=k; i<n; i++)
    {
        if(weight + rem_weight + w[i] > B)
            break;

        rem_val += p[i];
        rem_weight += w[i];
    }

    if(i < n)
        rem_val += (int)((p[i]/(float)w[i])*(B - (weight + rem_weight)));

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
    }

    /*int best_rem_val = (int)((B - weight)*(p[k]/(float)w[k]));
    if(val + best_rem_val <= best_val)
        return;*/
    if(!promising(k, n, d, B, p, w, c, sol, best_val, weight, val))
        return;

    for(int pick=1; pick>=0; pick--)
    {
        sol[k] = pick;
        //prv(sol);
        _bnb(k+1, n, d, B, p, w, c, sol, best, best_val, t);
    }
}

///
// Branch and Bound function
///
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

    sort(sort_map.begin(), sort_map.end(), comparator);
    p = map_to_indexes(p, sort_map);
    w = map_to_indexes(w, sort_map);
    c = map_to_indexes(c, sort_map);

    _bnb(0, n, d, B, p, w, c, sol, best, best_val, t);

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