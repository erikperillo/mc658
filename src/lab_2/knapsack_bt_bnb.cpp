/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/

///Preencher aqui para facilitar a correcao.
// Nome1: Erik de Godoy Perillo
// RA1: 135582
// Nome2:
// RA2:

#include "knapsack.h"
#include <iostream>
#include <algorithm>

using namespace std;

//debug flag
static bool debug = true;
//debug macros
#define str(x) #x
#define pr(x) ({if(debug) {cout << str(x)": " << x << endl;}})
#define prv(x) ({if(debug) {cout << str(x)": "; print_vec(x);}})

//counts iterations made in bnb
static int bnb_counter = 0;

///
// Bactracking function:
///
bool bt(int n, int d, int B, vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int t){
	return false;
}

//debug function
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

//comparator for sorting to be used in branch and bound
struct BBCompar
{
    BBCompar(const vector<int>& p, const vector<int>& w, const vector<int>& c):
        p(p), w(w), c(c) {;}

    bool operator()(size_t i, size_t j) const
    {
        return (p[i]/float(w[i])) > (p[j]/float(w[j]));
    }

    const vector<int>& p;
    const vector<int>& w;
    const vector<int>& c;
};

//sequence [start, end]
vector<int> seq(int start, int end)
{
    vector<int> sequence;
    for(int i=start; i<=end; i++)
        sequence.push_back(i);
    return sequence;
}

//swap
template <class T>
void swap(T* a, T* b)
{
    T aux = *a;
    *a = *b;
    *b = aux;
}

//gets new vector mapped to indexes
template<class T>
vector<T> map_to_indexes(const vector<T>& vec, const vector<int>& indexes)
{
    vector<T> mapped;
    for(int i=0; i<(int)vec.size(); i++)
        mapped.push_back(vec[indexes[i]]);
    return mapped;
}

//gets value of solution
int get_value(const vector<int>& sol, const vector<int>& p, int k)
{
    int val = 0;
    for(int i=0; i<=k; i++)
        val += sol[i]?p[i]:0;
    return val;
}

//gets weight of solution
int get_weight(int k, int d,
    const vector<int>& sol, const vector<int>& w, const vector<int>& c)
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

//returns true iff branch is promising
inline bool bnb_promising(int k, int n, int d, int B,
    const vector<int> &p, const vector<int> &w,
    int& weight, int& val,
    vector<int>& c_count, const vector<int>& c_hist, int& n_c,
    vector<int>& other_c, int& best_val)
{
    //counting classes used in/out current solution
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
    int j = 0;
    //using as much from already used classes as possible
    for(; (j+k)<n && j<used_c_rem; j++)
        if(weight + rem_weight + w[k+j] <= B)
        {
            rem_val += p[k+j];
            rem_weight += w[k+j];
        }
    //using from other classes
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

//branch and bound recursive method
void _bnb(
    int k, int n, int d, int B,
    const vector<int> &p, const vector<int> &w, const vector<int> &c,
    vector<int> &sol, int& weight, int& val,
    vector<int>& c_count, const vector<int>& c_hist, int& n_c,
    vector<int>& other_c,
    vector<int>& best, int& best_val, int t)
{
    bnb_counter++;

    //end of tree
    if(k == n)
        return;

    //too heavy to continue
    if(weight > B)
        return;

    //current solution is the best so far
    if(val > best_val)
    {
        best_val = val;
        best = sol;
        if(debug)
            cout << "BEST VAL = " << best_val << " ON LVL = " << k
                << " (weight = " << weight << ", get_weight() = "
                << get_weight(k, d, sol, w, c) << ")" << endl;
    }

    //checking if branch is promising
    if(!bnb_promising(k, n, d, B, p, w, weight, val,
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
bool bnb(
    int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, int t)
{
    bnb_counter = 0;
    vector<int> sort_map = seq(0, n-1);
    vector<int> best(n, 0);
    BBCompar comparator(p, w, c);
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

    //sorting by relative value
    sort(sort_map.begin(), sort_map.end(), comparator);
    p = map_to_indexes(p, sort_map);
    w = map_to_indexes(w, sort_map);
    c = map_to_indexes(c, sort_map);

    //calling recursive auxiliar method
    _bnb(0, n, d, B,
        p, w, c,
        sol, weight, val, c_count, c_hist, n_c, c_other,
        best, best_val, t);

    //setting solution
    sol = best;

    if(debug)
    {
        cout << "counter: " << bnb_counter << " n: " << n << ", d: " << d
            << ", B: " << B << endl;
        cout << "weight: " << get_weight(n-1, d, best, w, c)
            << ", val: " << get_value(sol, p, n-1) << endl;
    }

    return true;
}

/*brute force solution used for debug
static int bf_counter = 0;

void _bf(int k, int n, int d, int B,
    vector<int> &p, vector<int> &w, vector<int> &c,
    vector<int> &sol, vector<int>& best, int t)
{
    bf_counter++;

    if(k == n)
        return;

    int weight = get_weight(n-1, d, sol, w, c);
    int val = get_value(sol, p, n-1);
    int best_val = get_value(best, p, n-1);

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
    int val = get_value(best, p, n-1);
    cout << "weight: " << weight << ", val: " << val << endl;

    return true;
}*/
