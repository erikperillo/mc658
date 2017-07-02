/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/

#include "prizecollectingpath.h"
#include <limits>
#include <random>
#include <algorithm>

/*
 * Nome: Erik de Godoy Perillo
 * RA: 135582
 */

//return codes
enum
{
    //no solution found
    NO_SOL=0,
    //optimal solution found
    OPT_SOL=1,
    //heuristic solution found
    HEUR_SOL=2
};

typedef pair<ListDigraph::Arc, double> arc_cost_pair;

inline double arc_cost(
    const ListDigraph::Arc& a,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost)
{
    return prize[g.target(a)] - cost[a];
}

vector<arc_cost_pair> make_arc_cost_pairs(
    const ListDigraph::Node& v,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost)
{
    vector<arc_cost_pair> cands;
    for(ListDigraph::OutArcIt a(g, v); a!=INVALID; ++a)
        cands.push_back(make_pair(a, arc_cost(a, g, prize, cost)));
    return cands;
}

class ArcCompar
{
    public:
    ArcCompar(
        const ListDigraph& g,
        const ListDigraph::NodeMap<double>& prize,
        const ListDigraph::ArcMap<double>& cost):
        g(g), prize(prize), cost(cost)
    {;}

    bool operator()(const ListDigraph::Arc& a1, const ListDigraph::Arc& a2)
    {
        return arc_cost(a1, g, prize, cost) > arc_cost(a2, g, prize, cost);
    }

    private:
    const ListDigraph& g;
    const ListDigraph::NodeMap<double>& prize;
    const ListDigraph::ArcMap<double>& cost;
};

vector<ListDigraph::Arc> candidates_list(
    ListDigraph& g,
    ListDigraph::NodeMap<double>& prize, //prizes on nodes
    ListDigraph::ArcMap<double>& cost, //costs of arcs
    ListDigraph::Node v,
    int k)
{
    vector<ListDigraph::Arc> cands;
    for(ListDigraph::OutArcIt a(g, v); a!=INVALID; ++a)
        cands.push_back(a);
    sort(cands.begin(), cands.end(), ArcCompar(g, prize, cost));
    return vector<ListDigraph::Arc>(cands.begin(), cands.begin()+k);
}

int rand_greedy_choice(
    const vector<arc_cost_pair>& arc_cost_pairs, int start_index)
{
    static uniform_real_distribution<double> unif(0.0, 1.0);
    static default_random_engine rand;

    double cost_sum = 0.0;
    for(int i=start_index; i<static_cast<int>(arc_cost_pairs.size()); i++)
        cost_sum += arc_cost_pairs[i].second;

    double cum_prob = 0.0;
    double r = unif(rand);
    for(int i=start_index; i<static_cast<int>(arc_cost_pairs.size()); i++)
    {
        cum_prob += arc_cost_pairs[i].second/cost_sum;
        if(r <= cum_prob)
            return i;
    }

    return -1;
}

/*
 * Simple heuristic for lower-bounding PLI: DFS from s to t
 */
bool _rand_greedy_sol(
    ListDigraph& g, //directed graph
    ListDigraph::NodeMap<bool>& visited, //prizes on nodes
    ListDigraph::NodeMap<double>& prize, //prizes on nodes
    ListDigraph::ArcMap<double>& cost, //costs of arcs
    vector<ListDigraph::Arc>& path,
    ListDigraph::Node& s, ListDigraph::Node& t) //source/dest nodes
{
    //marking node as visited
    visited[s] = true;

    //found destiny
    if(s == t)
        return true;

    //making pairs arc-cost for each arc in out-adj list of s
    vector<arc_cost_pair> arcs_costs = make_arc_cost_pairs(s, g, prize, cost);

    //selecting order of adj nodes to be visited in random-greedy style:
    //an arc a is selected with (normalized) probability arc_cost(a)
    for(int i=0; i<static_cast<int>(arcs_costs.size()); i++)
    {
        int index = rand_greedy_choice(arcs_costs, i);
        ListDigraph::Node v = g.target(arcs_costs[index].first);

        //if true, a path from s-t was found
        if(!visited[v] && _rand_greedy_sol(g, visited, prize, cost, path, v, t))
        {
            path.push_back(arcs_costs[index].first);
            return true;
        }

        if(i != index)
            iter_swap(arcs_costs.begin()+i, arcs_costs.begin()+index);
    }

    return false;
}

vector<ListDigraph::Arc> rand_greedy_sol(
    ListDigraph& g,
    ListDigraph::NodeMap<double>& prize, //prizes on nodes
    ListDigraph::ArcMap<double>& cost, //costs of arcs
    ListDigraph::Node& s, ListDigraph::Node& t) //source/dest nodes
{
    ListDigraph::NodeMap<bool> visited(g);
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        visited[v] = false;

    vector<ListDigraph::Arc> path;
    _rand_greedy_sol(g, visited, prize, cost, path, s, t);
    reverse(path.begin(), path.end());

    return path;
}

double path_cost(
    const vector<ListDigraph::Arc>& path,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost)
{
    double p_cost = 0.0;

    //first node s
    if(path.size() > 0)
        p_cost += prize[g.source(path[0])];
    //adding costs including t
    for(const auto& a: path)
        p_cost += arc_cost(a, g, prize, cost);

    return p_cost;
}

/*void _dfs(
    ListDigraph& g, //directed graph
    ListDigraph::NodeMap<bool>& visited, //prizes on nodes
    ListDigraph::ArcMap<GRBVar>& x_a,
    ListDigraph::NodeMap<GRBVar>& x_v,
    ListDigraph::Node& s, ListDigraph::Node& t) //source/dest nodes
{
    //marking node as visited
    visited[s] = true;

    for(ListDigraph::OutArcIt a(g, s); a!=INVALID; ++a)
    {
        ListDigraph::Node v = g.target(a);
        if(!visited[v] && x_v[v].get(GRB_DoubleAttr_X) > 0.5
            && x_a[a].get(GRB_DoubleAttr_X) > 0.5)
            _dfs(g, visited, x_a, x_v, v, t);
    }
}

bool dfs(ListDigraph& g,
    ListDigraph::ArcMap<GRBVar>& x_a,
    ListDigraph::NodeMap<GRBVar>& x_v,
    ListDigraph::Node& s, ListDigraph::Node& t)
{
    ListDigraph::NodeMap<bool> visited(g);
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        visited[v] = false;
    _dfs(g, visited, x_a, x_v, s, t);
    int c = 0;
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        c += visited[v];
    cout << "c = " << c << endl;
    return visited[t];
}*/

vector<ListDigraph::Arc> arc_path_from_pli_sol(
    ListDigraph::ArcMap<GRBVar>& x_a,
    ListDigraph& g,
    ListDigraph::Node s, ListDigraph::Node t)
{
    vector<ListDigraph::Arc> path;
    ListDigraph::Node v = s;

    while(v != t)
        for(ListDigraph::OutArcIt a(g, v); a!=INVALID; ++a)
            if(x_a[a].get(GRB_DoubleAttr_X) > 0.5)
            {
                v = g.target(a);
                path.push_back(a);
                break;
            }

    reverse(path.begin(), path.end());
    return path;
}

vector<ListDigraph::Node> node_path_from_pli_sol(
    ListDigraph::ArcMap<GRBVar>& x_a,
    ListDigraph& g,
    ListDigraph::Node s, ListDigraph::Node t)
{
    vector<ListDigraph::Node> path;
    ListDigraph::Node v = s;

    path.push_back(v);
    while(v != t)
        for(ListDigraph::OutArcIt a(g, v); a!=INVALID; ++a)
            if(x_a[a].get(GRB_DoubleAttr_X) > 0.5)
            {
                v = g.target(a);
                path.push_back(v);
                break;
            }

    reverse(path.begin(), path.end());
    return path;
}

/*
 * Simple heuristic for lower-bounding PLI: DFS from s to t
 */
double _pli_cutoff(
    ListDigraph& g, //directed graph
    ListDigraph::NodeMap<bool>& visited, //prizes on nodes
    ListDigraph::NodeMap<double>& prize, //prizes on nodes
    ListDigraph::ArcMap<double>& cost, //costs of arcs
    ListDigraph::Node& s, ListDigraph::Node& t) //source/dest nodes
{
    visited[s] = true;

    if(s == t)
        return prize[t];

    double lb = -(numeric_limits<double>::max()/2);
    for(ListDigraph::OutArcIt a(g, s); a!=INVALID; ++a)
    {
        ListDigraph::Node v = g.target(a);
        if(!visited[v])
        {
            double v_lb = _pli_cutoff(g, visited, prize, cost, v, t);
            lb = max(lb, v_lb - cost[a] + prize[v]);
        }
    }

    return lb;
}

double pli_cutoff(
    ListDigraph& g, //directed graph
    ListDigraph::NodeMap<double>& prize, //prizes on nodes
    ListDigraph::ArcMap<double>& cost, //costs of arcs
    ListDigraph::Node& s, ListDigraph::Node& t) //source/dest nodes
{
    ListDigraph::NodeMap<bool> visited(g);
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        visited[v] = false;
    return max(0.0, _pli_cutoff(g, visited, prize, cost, s, t));
}

/*
 * PLI solution
 */
int prize_collecting_st_path_pli(
    ListDigraph& g, //directed graph
    ListDigraph::NodeMap<double>& prize, //prizes on nodes
    ListDigraph::ArcMap<double>& cost, //costs of arcs
    ListDigraph::Node s, ListDigraph::Node t, //source/dest nodes
    std::vector<ListDigraph::Node> &path, //path from s to t
    double &LB, double &UB, //lower/upper bounds
    int tMax) //time limit
{
    //creating and setting up env
	GRBEnv env = GRBEnv();
    env.set(GRB_DoubleParam_TimeLimit, (double)tMax);

    //creating and setting up model
	GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "prize_collecting_path_pli");

    //variables of lp representing arcs used
    ListDigraph::ArcMap<GRBVar> x_a(g);
    //variables of lp representing nodes used
    ListDigraph::NodeMap<GRBVar> x_v(g);

    int i = 0;
    //creating binary vars with costs of arcs
    for(ListDigraph::ArcIt a(g); a!=INVALID; ++a)
        x_a[a] = model.addVar(
            0, 1, 0, GRB_BINARY, "a_" + to_string(i++));

    //creating binary vars with prizes of nodes
    i = 0;
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        x_v[v] = model.addVar(
            0, 1, 0, GRB_BINARY, "v_" + to_string(i++));

    //constraint: only one arc exits source
    GRBLinExpr sum = 0;
    for(ListDigraph::OutArcIt a(g, s); a!=INVALID; ++a)
        sum += x_a[a];
    model.addConstr(sum == 1);

    //constraint: zero arcs enter source
    sum = 0;
    for(ListDigraph::InArcIt a(g, s); a!=INVALID; ++a)
        sum += x_a[a];
    model.addConstr(sum == 0);

    //constraint: only one arcs enters destiny
    sum = 0;
    for(ListDigraph::InArcIt a(g, t); a!=INVALID; ++a)
        sum += x_a[a];
    model.addConstr(sum == 1);

    //constraint: zero arcs exit destiny
    sum = 0;
    for(ListDigraph::OutArcIt a(g, t); a!=INVALID; ++a)
        sum += x_a[a];
    model.addConstr(sum == 0);

    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        if(v != s && v != t)
        {
            //constraint: only one arc leaves selected nodes
            sum = 0;
            for(ListDigraph::OutArcIt a(g, v); a!=INVALID; ++a)
                sum += x_a[a];
            model.addConstr(sum == x_v[v]);

            //constraint: only one arc enters selected nodes
            sum = 0;
            for(ListDigraph::InArcIt a(g, v); a!=INVALID; ++a)
                sum += x_a[a];
            model.addConstr(sum == x_v[v]);
        }

    //objetive function to be maximized
    GRBLinExpr obj = 0;
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        obj += x_v[v]*prize[v];
    for(ListDigraph::ArcIt a(g); a!=INVALID; ++a)
        obj -= x_a[a]*cost[a];
    model.setObjective(obj, GRB_MAXIMIZE);

    //setting lower bound with simple heuristic
    LB = pli_cutoff(g, prize, cost, s, t);
    //setting cutoff for lower bound
    model.set(GRB_DoubleParam_Cutoff, LB);

    //optimizing
    model.optimize();

    if(model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        return NO_SOL;

    //setting bounds
    UB = model.get(GRB_DoubleAttr_ObjBoundC);

    //debug
    cout << "LB: " << LB << " | UB: " << UB << endl;
    auto heur_path = rand_greedy_sol(g, prize, cost, s, t);
    double heur_path_cost = path_cost(heur_path, g, prize, cost);
    auto arc_path = arc_path_from_pli_sol(x_a, g, s, t);
    double p_cost = path_cost(arc_path, g, prize, cost);
    cout << "heur LB: " << heur_path_cost << " | path UB: " << p_cost << endl;

    path = node_path_from_pli_sol(x_a, g, s, t);

    //returning status
    if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        return OPT_SOL;
    else
        return HEUR_SOL;
}

/*
 * Heuristic
 */
int prize_collecting_st_path_heuristic(
    ListDigraph& g, //directed graph
    ListDigraph::NodeMap<double>& prize, //prizes on nodes
    ListDigraph::ArcMap<double> &cost, //costs of arcs
    ListDigraph::Node s, ListDigraph::Node t, //source/dest nodes
    std::vector<ListDigraph::Node> &path, //path from s to t
    double &LB, double &UB, //lower/upper bounds
    int tMax) //time limit
{
	return 0;
}
