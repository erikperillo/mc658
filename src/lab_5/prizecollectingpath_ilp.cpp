/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/

#include "prizecollectingpath.h"
#include <limits>
#include <random>
#include <algorithm>
#include <ctime>

/*
 * Nome: Erik de Godoy Perillo
 * RA: 135582
 */

//infinite
#define INF (numeric_limits<double>::max()/2);

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

//type of pair arc and cost, calculated by function arc_cost
typedef pair<ListDigraph::Arc, double> arc_cost_pair;

//maximum number of iterations
#define HEUR_MAX_N_ITS 128
//maximum number of iterations allowed for best solution not to improve
#define HEUR_MAX_N_ITS_NOIMPROVE 4
//number of iterations of local search
#define HEUR_MAX_N_LOC_SEARCH_ITS 3

/*
 * Cost of arc (u, v) = prize(v) - cost((u, v))
 */
double arc_cost(
    const ListDigraph::Arc& a,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost)
{
    return prize[g.target(a)] - cost[a];
}

/*
 * Makes pairs (arc, cost) for an adjacency list of node v
 */
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

/*
 * Computes cost of a given path s-t
 */
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

/*
 * Comparator class of arcs and their costs calculated with arc_cost method
 */
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

/*
 * Makes random greedy choice of next arc for DFSjjj
 * It performs DFS from s to t and the next node to be searched is selected
 * by a random-greedy method.
 */
int rand_greedy_choice(
    const vector<arc_cost_pair>& arc_cost_pairs, int start_index)
{
    //random generators
    static uniform_real_distribution<double> unif(0.0, 1.0);
    static default_random_engine rand;

    //roulette method for selecting
    double cost_sum = 0.0;
    for(int i=start_index; i<static_cast<int>(arc_cost_pairs.size()); i++)
        cost_sum += arc_cost_pairs[i].second;

    if(cost_sum == 0)
        return start_index;

    double cum_prob = 0.0;
    double r = unif(rand);
    for(int i=start_index; i<static_cast<int>(arc_cost_pairs.size()); i++)
    {
        cum_prob += arc_cost_pairs[i].second/cost_sum;
        if(r <= cum_prob)
            return i;
    }

    return start_index;
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

vector<ListDigraph::Arc> _local_search(
    const vector<ListDigraph::Arc>& sol,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost)
{
    ListDigraph::NodeMap<bool> used(g);
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        used[v] = false;
    for(const auto& a: sol)
    {
        ListDigraph::Node u = g.source(a);
        ListDigraph::Node v = g.target(a);
        used[u] = true;
        used[v] = true;
    }

    vector<ListDigraph::Arc> path;

    for(const auto& a: sol)
    {
        double best_cost = arc_cost(a, g, prize, cost);
        ListDigraph::Arc fst = a;
        ListDigraph::Arc snd = a;
        ListDigraph::Node u = g.source(a);
        ListDigraph::Node v = g.target(a);

        for(ListDigraph::OutArcIt au(g, u); au!=INVALID; ++au)
        {
            ListDigraph::Node w = g.target(au);

            if(used[w])
                continue;

            for(ListDigraph::OutArcIt aw(g, w); aw!=INVALID; ++aw)
            {
                ListDigraph::Node x = g.target(aw);

                if(x != v)
                    continue;

                double cst = arc_cost(au, g, prize, cost)
                    + arc_cost(aw, g, prize, cost);
                if(cst > best_cost)
                {
                    fst = au;
                    snd = aw;
                    best_cost = cst;
                }
            }
        }

        path.push_back(fst);
        if(fst != snd)
        {
            ListDigraph::Node u = g.source(snd);
            used[u] = true;
            path.push_back(snd);
        }
    }

    return path;
}

vector<ListDigraph::Arc> local_search(
    const vector<ListDigraph::Arc>& sol,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost,
    int n_its)
{
    double best_cost = -INF;
    vector<ListDigraph::Arc> ls_sol = sol;

    for(int i=0; i<n_its; i++)
    {
        ls_sol = _local_search(ls_sol, g, prize, cost);

        double cst = path_cost(ls_sol, g, prize, cost);
        if(cst > best_cost)
            best_cost = cst;
        else
            break;
        cout << "\tlocal_search: " << cst << endl;
    }

    return ls_sol;
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
    //LB = pli_cutoff(g, prize, cost, s, t);
    double _ignore;
    vector<ListDigraph::Node> heur_path;
    prize_collecting_st_path_heuristic(
        g, prize, cost, s, t, heur_path, _ignore, LB, tMax);

    try
    {
    //setting cutoff for lower bound
    model.set(GRB_DoubleParam_Cutoff, LB*0.99);

    //optimizing
    model.optimize();

    if(model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        return NO_SOL;

    //setting bounds
    UB = model.get(GRB_DoubleAttr_ObjBoundC);

    //debug
    cout << "LB: " << LB << " | UB: " << UB << endl;

    path = node_path_from_pli_sol(x_a, g, s, t);

    //returning status
    if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        return OPT_SOL;
    else
        return HEUR_SOL;
    }
    catch(GRBException e)
    {
        cout << "EITA\n";
        return NO_SOL;
    }
}

vector<ListDigraph::Node> node_path_from_arc_path(
    const vector<ListDigraph::Arc>& arc_path,
    const ListDigraph& g)
{
    vector<ListDigraph::Node> path;

    for(int i=0; i<static_cast<int>(arc_path.size()); i++)
    {
        if(i == 0)
        {
            ListDigraph::Node u = g.source(arc_path[i]);
            path.push_back(u);
        }

        ListDigraph::Node v = g.target(arc_path[i]);
        path.push_back(v);
    }

    return path;
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
    vector<ListDigraph::Arc> best_sol;
    vector<ListDigraph::Arc> sol;
    double best_sol_cost = -INF;

    clock_t start_t = clock();

    int no_improve_count = 0;

    for(int i=0; i<HEUR_MAX_N_ITS; i++)
    {
        //random greedy initial solution
        sol = rand_greedy_sol(g, prize, cost, s, t);

        cout << "rand_greedy: " << path_cost(sol, g, prize, cost) << endl;
        //local search optimization
        sol = local_search(sol, g, prize, cost, HEUR_MAX_N_LOC_SEARCH_ITS);

        //updating best solution if there's one
        double sol_cost = path_cost(sol, g, prize, cost);
        if(sol_cost >= best_sol_cost)
        {
            best_sol_cost = sol_cost;
            best_sol = sol;
            no_improve_count = 0;
        }
        else
            no_improve_count++;

        //break if timeout
        if(double(clock() - start_t)/CLOCKS_PER_SEC > double(tMax))
        {
            cout << "TIMEOUT" << endl;
            break;
        }
        //break if best solution did not improve for certain number of its
        if(no_improve_count > HEUR_MAX_N_ITS_NOIMPROVE)
        {
            cout << "NO_IMPR_COUNT" << endl;
            break;
        }
    }

    //setting lower/upper bounds
    LB = 0.0;
    UB = best_sol_cost;

    //setting path
    path = node_path_from_arc_path(best_sol, g);

	return HEUR_SOL;
}
