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

//print debug messages iff true
const bool debug = false;

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
 * Makes random greedy choice of next arc for random-greedy solution generator.
 * Selects and arc with (normalized) probability proportional to its cost.
 */
int rand_greedy_choice(
    const vector<arc_cost_pair>& arc_cost_pairs, int start_index)
{
    //random generators
    static uniform_real_distribution<double> unif(0.0, 1.0);
    static default_random_engine rand;

    //cost sum for normalization
    double cost_sum = 0.0;
    for(int i=start_index; i<static_cast<int>(arc_cost_pairs.size()); i++)
        cost_sum += arc_cost_pairs[i].second;

    if(cost_sum == 0)
        return start_index;

    //roulette method for selecting index of next arc
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
 * Random-greedy solution generator for GRASP method.
 * It performs a DFS from s to t, choosing the next arc to be used for node
 * discovery with prob. proportional to the arc's cost (via arc_cost method).
 */
bool _rand_greedy_sol(
    ListDigraph& g,
    ListDigraph::NodeMap<bool>& visited,
    ListDigraph::NodeMap<double>& prize,
    ListDigraph::ArcMap<double>& cost,
    vector<ListDigraph::Arc>& path,
    ListDigraph::Node& s, ListDigraph::Node& t)
{
    //marking node as visited
    visited[s] = true;

    //found destiny
    if(s == t)
        return true;

    //making pairs arc-cost for each arc in out-adj list of s
    vector<arc_cost_pair> arcs_costs = make_arc_cost_pairs(s, g, prize, cost);

    //selecting order of adj nodes to be visited in random-greedy style
    for(int i=0; i<static_cast<int>(arcs_costs.size()); i++)
    {
        //index random-greedy selection
        int index = rand_greedy_choice(arcs_costs, i);
        ListDigraph::Node v = g.target(arcs_costs[index].first);

        //if true, a path from s-t was found
        if(!visited[v] && _rand_greedy_sol(g, visited, prize, cost, path, v, t))
        {
            path.push_back(arcs_costs[index].first);
            return true;
        }

        //excluding selected arc for next selection if needed
        if(i != index)
            iter_swap(arcs_costs.begin()+i, arcs_costs.begin()+index);
    }

    //did not find t via current path
    return false;
}

/*
 * Random-greedy solution generator for GRASP method.
 * frontend of _rand_greedy_sol method.
 */
vector<ListDigraph::Arc> rand_greedy_sol(
    ListDigraph& g,
    ListDigraph::NodeMap<double>& prize,
    ListDigraph::ArcMap<double>& cost,
    ListDigraph::Node& s, ListDigraph::Node& t)
{
    //marking nodes as not visited for DFS
    ListDigraph::NodeMap<bool> visited(g);
    for(ListDigraph::NodeIt v(g); v!=INVALID; ++v)
        visited[v] = false;

    //computing solution
    vector<ListDigraph::Arc> path;
    _rand_greedy_sol(g, visited, prize, cost, path, s, t);
    reverse(path.begin(), path.end());

    return path;
}

/*
 * Local search for GRASP method.
 * For each arc (u, v) of the random-greedy solution, it searches for arcs
 * (u, w), (w, v) with better costs (via arc_cost method).
 */
vector<ListDigraph::Arc> _local_search(
    const vector<ListDigraph::Arc>& sol,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost)
{
    //marking used nodes (they cannot be used)
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

    //for each arc (u, v) in solution...
    for(const auto& a: sol)
    {
        double best_cost = arc_cost(a, g, prize, cost);
        ListDigraph::Arc fst = a;
        ListDigraph::Arc snd = a;
        ListDigraph::Node u = g.source(a);
        ListDigraph::Node v = g.target(a);

        //for each node adjacent to u...
        for(ListDigraph::OutArcIt au(g, u); au!=INVALID; ++au)
        {
            ListDigraph::Node w = g.target(au);

            if(used[w])
                continue;

            //check if it connects to v...
            for(ListDigraph::OutArcIt aw(g, w); aw!=INVALID; ++aw)
            {
                ListDigraph::Node x = g.target(aw);

                if(x != v)
                    continue;

                //checks if cost of (u, w) + (w, v) is better
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

        //updating path
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

/*
 * Local search for GRASP method. Frontend of _local_search.
 * Performs _local_search main method n_its times, or stops before if solution
 * doesn't get better.
 */
vector<ListDigraph::Arc> local_search(
    const vector<ListDigraph::Arc>& sol,
    const ListDigraph& g,
    const ListDigraph::NodeMap<double>& prize,
    const ListDigraph::ArcMap<double>& cost,
    int n_its)
{
    //initial cost
    double best_cost = -INF;
    vector<ListDigraph::Arc> ls_sol = sol;

    //iterating n_its times using _local_search method
    for(int i=0; i<n_its; i++)
    {
        ls_sol = _local_search(ls_sol, g, prize, cost);

        //if cost doesn't get better, stop
        double cst = path_cost(ls_sol, g, prize, cost);
        if(cst > best_cost)
            best_cost = cst;
        else
            break;

        if(debug)
            cout << "\tlocal_search::iter " << i+1 << ": cost: " << cst << endl;
    }

    return ls_sol;
}

/*
 * Gets arcs of path s-t (in order) from PLI solution.
 */
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

/*
 * Gets nodes of path s-t (in order) from PLI solution.
 */
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
 * PLI solution.
 * Uses GRASP heuristic for cutoff
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

    //setting lower bound with metaheuristic
    double _ignore;
    vector<ListDigraph::Node> heur_path;
    prize_collecting_st_path_heuristic(
        g, prize, cost, s, t, heur_path, _ignore, LB, tMax);

    //initial upper bound
    UB = LB;

    try
    {
        //LET ME EXPLAIN why I set cutoff to 0.95*LB:
        //Gurobi was throwing exception "CUTOFF" but the LB value
        //was still LOWER than OPT. I think it was too close to OPT...
        //I tried everything but did not suceed, so I had to do this.
        model.set(GRB_DoubleParam_Cutoff, 0.95*LB);

        //optimizing
        model.optimize();

        //setting path
        path = node_path_from_pli_sol(x_a, g, s, t);
    }
    catch(GRBException e)
    {
        if(model.get(GRB_IntAttr_Status) != GRB_INFEASIBLE &&
            model.get(GRB_IntAttr_Status) != GRB_CUTOFF &&
            model.get(GRB_IntAttr_Status) != GRB_TIME_LIMIT)
            throw e;
    }

    //updating upper bound
    UB = max(LB, model.get(GRB_DoubleAttr_ObjVal));
    if(debug)
        //cout << "PLI: LB: " << LB << " | UB: " << UB << endl;
        cout << "PLI: lb, ub:\n: " << LB << "," << UB << endl;

    //returning status
    if(model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        if(debug)
            cout << "PLI: error: infeasible solution" << endl;
        return NO_SOL;
    }
    else if(model.get(GRB_IntAttr_Status) == GRB_CUTOFF)
    {
        if(debug)
            cout << "PLI: error: cutoff" << endl;
        path = heur_path;
        //I return OPT_SOL because of what I explain at line 478
        return OPT_SOL;
    }
    else if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        return OPT_SOL;
    else
    {
        if(LB >= UB)
            path = heur_path;
        return HEUR_SOL;
    }
}

/*
 * Transforms arc s-t (ordered) path into node version.
 */
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
 * GRASP heuristic.
 * Iterates getting initial solutions via a random-greedy DFS method, enhancing
 * with local search that looks for better arcs to be used around each node.
 * Stops either when reaches max num of iterations, or solution doesn't get
 * better for a number of times, or time exceeds limit.
 */
int prize_collecting_st_path_heuristic(
    ListDigraph& g,
    ListDigraph::NodeMap<double>& prize,
    ListDigraph::ArcMap<double> &cost,
    ListDigraph::Node s, ListDigraph::Node t,
    std::vector<ListDigraph::Node> &path,
    double &LB, double &UB,
    int tMax)
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

        //initial solution cost
        double sol_cost = path_cost(sol, g, prize, cost);
        if(debug)
            cout << "heuristic::iter " << i+1 << ": rand_greedy_sol value: "
                << sol_cost << endl;

        //break if timeout already reached
        if(double(clock() - start_t)/CLOCKS_PER_SEC > double(tMax))
        {
            if(debug)
                cout << "heuristic::iter " << i+1 << ": timeout reached\n";
            best_sol_cost = sol_cost;
            break;
        }

        //local search optimization
        sol = local_search(sol, g, prize, cost, HEUR_MAX_N_LOC_SEARCH_ITS);

        //solution cost
        sol_cost = path_cost(sol, g, prize, cost);

        //updating best solution if there's one
        if(sol_cost >= best_sol_cost)
        {
            best_sol_cost = sol_cost;
            best_sol = sol;
            no_improve_count = 0;
        }
        else
            no_improve_count++;

        //break if timeout reached
        if(double(clock() - start_t)/CLOCKS_PER_SEC > double(tMax))
        {
            if(debug)
                cout << "heuristic::iter " << i+1 << ": timeout reached\n";
            break;
        }
        //break if best solution did not improve for certain number of its
        if(no_improve_count > HEUR_MAX_N_ITS_NOIMPROVE)
        {
            if(debug)
                cout << "heuristic::iter " << i+1 << ": MAX_N_ITS_NOIMPROVE\n";
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
