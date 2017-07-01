/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/

#include "prizecollectingpath.h"

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

    //optimizing
    model.optimize();

    if(model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        return NO_SOL;

    //setting bounds
    LB = 0;
    UB = model.get(GRB_DoubleAttr_ObjBoundC);

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
