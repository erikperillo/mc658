/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "transportation.h"
#include <lemon/lp.h>

/*
 * Nome: Erik de Godoy Perillo
 * RA:   135582
 * Versao do gurobi usada: 7.0.2
 */

using namespace lemon;
using namespace std;

/*
 * Solving the problem.
 */
bool pl(
    ListBpGraph& g,                 //graph
    ListBpGraph::EdgeMap<int> &c,   //costs of edges
    ListBpGraph::NodeMap<int> &v,   //requirements/capacities of nodes
    ListBpGraph::EdgeMap<int> &sol, //solution
    int tMax)                       //maximum time
{
    //creating and setting up env
	GRBEnv env = GRBEnv();
    env.set(GRB_DoubleParam_TimeLimit, (double)tMax);

    //creating and setting up model
	GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "pl_gurobi");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    //variables of lp, mapping edges to values used in solution
    ListBpGraph::EdgeMap<GRBVar> x(g);

    //creating vars with cost
    int i = 0;
    for(ListBpGraph::EdgeIt e(g); e!=INVALID; ++e)
        x[e] = model.addVar(
            0, GRB_INFINITY, c[e], GRB_CONTINUOUS, "E_" + to_string(i++));

    //setting constraints
    for(ListBpGraph::NodeIt n(g); n!=INVALID; ++n)
    {
        GRBLinExpr flow = 0;
        bool in_t = false;

        for(ListBpGraph::IncEdgeIt e(g, n); e!=INVALID; ++e)
        {
            in_t = g.u(e) == n;
            flow += x[e];
        }

        if(in_t)
            model.addConstr(flow >= v[n]);
        else
            model.addConstr(flow <= v[n]);
    }

    //optimizing
    model.optimize();

    //writing solution to sol if found one
    if(model.get(GRB_IntAttr_SolCount) > 0)
        for(ListBpGraph::EdgeIt e(g); e!=INVALID; ++e)
            sol[e] = (int)x[e].get(GRB_DoubleAttr_X);

    //returning status
    return (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL);
}
