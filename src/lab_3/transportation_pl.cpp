/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "transportation.h"
#include <lemon/lp.h>

/*
 * OBS: Fiz sozinho, entao fiz duas funcoes uma pro lemon e outra pro gurobi.
 * Nome: Erik de Godoy Perillo
 * RA:   135582
 */

using namespace lemon;
using namespace std;

//use lemon if true, use gurobi if false
const static bool use_lemon = true;

/*
 * Solving the problem using Gurobi.
 */
bool pl_gurobi(
    ListBpGraph& g,                 //graph
    ListBpGraph::EdgeMap<int> &c,   //costs of edges
    ListBpGraph::NodeMap<int> &v,   //requirements/capacities of nodes
    ListBpGraph::EdgeMap<int> &sol, //solution
    int tMax)                       //maximum time
{
    //creating env and model
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

    //setting up model/env
	env.set(GRB_DoubleParam_TimeLimit, tMax);
	model.set(GRB_StringAttr_ModelName, "pl_gurobi");
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    //creating vars with cost
    int i = 0;
    ListBpGraph::EdgeMap<GRBVar> x(g);
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

    //writing solution to sol if found one and returning
    if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        for(ListBpGraph::EdgeIt e(g); e!=INVALID; ++e)
            sol[e] = (int)x[e].get(GRB_DoubleAttr_X);
        return true;
    }
    else
        return false;
}

/*
 * Solving the problem using Lemon.
 */
bool pl_lemon(
    ListBpGraph& g,                 //graph
    ListBpGraph::EdgeMap<int> &c,   //costs of edges
    ListBpGraph::NodeMap<int> &v,   //requirements/capacities of nodes
    ListBpGraph::EdgeMap<int> &sol, //solution
    int tMax)                       //maximum time
{
    //model
    Lp lp;

    //objetive function
    Lp::Expr obj = 0;
    //vars with cost
    ListBpGraph::EdgeMap<Lp::Col> x(g);
    //creating vars and objective function
    for(ListBpGraph::EdgeIt e(g); e!=INVALID; ++e)
    {
        //setting up var
        x[e] = lp.addCol();
        lp.colLowerBound(x[e], 0);
        lp.colUpperBound(x[e], Lp::INF);
        //objective function
        obj += c[e]*x[e];
    }

    //setting constraints
    for(ListBpGraph::NodeIt n(g); n!=INVALID; ++n)
    {
        Lp::Expr flow = 0;
        bool in_t = false;

        for(ListBpGraph::IncEdgeIt e(g, n); e!=INVALID; ++e)
        {
            in_t = g.u(e) == n;
            flow += x[e];
        }

        if(in_t)
            lp.addRow(flow >= v[n]);
        else
            lp.addRow(flow <= v[n]);
    }

    //solving
    lp.min();
    lp.obj(obj);
    lp.solve();

    //writing solution to sol if found one and returning
    if(lp.primalType() == Lp::OPTIMAL)
    {
        for(ListBpGraph::EdgeIt e(g); e!=INVALID; ++e)
            sol[e] = (int)lp.primal(x[e]);
        return true;
    }
    else
        return false;
}

///
// PL function
///
bool pl(
    ListBpGraph& g,                 //graph
    ListBpGraph::EdgeMap<int> &c,   //costs of edges
    ListBpGraph::NodeMap<int> &v,   //requirements/capacities of nodes
    ListBpGraph::EdgeMap<int> &sol, //solution
    int tMax)                       //maximum time
{
    if(use_lemon)
        return pl_lemon(g, c, v, sol, tMax);
    else
        return pl_gurobi(g, c, v, sol, tMax);
}
