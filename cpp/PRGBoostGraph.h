//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#ifndef MAKE_PRG_PRGBOOSTGRAPH_H
#define MAKE_PRG_PRGBOOSTGRAPH_H

#include "includes.h"
#include "Subalignment.h"

#define NULL_VERTEX_ID UINT32_MAX

//vertex properties
struct VertexInfo {
public:
    //the subalignment represented by this vertex
    SubAlignment subalignments;
    IntervalType intervalType;
    VertexInfo() = default;
    //TODO: remove this, boost always use default constructor
    /*
    VertexInfo (const SubAlignment &subalignments, VertexType, vertexType) :
        subalignments{subalignments}, vertexType{vertexType} {}
    */
};

//edge properties
struct EdgeInfo {
    //empty for now
    EdgeInfo() = default;
};

//we want to add the index properties to vertices and edges
typedef boost::property<boost::vertex_index_t, uint32_t, VertexInfo> VertexProps;
typedef boost::property<boost::edge_index_t, uint32_t, EdgeInfo> EdgeProps;

//declare the Boost graph we want
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexProps, EdgeProps> BoostGraph;

#endif //MAKE_PRG_PRGBOOSTGRAPH_H
