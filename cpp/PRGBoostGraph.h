//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#ifndef MAKE_PRG_PRGBOOSTGRAPH_H
#define MAKE_PRG_PRGBOOSTGRAPH_H

#include "includes.h"
#include "Subalignment.h"

//vertex properties
struct VertexInfo {
public:
    //the subalignment represented by this vertex
    SubAlignment subAlignment;
    VertexInfo() = default;

    /*
     * TODO: do we need this?
     */
    friend std::ostream& operator<<(std::ostream& os, const VertexInfo& vertexInfo) {
        os << std::endl << "---------- VertexInfo [BEGIN] ---------" << std::endl << vertexInfo.subAlignment <<
                           "---------- VertexInfo [END] -----------";
        return os;
    }
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

//typedefs to help
typedef boost::graph_traits<BoostGraph>::vertex_descriptor VertexDescriptor;
typedef boost::graph_traits<BoostGraph>::edge_descriptor EdgeDescriptor;

#endif //MAKE_PRG_PRGBOOSTGRAPH_H
