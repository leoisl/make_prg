//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#ifndef MAKE_PRG_PRGBOOSTGRAPH_H
#define MAKE_PRG_PRGBOOSTGRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <string>

//vertex properties
struct VertexInfo {
    std::string seq;
    size_t begin;
    size_t end;
};

//edge properties
struct EdgeInfo {
    std::string seq;
};

//we want to add the index properties to vertices and edges
typedef boost::property<boost::vertex_index_t, std::size_t, VertexInfo> VertexProps;
typedef boost::property<boost::edge_index_t, std::size_t, EdgeInfo> EdgeProps;



typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProps, EdgeProps> BoostGraph;

#endif //MAKE_PRG_PRGBOOSTGRAPH_H
