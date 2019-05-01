//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#ifndef MAKE_PRG_BUILDPRG_H
#define MAKE_PRG_BUILDPRG_H

#include "includes.h"
#include "PRGBoostGraph.h"
#include "Subalignment.h"
#include "Utils.h"

class BuildPRG {
private:
    std::vector<std::string> MSA;
    uint32_t MSA_NbColumns;
    uint32_t k;
    BoostGraph graph;
    std::string sep;
    std::string prg;

    //read the MSA from a fasta file and put it in MSA
    void readMSAFromFastaFile(const std::string &filepath);

    void recursivelyBuildGraph(const SubAlignment &subAlignment, uint32_t parentVertexId);

    //Build a graph representing this MSA
    //Common regions >= k compose the node of the graph
    //Arcs between the nodes are the sequences between these common regions
    inline void buildGraph() {
      recursivelyBuildGraph(SubAlignment(0, MSA.size(), 0, MSA_NbColumns, &MSA), NULL_VERTEX_ID);
    }

    VertexDescriptor createVertex(const SubAlignment &subAlignment, const IntervalType &intervalType);

public:
    BuildPRG(const std::string &filepath, uint32_t k=3, std::string sep=" ") : MSA{}, k{k}, graph{}, sep{sep} {
      readMSAFromFastaFile(filepath);
      buildGraph();
    }

    //build the PRG traversing the graph
    //Unsure if this is needed now
      /*
    std::string getPRG() const {
      auto currentAndEndNodeItPair = vertices(graph);
      decltype(currentAndEndNodeItPair.first) currentNodeIt, endNodeIt;
      std::tie(currentNodeIt, endNodeIt) = currentAndEndNodeItPair;
      bool startNodeIsDummy = itIsATrulyDummyStartNode();
      currentNodeIt+=int(startNodeIsDummy);
      bool endNodeIsDummy = itIsATrulyDummyEndNode();
      endNodeIt-=int(endNodeIsDummy);

      int varSiteMarker = 5;
      int alleleSiteMarker = 6;
      std::stringstream ss;

      //goes through all nodes, from left to right
      for (; currentNodeIt != endNodeIt; ++currentNodeIt ) {

        //should I explore this node?
        if (out_degree(*currentNodeIt, graph)>0 && (currentNodeIt!=endNodeIt-1 || (!endNodeIsDummy && currentNodeIt==endNodeIt-1))) {
          //yeah
          auto currentEdgeItAndEndEdgeIt = out_edges(*currentNodeIt, graph);
          decltype(currentEdgeItAndEndEdgeIt.first) currentEdgeIt, endEdgeIt;

          //print the node seq
          ss << graph[*currentNodeIt].seq << sep << varSiteMarker;
          bool printAlleleSiteMarker=false;

          //explore edge
          for (std::tie(currentEdgeIt, endEdgeIt) = currentEdgeItAndEndEdgeIt;
               currentEdgeIt<endEdgeIt; ++currentEdgeIt) {
            //print each out-edge as an allele
            if (printAlleleSiteMarker)
              ss << sep << alleleSiteMarker;
            else
              printAlleleSiteMarker = true;
            ss << sep << graph[*currentEdgeIt].seq;
          }
          ss << sep << varSiteMarker << sep;
          varSiteMarker+=2;
          alleleSiteMarker+=2;
        }

        //print the last one
        if (currentNodeIt == endNodeIt-1)
          ss << graph[*currentNodeIt].seq;
      }

      return ss.str();
    }
       */
};


#endif //MAKE_PRG_BUILDPRG_H
