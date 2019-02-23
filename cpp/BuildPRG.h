//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#ifndef MAKE_PRG_BUILDPRG_H
#define MAKE_PRG_BUILDPRG_H

#include <vector>
#include <string>
#include <iostream>
#include "PRGBoostGraph.h"
#include "Utils.h"
#include <sstream>

class BuildPRG {
private:
    std::vector<std::string> MSA;
    int k;
    std::string consensusString;
    BoostGraph graph;
    std::string sep;
    std::string prg;

    //read the MSA from a fasta file and put it in MSA
    void readMSAFromFastaFile(const std::string &filepath) {
      auto lines = Utils::getVectorStringFromFile(filepath);
      for (const auto &line : lines) {
        if (line[0]!='>')
          MSA.push_back(line);
      }
    }

    //builds the consensus string
    std::string buildConsensusString() {
      size_t nbColumns = MSA[0].size();
      for (size_t j=0; j<nbColumns; ++j) {
        char consensusBase=MSA[0][j];
        for (size_t i=1; i<MSA.size(); ++i) {
          //check
          if (consensusBase!=MSA[i][j]) { //inconsistency, we are done
            consensusBase = '?';
            break; //go to next column
          }
        }
        consensusString+=consensusBase;
      }
      //std::cout << "consensus: " << consensusString << std::endl;
      return consensusString;
    }

    //Build a graph representing this MSA
    //Common regions >= k compose the node of the graph
    //Arcs between the nodes are the sequences between these common regions
    void buildGraph() {
      //build first the nodes
      //let us add a dummy start node - it makes things easier
      {
        auto vertexDescriptor = add_vertex(graph);
        graph[vertexDescriptor].begin = graph[vertexDescriptor].end = 0;
      }

      //build the real nodes
      for (size_t i=0; i<consensusString.size(); ++i) {
        if (consensusString[i]!='?') {
          size_t j;
          for (j=i+1; j<consensusString.size() && consensusString[j]!='?'; ++j);
          if (j-i>=k){
            //std::cout << "new node: " << consensusString.substr(i, j-i) << std::endl;
            auto vertexDescriptor = add_vertex(graph);
            graph[vertexDescriptor].seq = consensusString.substr(i, j-i);
            graph[vertexDescriptor].begin = i;
            graph[vertexDescriptor].end = j;
          }
          i=j-1;
        }
      }
      //let us add a dummy end node - it makes things easier
      {
        auto vertexDescriptor = add_vertex(graph);
        graph[vertexDescriptor].begin = graph[vertexDescriptor].end = consensusString.size()+1;
      }

      //build the arcs now
      {
        auto currentAndEndNodeItPair = vertices(graph);
        decltype(currentAndEndNodeItPair.first) currentNodeIt, endNodeIt;
        std::tie(currentNodeIt, endNodeIt) = currentAndEndNodeItPair;

        auto lastNodeIt = currentNodeIt; //lastNode == first node
        ++currentNodeIt; //currentNode == second  node
        for (; currentNodeIt != endNodeIt; ++currentNodeIt ) {
          //check if we need to build an arc
          if (graph[*lastNodeIt].end <= graph[*currentNodeIt].begin) {
            //yeah - add one edge for each sequence alignment
            for (const auto &SA : MSA) {
              auto edgeDescriptor = add_edge(*lastNodeIt, *currentNodeIt, graph).first;
              auto edgeSeq = SA.substr(graph[*lastNodeIt].end, graph[*currentNodeIt].begin-graph[*lastNodeIt].end);
              //std::cout << "[edgeSeq] = " << edgeSeq << std::endl;
              graph[edgeDescriptor].seq = edgeSeq;
            }
          }
          lastNodeIt = currentNodeIt;
        }
      }
    }

    //returns true if the start node is truly a dummy node, or it has a purpose
    bool itIsATrulyDummyStartNode() const {
      auto currentEdgeItAndEndEdgeIt = out_edges(0, graph);
      decltype(currentEdgeItAndEndEdgeIt.first) currentEdgeIt, endEdgeIt;
      for (std::tie(currentEdgeIt, endEdgeIt) = currentEdgeItAndEndEdgeIt;
           currentEdgeIt<endEdgeIt; ++currentEdgeIt) {
        if (graph[*currentEdgeIt].seq != "") {
          return false;
        }
      }
      return true;
    }

    //returns true if the end node is truly a dummy node, or it has a purpose
    bool itIsATrulyDummyEndNode() const {
      auto currentEdgeItAndEndEdgeIt = in_edges(int(num_vertices(graph))-1, graph);
      decltype(currentEdgeItAndEndEdgeIt.first) currentEdgeIt, endEdgeIt;
      for (std::tie(currentEdgeIt, endEdgeIt) = currentEdgeItAndEndEdgeIt;
           currentEdgeIt<endEdgeIt; ++currentEdgeIt) {
        if (graph[*currentEdgeIt].seq != "") {
          return false;
        }
      }
      return true;
    }



public:
    BuildPRG(const std::string &filepath, int k=3, std::string sep=" ") : MSA(), k(k), consensusString(), graph(), sep(sep) {
      readMSAFromFastaFile(filepath);
      buildConsensusString();
      buildGraph();
    }

    //build the PRG traversing the graph
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
};


#endif //MAKE_PRG_BUILDPRG_H
