//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#include "BuildPRG.h"


//read the MSA from a fasta file and put it in MSA
void BuildPRG::readMSAFromFastaFile(const std::string &filepath) {
  BOOST_LOG_TRIVIAL(info) << "Reading MSA from " << filepath << "...";
  auto lines = Utils::getVectorStringFromFile(filepath);
  for (const auto &line : lines) {
    if (line[0]!='>')
      MSA.push_back(line);
  }
  BOOST_LOG_TRIVIAL(info) << "Done!";
  MSA_NbColumns = MSA[0].size();
}

void BuildPRG::recursivelyBuildGraph(const SubAlignment &subAlignment, uint32_t parentVertexId) {
  BOOST_LOG_TRIVIAL(debug) << "@ BuildPRG::recursivelyBuildGraph: " << std::endl
                           << "subAlignment:" << std::endl
                           << subAlignment << std::endl;
  auto parentVertexDescriptor = add_vertex(graph);
  graph[parentVertexDescriptor].subalignments = subAlignment; //TODO: use move semantics here
  graph[parentVertexDescriptor].intervalType = NONMATCH;

  if (parentVertexId != NULL_VERTEX_ID) {
    //TODO: build the edge from parent to this node
  }

  //find the match/non-match regions and build the children of this node
  {
    //1. get the match and non-match regions
    std::vector <Interval> intervals = subAlignment.getIntervals(k);

    //2. create the nodes for each interval, and edges from the parent to these nodes
    std::vector<decltype(parentVertexDescriptor)> nonMatchChildDescriptors;
    for (const auto &interval : intervals) {
      //create the node
      auto childVertexDescriptor = add_vertex(graph);
      //configure new vertex
      graph[childVertexDescriptor].subalignments = SubAlignment(subAlignment.getSequencesNumbers(), interval.start,
                                                                interval.end, &MSA); //TODO: use move semantics here
      graph[childVertexDescriptor].intervalType = interval.intervalType;

      if (graph[childVertexDescriptor].intervalType == NONMATCH)
        nonMatchChildDescriptors.push_back(childVertexDescriptor);

      //create the edge
      add_edge(parentVertexDescriptor, childVertexDescriptor, graph); //TODO: check return value?
    }

    //3. for each non-match child, cluster the sequences and create non match children
    //TODO
  }
}