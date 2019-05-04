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


  //lets just use uppercase bases, to simplify downstream process
  for (std::string &line : lines)
      boost::to_upper(line);

  MSA_NbColumns = MSA[0].size();
}

VertexDescriptor BuildPRG::createVertex(const SubAlignment &subAlignment) {
  VertexDescriptor vertexDescriptor = add_vertex(graph);
  graph[vertexDescriptor].subAlignment = subAlignment;
  BOOST_LOG_TRIVIAL(debug) << "Created vertex: " << graph[vertexDescriptor];
  return vertexDescriptor;
}

void BuildPRG::recursivelyBuildGraph(const SubAlignment &subAlignment, uint32_t parentVertexId) {
  BOOST_LOG_TRIVIAL(debug) << "@ BuildPRG::recursivelyBuildGraph: " << std::endl
                           << "subAlignment:" << std::endl
                           << subAlignment << std::endl;

  //create the vertex representing this subalignment
  VertexDescriptor parentVertexDescriptor = createVertex(subAlignment);

  if (parentVertexId != NULL_VERTEX_ID) {
    //TODO: build the edge from parent to this node
  }

  //find the match/non-match regions and build the children of this node
  {
    //1. get the match and non-match regions
    std::vector <Interval> intervals = subAlignment.getMatchAndNonMatchIntervals(k);

    //2. create the nodes for each interval, and edges from the parent to these nodes
    std::vector<VertexDescriptor> nonMatchChildDescriptors; //will store the vertices that are non-match
    for (const auto &interval : intervals) {
      //create the node
      VertexDescriptor childVertexDescriptor = createVertex(SubAlignment(subAlignment.getSequencesNumbers(), interval, &MSA));

      //record the non-match
      if (graph[childVertexDescriptor].subAlignment.getInterval().intervalType == NONMATCH)
        nonMatchChildDescriptors.push_back(childVertexDescriptor);

      //create the edge
      add_edge(parentVertexDescriptor, childVertexDescriptor, graph); //TODO: check return value?
    }

    //3. for each non-match child, cluster the sequences and create non match children
    //TODO
  }
}