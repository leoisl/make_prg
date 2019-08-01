//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#include "BuildPRG.h"


//read the MSA from a fasta file and put it in MSA
void BuildPRG::readMSAFromFastaFile(const std::string &filepath) {
    BOOST_LOG_TRIVIAL(info) << "Reading MSA from " << filepath << "...";
    auto lines = Utils::getVectorStringFromFile(filepath);
    for (const auto &line : lines) {
        if (line[0] != '>')
            MSA.push_back(line);
    }
    BOOST_LOG_TRIVIAL(info) << "Done!";


    //lets just use uppercase bases, to simplify downstream process
    for (std::string &line : lines)
        boost::to_upper(line);

    MSA_NbColumns = MSA[0].size();
}

VertexDescriptor BuildPRG::createVertex(const SubAlignment &subAlignment, uint32_t nestingLevel) {
    VertexDescriptor vertexDescriptor = add_vertex(graph);
    graph[vertexDescriptor].subAlignment = subAlignment;
    graph[vertexDescriptor].nestingLevel = nestingLevel;
    BOOST_LOG_TRIVIAL(debug) << "Created vertex: " << graph[vertexDescriptor];
    return vertexDescriptor;
}

void BuildPRG::recursivelyBuildGraph(const SubAlignment &subAlignment, uint32_t nestingLevel) {
    BOOST_LOG_TRIVIAL(debug) << "@ BuildPRG::recursivelyBuildGraph: " << std::endl
                             << "subAlignment:" << std::endl
                             << subAlignment
                             << "Nesting level: " << nestingLevel << std::endl;

    //create the vertex representing this subalignment
    VertexDescriptor parentVertexDescriptor = createVertex(subAlignment, nestingLevel);

    //find the match/non-match regions and build the children of this node
    {
        //1. get the match and non-match regions
        //note that the intervals returned here can be MATCH, NONMATCH or TOO_SHORT
        std::vector<Interval> intervals = subAlignment.getMatchAndNonMatchIntervals(k);

        //2. create the nodes for each interval, and edges from the parent to these nodes
        std::vector<VertexDescriptor> nonMatchChildDescriptors; //will store the vertices that are non-match to cluster after
        for (const auto &interval : intervals) {
            //create the node
            VertexDescriptor childVertexDescriptor = createVertex(
                    SubAlignment(subAlignment.getSequencesNumbers(), interval, &MSA), nestingLevel); //it is still the same nesting level, as these are the match/non-match nodes (parent is unclassified)

            //record the non-match
            if (graph[childVertexDescriptor].subAlignment.getInterval().intervalType == NONMATCH)
                nonMatchChildDescriptors.push_back(childVertexDescriptor);

            //create the edge
            add_edge(parentVertexDescriptor, childVertexDescriptor, graph); //TODO: check return value?
        }

        //3. each non-match child defines a new site, but first we need to cluster
        for (const VertexDescriptor &nonMatchVertex : nonMatchChildDescriptors) {
            //3.1. checks if the max nesting level has been reached
            if (nestingLevel == maxNestingLevel) {
                //yes, we should not proceed
                //flag this
                graph[nonMatchVertex].subAlignment.getInterval().intervalType = IntervalType::NONMATCH_MAX_NESTING_LEVEL;
            }else {
                //ok, we can proceed
                //3.2 cluster the sequences in these non-match intervals
                SubAlignment::Clusters clusters = graph[nonMatchVertex].subAlignment.kMeansCluster(k);

                //recursion on each cluster
                for (const auto &subAlignment : clusters)
                    recursivelyBuildGraph(subAlignment, nestingLevel+1);
            }
        }
    }
}