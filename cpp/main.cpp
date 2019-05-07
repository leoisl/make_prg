//
// Created by Leandro Ishi Soares de Lima on 22/02/2019.
//

#include "includes.h"
#include "BuildPRG.h"

int main(void) {
    //initialize logger
    logging::core::get()->set_filter
            (
                    logging::trivial::severity >= logging::trivial::debug
            );

    //initialize random seed
    std::srand (std::time(NULL));


    /*
    //same tests from test_make_prg.py
    std::vector <std::string> files = {"../test/match.fa",
                                       "../test/nonmatch.fa",
                                       "../test/match.nonmatch.fa",
                                       "../test/nonmatch.match.fa",
                                       "../test/match.nonmatch.match.fa",
                                       "../test/shortmatch.nonmatch.match.fa",
                                       "../test/match.nonmatch.shortmatch.fa",
                                       "../test/match.staggereddash.fa",
                                       "../test/contains_n.fa",
                                       "../test/contains_RYKMSW.fa",
                                       "../test/contains_n_and_RYKMSW.fa",
                                       "../test/contains_n_and_RYKMSW_no_variants.fa"};
    std::vector <std::string> answers = {"ACGTGTTTTGTAACTGTGCCACACTCTCGAGACTGCATATGTGTC",
                                         " 5 AAACGTGGTT 6 CCCCCCCCCC 5 ",
                                         "AAACG 5 TGGTT 6 CCCCC 5 ",
                                         " 5 AAACGT 6 CCCCCC 5 GGTT",
                                         "AAACG 5 T 6 C 5 GGTT",
                                         " 5 AAACGT 6 ATTTTC 5 GGTT",
                                         "AAAC 5 GTGGTT 6 CCCCCT 5 ",
                                         "AAACGTGGTT",
                                         "AAACG 5 T 6 C 5 GGTT",
                                         "AAACG 5 T 6 C 5 GGTT",
                                         "AAACG 5 T 6 C 5 GGTT",
                                         "AAACGTGGTT"};


    for (int i = 0; i < files.size(); i++) {
      BuildPRG buildPRG(files[i]);

      auto prg = buildPRG.getPRG();
      std::cout << "[TEST " << i << "]: ";
      if (prg == answers[i])
        std::cout << "OK" << std::endl;
      else
        std::cout << "***FAILED***" << std::endl;

      std::cout << "File: " << files[i] << std::endl;
      std::cout << "COMPUTED PRG: " << prg << std::endl;
      std::cout << "CORRECT  PRG: " << answers[i] << std::endl;
      std::cout << "===================================================" << std::endl;
    }
     */


    BuildPRG buildPRG("../test/nested.fa");
}