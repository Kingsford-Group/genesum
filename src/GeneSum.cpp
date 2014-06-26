#include <iostream>
#include <cstring>
#include <stdexcept> // gcc < 4.8 complains without this

#include "TranscriptGeneMap.hpp"
#include "gff.h"
#include "ezOptionParser.hpp"

void usage(ez::ezOptionParser& opt) {
    std::string usage;
    opt.getUsage(usage);
    std::cout << usage;
}


class ExpressionRecord {
    public:
        ExpressionRecord(const std::string& targetIn, uint32_t lengthIn,
                         std::vector<double>& expValsIn) :
            target(targetIn), length(lengthIn), expVals(expValsIn) {}

        ExpressionRecord( ExpressionRecord&& other ) {
            std::swap(target, other.target);
            length = other.length;
            std::swap(expVals, other.expVals);
        }

        ExpressionRecord(std::vector<std::string>& inputLine) {
            if (inputLine.size() < 3) {
                std::string err ("Any expression line must contain at least 3 tokens");
                throw std::invalid_argument(err);
            } else {
                auto it = inputLine.begin();
                target = *it; ++it;
                length = std::stoi(*it); ++it;
                for (; it != inputLine.end(); ++it) {
                    expVals.push_back(std::stod(*it));
                }
            }
        }

        std::string target;
        uint32_t length;
        std::vector<double> expVals;

};

// From : http://stackoverflow.com/questions/9435385/split-a-string-using-c11
std::vector<std::string> split(const std::string& str, int delimiter(int) = ::isspace){
    using namespace std;
    vector<string> result;
    auto e=str.end();
    auto i=str.begin();
    while (i != e) {
        i = find_if_not(i,e, delimiter);
        if (i == e) break;
        auto j = find_if(i,e, delimiter);
        result.push_back(string(i,j));
        i = j;
    }
    return result;
}

int main(int argc, const char* argv[]) {
    using namespace ez;
    using namespace std;

    ezOptionParser opt;
    opt.overview = "Summarize expression estimates at the gene level.";
    opt.syntax = "genesum [option]";
    opt.example = "genesum -g foo.gtf -e bar.txt -o baz.txt\n\n";
    opt.footer = "genesum v0.1.0 Copyright (C) 2014 Rob Patro\nThis program is free and without warranty.\n";

    opt.add( "", // Default.
            0, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Print help/usage message", // Help description.
            "-h",     // Flag token.
            "-help", // Flag token.
            "--help", // Flag token.
            "--usage" // Flag token.
           );

    opt.add("", 1, 1, 0, "annotation file", "-g", "--gtf");
    opt.add("", 1, 1, 0, "expression file", "-e", "--exp");
    opt.add("", 1, 1, 0, "output file", "-o", "--out");
    opt.add("gene_name", 0, 1, 0, "aggregation key", "-k", "--key");

    opt.parse(argc, argv);

    if (opt.isSet("-h")) {
        usage(opt);
        return 1;
    }

    vector<string> badOptions;
    uint32_t i;
    if(!opt.gotRequired(badOptions)) {
        for(i=0; i < badOptions.size(); ++i) {
            cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";
        }
        usage(opt);
        return 1;
    }

    string gtfFile;
    opt.get("-g")->getString(gtfFile);

    string expFname;
    opt.get("-e")->getString(expFname);
    ifstream expFile(expFname);

    string outFname;
    opt.get("-o")->getString(outFname);

    string key;
    opt.get("-k")->getString(key);

    cerr << "Aggregating estimates using key [" << key << "]\n";
    cerr << "Parsing GTF/GFF [" << gtfFile << "] . . .";
    auto tgm = TranscriptGeneMap::fromGTF(gtfFile, key);
    cerr << " done\n";
    cerr << "File contained: " << tgm.numGenes() << " genes and "
              << tgm.numTranscripts() << " transcripts\n";

    cerr << "Parsing input expression file\n";
    vector<string> comments;
    unordered_map<string, vector<ExpressionRecord>> geneExps;
    string l;
    size_t ln{0};

    while (getline(expFile, l)) {
        if (++ln % 1000 == 0) {
            std::cerr << "\r\rParsed " << ln << " expression lines";
        }
        auto it = find_if(l.begin(), l.end(),
                    [](char c) -> bool {return !isspace(c);});
        if (it != l.end()) {
            if (*it == '#') {
                comments.push_back(l);
            } else {
                vector<string> toks = split(l);
                ExpressionRecord er(toks);
                auto gn = tgm.geneName(er.target);
                geneExps[gn].push_back(move(er));
            }
        }
    }
    cerr << "\ndone\n";
    expFile.close();

    cerr << "Aggregating expressions to gene level . . .";
    ofstream outFile(outFname);

    // preserve any comments in the output
    for (auto& c : comments) {
        outFile << c << '\n';
    }

    for (auto& kv : geneExps) {
        auto& gn = kv.first;

        uint32_t geneLength{kv.second.front().length};
        vector<double> expVals(kv.second.front().expVals.size(), 0);
        const size_t NE{expVals.size()};

        for (auto& tranExp : kv.second) {
            geneLength = max(geneLength, tranExp.length);
            for (size_t i = 0; i < NE; ++i) { expVals[i] += tranExp.expVals[i]; }
        }

        outFile << gn << '\t' << geneLength;
        for (size_t i = 0; i < NE; ++i) {
            outFile << '\t' << expVals[i];
        }
        outFile << '\n';
    }

    outFile.close();
    cerr << " done\n";

    return 0;
}

