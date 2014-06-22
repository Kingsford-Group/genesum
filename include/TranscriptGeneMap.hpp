/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#ifndef TRANSCRIPT_GENE_MAP_HPP
#define TRANSCRIPT_GENE_MAP_HPP

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstring>

#include "gff.h"

class TranscriptGeneMap {
using Index = size_t;
using Size = size_t;
using NameVector = std::vector<std::string>;
using IndexVector = std::vector<size_t>;
using IndexVectorList = std::vector<std::vector<size_t>>;

private:
    NameVector _transcriptNames;
    NameVector _geneNames;
    IndexVector _transcriptsToGenes;
    IndexVectorList _genesToTranscripts;
    bool _haveReverseMap;

    void _computeReverseMap() {

        _genesToTranscripts.resize( _geneNames.size(), {});

        Index transcriptID = 0;
        size_t maxNumTrans = 0;
        Index maxGene;
        for ( transcriptID = 0; transcriptID < _transcriptsToGenes.size(); ++transcriptID ) {
            _genesToTranscripts[ _transcriptsToGenes[transcriptID] ].push_back( transcriptID );
            if ( maxNumTrans < _genesToTranscripts[ _transcriptsToGenes[transcriptID] ].size() ) {
                maxNumTrans = _genesToTranscripts[ _transcriptsToGenes[transcriptID] ].size();
                maxGene = _transcriptsToGenes[transcriptID];
            }
        }
        std::cerr << "max # of transcripts in a gene was " << maxNumTrans << " in gene " << _geneNames[maxGene] << "\n";
    }

    /*
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & _transcriptNames;
        ar & _geneNames;
        ar & _transcriptsToGenes;
        ar & _genesToTranscripts;
        ar & _haveReverseMap;
    }
    */

public:
    TranscriptGeneMap() :
        _transcriptNames(NameVector()), _geneNames(NameVector()),
        _transcriptsToGenes(IndexVector()), _haveReverseMap(false) {}


    TranscriptGeneMap( const NameVector &transcriptNames,
                       const NameVector &geneNames,
                       const IndexVector &transcriptsToGenes ) :
        _transcriptNames(transcriptNames), _geneNames(geneNames),
        _transcriptsToGenes(transcriptsToGenes), _haveReverseMap(false) {}


    TranscriptGeneMap(const TranscriptGeneMap& other) = default;
    TranscriptGeneMap& operator=(const TranscriptGeneMap& other) = default;


    static TranscriptGeneMap fromGTF(const std::string& fname, std::string key) {

        using std::unordered_set;
        using std::unordered_map;
        using std::vector;
        using std::tuple;
        using std::string;
        using std::get;

        // Use GffReader to read the file
        GffReader reader(const_cast<char*>(fname.c_str()));
        // Remember the optional attributes
        reader.readAll(true);

        struct TranscriptKeyPair {
            const char* transcript_id;
            const char* key;
            TranscriptKeyPair(const char* t, const char* k) :
                transcript_id(t), key(k) {}
        };

        // The user can group transcripts by gene_id, gene_name, or
        // an optinal attribute that they provide as a string.
        enum class TranscriptKey { GENE_ID, GENE_NAME, DYNAMIC };

        // Select the proper attribute by which to group
        TranscriptKey tkey = TranscriptKey::GENE_ID;

        if (key == "gene_id") {
        } else if (key == "gene_name") {
            tkey = TranscriptKey::GENE_NAME;
        } else {
            tkey = TranscriptKey::DYNAMIC;
        }

        // Iterate over all transcript features and build the
        // transcript <-> key vector.
        auto nfeat = reader.gflst.Count();
        std::vector<TranscriptKeyPair> feats;
        for (int i=0; i < nfeat; ++i) {
            auto f = reader.gflst[i];
            if (f->isTranscript()) {
                const char* keyStr;
                switch (tkey) {
                    case TranscriptKey::GENE_ID:
                        keyStr = f->getGeneID();
                        break;
                    case TranscriptKey::GENE_NAME:
                        keyStr = f->getGeneName();
                        break;
                    case TranscriptKey::DYNAMIC:
                        keyStr = f->getAttr(key.c_str());
                        break;
                }
                feats.emplace_back(f->getID(), keyStr);
            }
        }

        // Given the transcript <-> key vector, build the
        // TranscriptGeneMap.

        IndexVector t2g;
        NameVector transcriptNames;
        NameVector geneNames;

        // holds the mapping from transcript ID to gene ID
        IndexVector t2gUnordered;
        // holds the set of gene IDs
        unordered_map<string, size_t> geneNameToID;

        // To read the input and assign ids
        size_t geneCounter = 0;
        string transcript;
        string gene;

        std::sort( feats.begin(), feats.end(),
                []( const TranscriptKeyPair & a, const TranscriptKeyPair & b) -> bool {
                return std::strcmp(a.transcript_id, b.transcript_id) < 0;
                } );

        std::string currentTranscript = "";
        for ( auto & feat : feats ) {

            std::string gene(feat.key);
            std::string transcript(feat.transcript_id);

            if ( transcript != currentTranscript ) {
                auto geneIt = geneNameToID.find(gene);
                size_t geneID = 0;

                if ( geneIt == geneNameToID.end() ) {
                    // If we haven't seen this gene yet, give it a new ID
                    geneNameToID[gene] = geneCounter;
                    geneID = geneCounter;
                    geneNames.push_back(gene);
                    ++geneCounter;
                } else {
                    // Otherwise lookup the ID
                    geneID = geneIt->second;
                }

                transcriptNames.push_back(transcript);
                t2g.push_back(geneID);

                //++transcriptID;
                currentTranscript = transcript;
            }

        }

        return TranscriptGeneMap(transcriptNames, geneNames, t2g);

    }



    Index INVALID { std::numeric_limits<Index>::max() };

    Index findTranscriptID( const std::string &tname ) {
        using std::distance;
        using std::lower_bound;
        auto it = lower_bound( _transcriptNames.begin(), _transcriptNames.end(), tname );
        return ( it == _transcriptNames.end() ) ? INVALID : ( distance(_transcriptNames.begin(), it) );
    }

    Size numTranscripts() {
        return _transcriptNames.size();
    }
    Size numGenes() {
        return _geneNames.size();
    }

    bool needReverse() {
        if ( _haveReverseMap ) {
            return false;
        } else {
            _computeReverseMap();
            return true;
        }
    }

    const IndexVector &transcriptsForGene( Index geneID ) {
        return _genesToTranscripts[geneID];
    }

    inline std::string nameFromGeneID( Index geneID ) {
        return _geneNames[geneID];
    }

    inline Index gene( Index transcriptID ) {
        return _transcriptsToGenes[transcriptID];
    }

    inline std::string geneName( Index transcriptID ) {
        return _geneNames[_transcriptsToGenes[transcriptID]];
    }

    inline std::string geneName (const std::string& transcriptName,
                                 bool complain=true) {
        auto tid = findTranscriptID(transcriptName);
        if (tid != INVALID) {
            return geneName(tid);
        } else {
            std::cerr << "WARNING: couldn't find transcript named ["
                      << transcriptName << "]; returning transcript "
                      << " as it's own gene\n";
            return transcriptName;
        }
    }

    inline std::string transcriptName( Index transcriptID ) {
        return _transcriptNames[transcriptID];
    }
};

#endif // TRANSCRIPT_GENE_MAP_HPP
