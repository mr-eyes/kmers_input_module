#include "kmerDecoder.hpp"

void kmerDecoder::next_chunk(){
    this->kmers.clear();
    this->extractKmers();
}

flat_hash_map<std::string,std::vector<kmer_row>>* kmerDecoder::getKmers(){
    return &this->kmers;
}


void kmerDecoder::initialize_kSeq(){
    fp = gzopen(this->fileName.c_str(), "r");
    kseqObj = kseq_init(fp);
    KS_FULL_COMMENT = true; // Allow retrieving full header from kseq
}

bool kmerDecoder::end() const{
    return this->FILE_END;
}

std::string kmerDecoder::get_filename(){
    return this->fileName;
}

Hasher* kmerDecoder::initHasher(hashingModes HM, int kSize)
{
    switch (HM) {
        case integer_hasher:return new IntegerHasher(kSize);
        case mumur_hasher:return new MumurHasher(42);
        case TwoBits_hasher:return new TwoBitsHasher(kSize);
        case nonCanonicalInteger_Hasher:return new noncanonical_IntegerHasher(kSize);
        case protein_hasher:return new aaHasher_default(kSize);
        case proteinDayhoff_hasher:return new aaHasher_dayhoff(kSize);
        default: throw "unknown hashing mode";
    }
}

// -------------------------------------

inline flat_hash_map<std::tuple<enum readingModes, enum hashingModes>, bool> allowed_modes ({
    {{KMERS, mumur_hasher}, true},
    {{SKIPMERS, mumur_hasher}, true},
    {{MINIMIZERS, mumur_hasher}, true},
    {{PROTEIN, mumur_hasher}, false},

    {{KMERS, integer_hasher}, true},
    {{SKIPMERS, integer_hasher}, true},
    {{MINIMIZERS, integer_hasher}, true},
    {{PROTEIN, integer_hasher}, false},

    {{KMERS, TwoBits_hasher}, true},
    {{SKIPMERS, TwoBits_hasher}, true},
    {{MINIMIZERS, TwoBits_hasher}, true},
    {{PROTEIN, TwoBits_hasher}, false},

    {{KMERS, nonCanonicalInteger_Hasher}, true},
    {{SKIPMERS, nonCanonicalInteger_Hasher}, true},
    {{MINIMIZERS, nonCanonicalInteger_Hasher}, true},
    {{PROTEIN, nonCanonicalInteger_Hasher}, false},

    {{KMERS, protein_hasher}, false},
    {{SKIPMERS, protein_hasher}, false},
    {{MINIMIZERS, protein_hasher}, false},
    {{PROTEIN, protein_hasher}, true},

    {{KMERS, proteinDayhoff_hasher}, false},
    {{SKIPMERS, proteinDayhoff_hasher}, false},
    {{MINIMIZERS, proteinDayhoff_hasher}, false},
    {{PROTEIN, proteinDayhoff_hasher}, true},
});

 kmerDecoder * kmerDecoder::getInstance(readingModes RM, hashingModes HM, map<string, int> params) {
      if (!allowed_modes[{RM, HM}]) {
        cerr << "incompatible reading and hashing modes" << endl;
        exit(1);
        }

      // TODO validate parameters

      switch (RM){
        case(KMERS):
            return (new Kmers(params["kSize"], HM)); break;
        case SKIPMERS:
            return (new Skipmers(params["m"], params["n"], params["k"])); break;
        case MINIMIZERS:
            return new Minimizers(params["k"], params["w"]); break;
        case PROTEIN:
            return new aaKmers(params["kSize"], HM); break;
      }

    }

    kmerDecoder * kmerDecoder::getInstance(string fileName, int chunkSize, readingModes RM, hashingModes HM, map<string, int> params){
        auto * KD = getInstance(RM, HM, std::move(params));
        KD->fileName = fileName;
        KD->chunk_size = chunkSize;
        KD->initialize_kSeq();
        return KD;
  }


  std::map<std::string, int> kmerDecoder::string_to_params(std::string const& s){
      // adapted from https://stackoverflow.com/a/38814150/3371177
      std::map<std::string, int> m;

    std::string::size_type key_pos = 0;
    std::string::size_type key_end;
    std::string::size_type val_pos;
    std::string::size_type val_end;

    while((key_end = s.find(':', key_pos)) != std::string::npos)
    {
        if((val_pos = s.find_first_not_of(":", key_end)) == std::string::npos)
            break;

        val_end = s.find(',', val_pos);
        m.emplace(s.substr(key_pos, key_end - key_pos), stoi(s.substr(val_pos, val_end - val_pos)));

        key_pos = val_end;
        if(key_pos != std::string::npos)
            ++key_pos;
    }

    return m;
  }


bool valid_kmer(const string& kmer) {
    for (int i = 0; i < kmer.length(); i++) {
        char c = kmer[i];
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C') { return false; }
    }
    return true;
}