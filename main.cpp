/********************************************************
 * main.cpp  –  Self-contained Bamboo Filter demonstration
 *
 * Implements:
 *   • 16-bit fingerprint Bamboo Filter with
 *       – smooth, segment-by-segment resizing (Wang et al.)
 *       – Cuckoo-style eviction for collisions (Fan et al.)
 *   • command-line interface, FASTA reader, random k-mer tester
 *
 * All standard-library; no external deps.
 * MIT-licensed – see LICENSE.
 *
 * Written 2025 by The Bamboo Filter Contributors
 ********************************************************/

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <random>

using namespace std;

static mt19937 rng(random_device{}());

inline uint32_t rnd32() { return rng(); }

/*----------------------------------------------------
  Zamjena za uniformno <0,n)
----------------------------------------------------*/
inline size_t rndIndex(size_t n) { return rnd32() % n; }

/*-------------------------------------------------------
  Simple struct for CLI parameters
 -------------------------------------------------------*/
struct Config
{
    string genomeFile = "";
    string outputFile = "results.txt";
    size_t capacity = 4096;
    size_t bucketSize = 4;
    float loadFactor = 0.8f;
    size_t maxIter = 2000;
    int k = 10;
    size_t numKmers = 1000;
    size_t segmentSize = 32;
};

/*-------------------------------------------------------
  CLI parsing (all flags have "--flag=value" form)
 -------------------------------------------------------*/
Config parseArgs(int argc, char **argv)
{
    Config c;
    auto getVal = [&](const char *arg, const char *name) -> const char *
    {
        std::size_t n = std::strlen(name);
        return std::strncmp(arg, name, n) == 0 ? arg + n : nullptr;
    };
    for (int i = 1; i < argc; ++i)
    {
        if (auto v = getVal(argv[i], "--genome="))
            c.genomeFile = v;
        else if (auto v = getVal(argv[i], "--out="))
            c.outputFile = v;
        else if (auto v = getVal(argv[i], "--capacity="))
            c.capacity = std::stoul(v);
        else if (auto v = getVal(argv[i], "--bucketSize="))
            c.bucketSize = std::stoul(v);
        else if (auto v = getVal(argv[i], "--loadFactor="))
            c.loadFactor = std::stof(v);
        else if (auto v = getVal(argv[i], "--maxIter="))
            c.maxIter = std::stoul(v);
        else if (auto v = getVal(argv[i], "--kmer="))
            c.k = std::stoi(v);
        else if (auto v = getVal(argv[i], "--numKmers="))
            c.numKmers = std::stoul(v);
        else if (auto v = getVal(argv[i], "--segmentSize="))
            c.segmentSize = std::stoul(v);
        else if (std::strcmp(argv[i], "--help") == 0 ||
                 std::strcmp(argv[i], "-h") == 0)
        {
            cout << "Run with --help to see usage in README.md\n";
            exit(0);
        }
        else
        {
            std::cerr << "Unknown arg: " << argv[i] << '\n';
            exit(1);
        }
    }
    return c;
}

/*-------------------------------------------------------
  FASTA / plain genome loader
 -------------------------------------------------------*/
string readGenome(const string &path)
{
    if (path.empty())
        return "";
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open genome file: " + path);
    string g, line;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.front() == '>')
            continue; // FASTA header
        g += line;
    }
    if (g.empty())
        throw std::runtime_error("Genome file empty/invalid");
    return g;
}

vector<string> sampleKmers(const string &g, int k, size_t n)
{
    vector<string> v;
    if (g.size() < (size_t)k)
        return v;
    v.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        size_t pos = rndIndex(g.size() - k + 1);
        v.emplace_back(g.substr(pos, k));
    }
    return v;
}

vector<string> randomStrings(int k, size_t n)
{
    static const char A[5] = {'A', 'C', 'G', 'T', 'N'};
    vector<string> v;
    v.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        string s(k, 'A');
        for (int j = 0; j < k; ++j)
            s[j] = A[rndIndex(5)];
        v.push_back(std::move(s));
    }
    return v;
}

/*-------------------------------------------------------
  16-bit fingerprint Bamboo Filter with
  – Cuckoo eviction
  – smooth segment migration
 -------------------------------------------------------*/
class BambooFilter
{
    using Fp = std::uint32_t;

public:
    BambooFilter(size_t cap, size_t bucketSz,
                 float loadFac, size_t maxEvict,
                 size_t segSz)
        : oldCap_(cap), B_(bucketSz), loadF_(loadFac),
          maxEvict_(maxEvict), segSz_(segSz)
    {
        old_.resize(oldCap_);
    }

    bool insert(const string &key)
    {
        if (contains(key))
            return true;
        maybeExpand();
        Fp fp = makeFp(key);
        size_t i1 = idxHash(key) % oldCap_;
        if (tryPut(old_, i1, fp) || cuckoo(old_, i1, fp, 0))
        {
            ++size_;
            return true;
        }
        // second candidate
        size_t i2 = altIndex(i1, fp, oldCap_);
        if (tryPut(old_, i2, fp) || cuckoo(old_, i2, fp, 0))
        {
            ++size_;
            return true;
        }
        return false; // rare
    }

    bool contains(const string &key)
    {
        Fp fp = makeFp(key);
        size_t i1 = idxHash(key) % oldCap_;
        size_t i2 = altIndex(i1, fp, oldCap_);
        if (has(old_, i1, fp) || has(old_, i2, fp))
            return true;
        if (expanding_)
        {
            size_t ni1 = i1 % newCap_;
            size_t ni2 = altIndex(ni1, fp, newCap_);
            return has(new_, ni1, fp) || has(new_, ni2, fp);
        }
        return false;
    }

    size_t size() const { return size_; }
    size_t capacity() const { return expanding_ ? oldCap_ + newCap_ : oldCap_; }

private:
    /* utilities */
    static size_t idxHash(const string &s)
    {
        static std::hash<string> H;
        return H(s);
    }
    static Fp makeFp(const string &s)
    {
        static std::hash<string> H;
        return static_cast<Fp>(H(s) & 0xffffu);
    }
    static size_t altIndex(size_t i, Fp fp, size_t cap)
    {
        return (i ^ (static_cast<size_t>(fp) * 0x5bd1e995)) % cap;
    }
    bool has(const vector<vector<Fp>> &arr,
             size_t idx, Fp fp) const
    {
        if (idx >= arr.size())
            return false;
        return std::find(arr[idx].begin(), arr[idx].end(), fp) != arr[idx].end();
    }
    bool tryPut(vector<vector<Fp>> &arr,
                size_t idx, Fp fp)
    {
        auto &b = arr[idx];
        if (b.size() < B_)
        {
            b.push_back(fp);
            return true;
        }
        return false;
    }

    bool cuckoo(vector<vector<Fp>> &arr,
                size_t idx, Fp fp, size_t depth)
    {
        if (depth >= maxEvict_)
            return false;
        auto &bucket = arr[idx];
        if (bucket.empty())
            return false;

        size_t victim = rndIndex(bucket.size()); // nasumična pozicija
        swap(fp, bucket[victim]);

        size_t alt = altIndex(idx, fp, arr.size());
        if (tryPut(arr, alt, fp))
            return true;
        return cuckoo(arr, alt, fp, depth + 1);
    }

    /* smooth expansion */
    void maybeExpand()
    {
        if (!expanding_)
        {
            float lf = float(size_) / float(oldCap_ * B_);
            if (lf <= loadF_)
                return;
            newCap_ = oldCap_ * 2;
            new_.resize(newCap_);
            expanding_ = true;
            migrateCursor_ = 0;
        }
        migrateSegment();
    }
    void migrateSegment()
    {
        size_t end = std::min(migrateCursor_ + segSz_, oldCap_);
        for (; migrateCursor_ < end; ++migrateCursor_)
        {
            for (Fp fp : old_[migrateCursor_])
            {
                size_t i1 = migrateCursor_ % newCap_;
                tryPut(new_, i1, fp) ||
                    cuckoo(new_, i1, fp, 0) ||
                    cuckoo(new_, altIndex(i1, fp, newCap_), fp, 0);
            }
            old_[migrateCursor_].clear();
        }
        if (migrateCursor_ == oldCap_)
        {
            old_ = std::move(new_);
            oldCap_ = newCap_;
            new_.clear();
            newCap_ = 0;
            expanding_ = false;
        }
    }

    /* data members */
    vector<vector<Fp>> old_, new_;
    size_t oldCap_, newCap_{0}, B_;
    float loadF_;
    size_t maxEvict_, segSz_;
    bool expanding_{false};
    size_t migrateCursor_{0};
    size_t size_{0};
};

/*-------------------------------------------------------
  Memory usage estimate (rough)
 -------------------------------------------------------*/
size_t memBytes(size_t buckets, size_t slots, bool fp16 = true)
{
    size_t perSlot = fp16 ? 2 : 1;
    return buckets * (24 + slots * perSlot); // 24B vector overhead ≈
}

/*-------------------------------------------------------
  MAIN
 -------------------------------------------------------*/
int main(int argc, char **argv)
{
    Config cfg = parseArgs(argc, argv);

    string genome;
    try
    {
        genome = readGenome(cfg.genomeFile);
    }
    catch (const std::exception &e)
    {
        if (!cfg.genomeFile.empty())
        {
            std::cerr << e.what() << '\n';
            return 1;
        }
    }

    auto positives = genome.empty()
                         ? randomStrings(cfg.k, cfg.numKmers)
                         : sampleKmers(genome, cfg.k, cfg.numKmers);
    auto negatives = randomStrings(cfg.k, cfg.numKmers);

    BambooFilter bf(cfg.capacity, cfg.bucketSize,
                    cfg.loadFactor, cfg.maxIter, cfg.segmentSize);

    auto tic = std::chrono::steady_clock::now();
    for (auto &s : positives)
        bf.insert(s);
    auto toc = std::chrono::steady_clock::now();
    double insertMs = std::chrono::duration<double, std::milli>(toc - tic).count();

    size_t tp = 0, fp = 0;
    tic = std::chrono::steady_clock::now();
    for (auto &s : positives)
        if (bf.contains(s))
            ++tp;
    for (auto &s : negatives)
        if (bf.contains(s))
            ++fp;
    toc = std::chrono::steady_clock::now();
    double queryMs = std::chrono::duration<double, std::milli>(toc - tic).count();

    std::ofstream out(cfg.outputFile);
    if (!out)
    {
        std::cerr << "Cannot open " << cfg.outputFile << '\n';
        return 1;
    }

    out << "# Bamboo Filter report\n";
    out << "Items inserted        : " << positives.size() << '\n';
    out << "Buckets total         : " << bf.capacity() << '\n';
    out << "Approx memory (bytes) : " << memBytes(bf.capacity(), cfg.bucketSize) << '\n';
    out << "Insertion time (ms)   : " << insertMs << '\n';
    out << "Query time (ms)       : " << queryMs << '\n';
    out << "True positives        : " << tp << "/" << positives.size() << '\n';
    out << "False positives       : " << fp << "/" << negatives.size() << '\n';
    out << "FP-rate               : " << (negatives.empty() ? 0.0 : double(fp) / negatives.size()) << '\n';

    cout << "Done. Results in " << cfg.outputFile << '\n';
    return 0;
}
