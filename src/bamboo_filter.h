#ifndef MY_BAMBOO_FILTER_H
#define MY_BAMBOO_FILTER_H

#include <string>
#include <vector>
#include <cstdint>
#include <utility> // For std::pair

class MyBambooFilter {
public:
    using Fp = std::uint16_t;
    using Slot = std::pair<Fp, std::uint64_t>;

    MyBambooFilter(std::size_t initial_num_buckets, std::size_t slots_per_bucket, float load_factor_threshold, std::size_t max_cuckoo_kicks);

    void insert(const std::string& key);

    bool contains(const std::string& key) const;

    std::size_t size() const;

    std::size_t capacity_buckets() const;

    float loadFactor() const;

    std::size_t memoryUsage() const;

private:
    std::vector<std::vector<Slot>> table_;

    std::size_t num_buckets_;
    std::size_t slots_per_bucket_;
    float max_load_factor_;
    std::size_t max_cuckoo_kicks_;
    std::size_t current_items_count_{0};

    void _attempt_insert_or_kick(std::uint64_t original_hash_of_item);

    void maybe_expand();

    void rebuild_table();

    // Hashing utility methods
    static std::uint64_t fnv1a_hash_str(const void* data, std::size_t len);
    static Fp fingerprint_from_hash_val(std::uint64_t h);
    static std::size_t index_from_hash_val(std::uint64_t h, std::size_t num_buckets_param);
    static std::size_t alt_index_from_fp_val(std::size_t primary_idx, Fp fp, std::size_t num_buckets_param);
};

#endif // MY_BAMBOO_FILTER_H