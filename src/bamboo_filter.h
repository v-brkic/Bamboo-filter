#ifndef MY_BAMBOO_FILTER_H
#define MY_BAMBOO_FILTER_H

#include <string>
#include <vector>
#include <cstdint>
#include <utility> // For std::pair

/**
 * @file bamboo_filter.h
 * @brief Defines the MyBambooFilter class, a Cuckoo-style filter with expansion.
 *
 * This filter uses Cuckoo hashing for item placement and a table rebuilding
 * mechanism for expansion when the load factor exceeds a defined threshold.
 * It stores a 16-bit fingerprint along with the full 64-bit hash of the item
 * to ensure correct rebuilding and to aid in certain Cuckoo eviction scenarios.
 */
class MyBambooFilter {
public:
    /** @brief Type alias for the fingerprint (tag). */
    using Fp = std::uint16_t;
    /** @brief Type alias for a slot in the filter, storing fingerprint and full hash. */
    using Slot = std::pair<Fp, std::uint64_t>;

    /**
     * @brief Constructs a MyBambooFilter.
     * @param initial_num_buckets The initial number of buckets in the filter.
     * @param slots_per_bucket The number of slots (items) each bucket can hold before Cuckoo eviction or stashing.
     * @param load_factor_threshold The load factor at which the filter table rebuilds and expands.
     * @param max_cuckoo_kicks The maximum number of displacements allowed during a Cuckoo hashing attempt.
     */
    MyBambooFilter(std::size_t initial_num_buckets, std::size_t slots_per_bucket, float load_factor_threshold, std::size_t max_cuckoo_kicks);

    /**
     * @brief Inserts a key into the filter.
     * If the key is already likely present (based on a `contains` check),
     * the insertion might be skipped. This behavior is a trade-off for speed
     * when dealing with inputs containing many duplicates.
     * @param key The key to insert.
     */
    void insert(const std::string& key);

    /**
     * @brief Checks if a key is possibly in the filter.
     * This is a probabilistic check:
     * - If it returns false, the key is definitely not in the filter.
     * - If it returns true, the key may be in the filter (or it could be a false positive).
     * @param key The key to check.
     * @return True if the key might be in the filter, false otherwise.
     */
    bool contains(const std::string& key) const;

    /**
     * @brief Returns the number of items currently estimated to be in the filter.
     * This count reflects items successfully passed to the insertion logic.
     * @return The number of items.
     */
    std::size_t size() const;

    /**
     * @brief Returns the current number of buckets in the filter table.
     * This value increases when the filter expands.
     * @return The current number of buckets.
     */
    std::size_t capacity_buckets() const;

    /**
     * @brief Calculates the current load factor of the filter.
     * Load factor = (number of items) / (total number of slots).
     * @return The current load factor.
     */
    float loadFactor() const;

    /**
     * @brief Estimates the current memory usage of the filter in bytes.
     * This includes the main table structure and the stored slots.
     * @return Approximate memory usage in bytes.
     */
    std::size_t memoryUsage() const;

private:
    /** @brief The main table storing buckets, where each bucket is a vector of Slots. */
    std::vector<std::vector<Slot>> table_;

    /** @brief Current number of buckets in the filter. */
    std::size_t num_buckets_;
    /** @brief Number of slots each bucket can hold before Cuckoo/stash. */
    std::size_t slots_per_bucket_;
    /** @brief Load factor threshold that triggers table expansion. */
    float max_load_factor_;
    /** @brief Maximum number of kicks in a Cuckoo path before stashing. */
    std::size_t max_cuckoo_kicks_;
    /** @brief Number of items currently in the filter. */
    std::size_t current_items_count_{0};

    /**
     * @brief Internal method to perform the actual insertion logic (Cuckoo hashing, stashing).
     * This is called by both `insert()` and `rebuild_table()`.
     * @param original_hash_of_item The full 64-bit hash of the item to insert.
     */
    void _attempt_insert_or_kick(std::uint64_t original_hash_of_item);

    /**
     * @brief Checks if the filter needs to expand based on the current load factor
     * and triggers a rebuild if necessary.
     */
    void maybe_expand();

    /**
     * @brief Rebuilds the filter table, typically by doubling its capacity,
     * and re-inserts all existing items. This is a "stop-the-world" operation
     * that ensures all items are correctly placed after expansion using their full hashes.
     */
    void rebuild_table();

    // Hashing utility methods
    /**
     * @brief Computes a 64-bit hash (FNV-1a variant) for given data.
     * @param data Pointer to the data to hash.
     * @param len Length of the data in bytes.
     * @return 64-bit hash value.
     */
    static std::uint64_t fnv1a_hash_str(const void* data, std::size_t len);
    /**
     * @brief Extracts a 16-bit fingerprint from a 64-bit hash.
     * Ensures the fingerprint is non-zero.
     * @param h The 64-bit hash.
     * @return 16-bit fingerprint.
     */
    static Fp fingerprint_from_hash_val(std::uint64_t h);
    /**
     * @brief Calculates the primary bucket index from a 64-bit hash.
     * @param h The 64-bit hash.
     * @param num_buckets_param Current number of buckets in the table.
     * @return Bucket index.
     */
    static std::size_t index_from_hash_val(std::uint64_t h, std::size_t num_buckets_param);
    /**
     * @brief Calculates the alternate bucket index for Cuckoo hashing.
     * @param primary_idx The primary bucket index.
     * @param fp The fingerprint of the item.
     * @param num_buckets_param Current number of buckets in the table.
     * @return Alternate bucket index.
     */
    static std::size_t alt_index_from_fp_val(std::size_t primary_idx, Fp fp, std::size_t num_buckets_param);
};

#endif // MY_BAMBOO_FILTER_H