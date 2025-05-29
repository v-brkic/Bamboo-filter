#include "bamboo_filter.h"
#include <random>
#include <algorithm>
#include <stdexcept>    // For std::invalid_argument

// FNV-1a constants for 64-bit hash
constexpr std::uint64_t FNV_PRIME_64 = 0x100000001b3ULL;
constexpr std::uint64_t FNV_OFFSET_BASIS_64 = 0xcbf29ce484222325ULL;

//================================================================================
// Hashing Implementation
//================================================================================

std::uint64_t MyBambooFilter::fnv1a_hash_str(const void* data, std::size_t len) {
    auto p = static_cast<const unsigned char*>(data);
    std::uint64_t h = FNV_OFFSET_BASIS_64;
    for (size_t i = 0; i < len; ++i) {
        h = (h ^ p[i]) * FNV_PRIME_64;
    }
    return h;
}

MyBambooFilter::Fp MyBambooFilter::fingerprint_from_hash_val(std::uint64_t h) {
    Fp fp = h & 0xFFFF; // Use lower 16 bits for fingerprint
    return fp == 0 ? 1 : fp;
}

std::size_t MyBambooFilter::index_from_hash_val(std::uint64_t h, std::size_t num_buckets_param) {
    return (h >> 16) % num_buckets_param;
}

std::size_t MyBambooFilter::alt_index_from_fp_val(std::size_t primary_idx, Fp fp, std::size_t num_buckets_param) {
    std::uint64_t fp_intermediate_hash = static_cast<std::uint64_t>(fp) * 0x5bd1e995ULL; // Magic constant from MurmurHash
    return (primary_idx ^ fp_intermediate_hash) % num_buckets_param;
}

//================================================================================
// Constructor
//================================================================================

MyBambooFilter::MyBambooFilter(std::size_t initial_num_buckets_param, std::size_t slots_per_bucket_param,
                               float load_factor_threshold, std::size_t max_cuckoo_kicks_param)
  : num_buckets_(initial_num_buckets_param),
    slots_per_bucket_(slots_per_bucket_param),
    max_load_factor_(load_factor_threshold),
    max_cuckoo_kicks_(max_cuckoo_kicks_param),
    current_items_count_(0) {
    if (num_buckets_ == 0 || slots_per_bucket_ == 0) {
        throw std::invalid_argument("Number of buckets and slots per bucket must be greater than 0.");
    }
    table_.resize(num_buckets_); // Initialize the outer vector with num_buckets_ empty inner vectors
}

//================================================================================
// Public Methods: contains and insert
//================================================================================

bool MyBambooFilter::contains(const std::string& key) const {
    if (num_buckets_ == 0) return false; // Should not happen if constructor validation works

    const std::uint64_t h = fnv1a_hash_str(key.data(), key.length());
    const Fp fp_to_find = fingerprint_from_hash_val(h);
    const std::size_t i1 = index_from_hash_val(h, num_buckets_);

    if (i1 >= table_.size()) return false;

    for (const auto& slot : table_[i1]) {
        if (slot.first == fp_to_find) return true;
    }

    const std::size_t i2 = alt_index_from_fp_val(i1, fp_to_find, num_buckets_);
    if (i2 >= table_.size()) return false; // Defensive check

    for (const auto& slot : table_[i2]) {
        if (slot.first == fp_to_find) return true;
    }
    return false;
}

void MyBambooFilter::insert(const std::string& key) {
    if (contains(key)) {
        return;
    }

    maybe_expand();

    _attempt_insert_or_kick(fnv1a_hash_str(key.data(), key.length()));
    current_items_count_++;
}

//================================================================================
// Private Method: _attempt_insert_or_kick
//================================================================================

void MyBambooFilter::_attempt_insert_or_kick(std::uint64_t original_hash_of_item) {
    Slot slot_to_place = {fingerprint_from_hash_val(original_hash_of_item), original_hash_of_item};
    std::size_t i1 = index_from_hash_val(original_hash_of_item, num_buckets_);

    // Attempt to place in the primary bucket
    if (table_[i1].size() < slots_per_bucket_) {
        table_[i1].push_back(slot_to_place);
        return;
    }

    // Attempt to place in the alternate bucket
    std::size_t i2 = alt_index_from_fp_val(i1, slot_to_place.first, num_buckets_);
    if (table_[i2].size() < slots_per_bucket_) {
        table_[i2].push_back(slot_to_place);
        return;
    }

    static thread_local std::mt19937 rng(std::random_device{}());
    std::size_t current_bucket_idx = (rng() % 2 == 0) ? i1 : i2;

    for (std::size_t kick_count = 0; kick_count < max_cuckoo_kicks_; ++kick_count) {
        if (table_[current_bucket_idx].empty()) {
            table_[current_bucket_idx].push_back(slot_to_place);
            return;
        }

        std::uniform_int_distribution<std::size_t> dist(0, table_[current_bucket_idx].size() - 1);
        std::size_t victim_slot_in_bucket_offset = dist(rng);

        Slot temp_victim_slot = table_[current_bucket_idx][victim_slot_in_bucket_offset];
        table_[current_bucket_idx][victim_slot_in_bucket_offset] = slot_to_place;
        slot_to_place = temp_victim_slot;

        std::size_t victim_original_primary_idx = index_from_hash_val(slot_to_place.second, num_buckets_);
        if (current_bucket_idx == victim_original_primary_idx) {
            current_bucket_idx = alt_index_from_fp_val(victim_original_primary_idx, slot_to_place.first, num_buckets_);
        } else {
            current_bucket_idx = victim_original_primary_idx;
        }

        if (table_[current_bucket_idx].size() < slots_per_bucket_) {
            table_[current_bucket_idx].push_back(slot_to_place);
            return; // Successfully placed the kicked item
        }
    }

    table_[current_bucket_idx].push_back(slot_to_place);
}

//================================================================================
// Expansion Logic
//================================================================================

void MyBambooFilter::maybe_expand() {
    if (loadFactor() >= max_load_factor_) {
        rebuild_table();
    }
}

void MyBambooFilter::rebuild_table() {
    std::vector<std::uint64_t> all_original_hashes;
    all_original_hashes.reserve(current_items_count_);

    for (const auto& bucket : table_) {
        for (const auto& slot_item : bucket) {
            if (slot_item.first != 0) {
                all_original_hashes.push_back(slot_item.second);
            }
        }
    }

    num_buckets_ *= 2;
    table_.assign(num_buckets_, std::vector<Slot>()); // Create new, empty buckets

    current_items_count_ = 0;

    for (const auto& original_hash_to_reinsert : all_original_hashes) {
        _attempt_insert_or_kick(original_hash_to_reinsert);
        current_items_count_++; // Increment count for each successfully re-inserted item
    }
}

//================================================================================
// Utility Public Methods
//================================================================================

std::size_t MyBambooFilter::size() const {
    return current_items_count_;
}

std::size_t MyBambooFilter::capacity_buckets() const {
    return num_buckets_;
}

float MyBambooFilter::loadFactor() const {
    const std::size_t total_physical_slots = num_buckets_ * slots_per_bucket_;
    if (total_physical_slots == 0) return 0.0f;
    return static_cast<float>(current_items_count_) / total_physical_slots;
}

std::size_t MyBambooFilter::memoryUsage() const {
    std::size_t total_mem = sizeof(table_);
    for (const auto& bucket : table_) {
        total_mem += sizeof(std::vector<Slot>) + (bucket.capacity() * sizeof(Slot));
    }
    return total_mem;
}