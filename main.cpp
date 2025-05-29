#include <iostream>
#include <vector>
#include <string>
#include "src/bamboo_filter.h"

int main() {
    std::size_t initial_buckets = 1024;
    std::size_t slots_per_bucket = 4;
    float load_factor = 0.95f;
    std::size_t max_kicks = 500;

    MyBambooFilter filter(initial_buckets, slots_per_bucket, load_factor, max_kicks);

    std::cout << "Kreiran MyBambooFilter s početnim kapacitetom od " << filter.capacity_buckets() << " bucketa." << std::endl;

    std::vector<std::string> items_to_insert = {"apple", "banana", "cherry", "date", "elderberry"};
    std::vector<std::string> items_to_check = {"apple", "banana", "fig", "grape"};

    std::cout << "\nUmetanje elemenata..." << std::endl;
    for (const auto& item : items_to_insert) {
        filter.insert(item);
        std::cout << "Umetnut: " << item << ", Trenutna veličina: " << filter.size() << ", Faktor popunjenosti: " << filter.loadFactor() << std::endl;
    }

    std::cout << "\nProvjera elemenata..." << std::endl;
    for (const auto& item : items_to_check) {
        if (filter.contains(item)) {
            std::cout << "'" << item << "' je vjerojatno u filteru." << std::endl;
        } else {
            std::cout << "'" << item << "' definitivno NIJE u filteru." << std::endl;
        }
    }
    
    std::cout << "\nFinalna veličina filtera: " << filter.size() << " elemenata." << std::endl;
    std::cout << "Finalni broj bucketa: " << filter.capacity_buckets() << std::endl;
    std::cout << "Procijenjena memorija: " << filter.memoryUsage() / 1024.0 << " KB" << std::endl;

    return 0;
}