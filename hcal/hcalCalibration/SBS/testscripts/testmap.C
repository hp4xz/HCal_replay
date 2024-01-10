#include <iostream>
#include <map>
#include <vector>

int testmap() {
    // Original set of indices
    std::vector<int> originalIndices = {31, 45, 49, 0, 2};
    std::vector<int> newids = {0,0,1,4,7};

    // Perform some operations to reorder the original indices
    // For example, let's reverse the order
    //std::reverse(originalIndices.begin(), originalIndices.end());

    // Create a map to relate original indices to reordered indices
    std::map<int, int> indexMap;

    // Populate the map with original indices as keys and reordered indices as values
    for (int i = 0; i < originalIndices.size(); ++i) {
        indexMap[originalIndices[i]] = newids[i];
    }

    // Now you can use the map to look up reordered indices for original indices
    // For example, let's look up the reordered index for the original index 4
    int originalIndexToLookup = 0;
    int reorderedIndex = indexMap[originalIndexToLookup];

    std::cout << "Original Index: " << originalIndexToLookup << " Reordered Index: " << reorderedIndex << std::endl;

    // You can iterate over the map to access all original-to-reordered index mappings
    std::cout << "Original to Reordered Index Map:" << std::endl;
    for (const auto& pair : indexMap) {
        std::cout << "Original Index: " << pair.first << " newid Index: " << pair.second << std::endl;
    }
    

    return 0;
}
