#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

void mapIndices( std::vector<double>& primary, std::vector<double>& secondary, std::vector<int>& indices) {
    // Step 1: Duplicate the primary vector (secondary vector)
    //secondary = primary;
    
    // Step 2: Create a map between the original vector and its index
    // std::map<double, int> mapIndex;
    // for (int i = 0; i < primary.size(); ++i) {
    //     mapIndex[primary[i]] = i;
    // }
    
    std::vector<std::pair<double, int>> indexedValues;

    for( size_t i=0; i<primary.size(); ++i ){
      indexedValues.push_back(std::make_pair(primary[i], i));

    }

    // for (const auto& pair : indexedValues) {
    //   cout << pair.first << " ";
    // }
    // cout << endl;

    // Sort indexedValues based on the values in descending order
    std::sort(indexedValues.begin(), indexedValues.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first > rhs.first;
    });

    for (const auto& pair : indexedValues) {
      secondary.push_back(pair.first);
      indices.push_back(pair.second);
    }

    // for (const auto& pair : indexedValues) {
    //   cout << pair.first << " ";
    // }
    cout << endl;

    // Step 4: Create a third vector with primary vector index positions of each element in the secondary vector
    // indices.resize(secondary.size());
    // for (int i = 0; i < secondary.size(); ++i) {
    //     indices[i] = mapIndex[secondary[i]];
    // }
}

int sortmap() {
    // Create and populate the primary vector with unique doubles
    std::vector<double> primary = {3.0, 5.0, 4.0, 2.0, 6.0};
    
    // Declare the secondary vector and the vector for indices
    std::vector<double> secondary;
    std::vector<int> indices;
    
    // Call the function to map indices and sort
    mapIndices(primary, secondary, indices);
    
    // Print the results
    std::cout << "Primary Vector: ";
    for (const double& val : primary) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Secondary Vector: ";
    for (const double& val : secondary) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Indices Vector: ";
    for (const int& index : indices) {
        std::cout << index << " ";
    }
    std::cout << std::endl;
    
    return 0;
}
