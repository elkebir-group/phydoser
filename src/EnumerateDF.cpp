/*
 * enumerateDF.cpp
 *
 *  Created on: 20-mar-2020
 *      Author: N. Aguse
 */

#include "Rcpp.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <set>
#include <algorithm>
#include <sstream>
#include <regex>
#include <string.h>
using namespace Rcpp;

typedef std::map<std::string,std::string> StringMap;
typedef std::map<std::string,int> StringIntMap;
typedef std::map<int,std::string> IntStringMap;
typedef std::vector<StringMap> StringMapVector;
typedef std::set<std::string> StringSet;
typedef std::vector<StringSet> StringSetVector;
typedef std::set<StringSet> StringSetSet;
typedef std::bitset<10000> Subset;
typedef std::vector<Subset> Family;
typedef std::vector<int> Feature;
typedef std::vector<Feature> FeatureFamily;
typedef std::vector<int> IntVector;
typedef std::set<int> IntSet;
typedef std::set<IntSet> IntSetSet;
typedef std::vector<IntSet> IntSetVector;
typedef std::vector<std::string> StringVector;

bool is_cover(Feature pi, Family& Fam, Subset Universe){
    Subset cover;
    for (int mut : pi){
        cover |= Fam[mut-1];
    }
    return cover == Universe;
}

void getDFF(Family& fam, Subset U, int m, FeatureFamily& Phi, FeatureFamily& allFeatures){
    //  FeatureFamily Phi;
    for (int i = 0; i < allFeatures.size(); ++i){
        Feature afeature = allFeatures[i];
        bool c = is_cover(afeature, fam, U);
        if (c){
            bool to_add = true;
            for (int j = 0; j < Phi.size(); ++j){
                Feature feature = Phi[j];
                if (std::includes(afeature.begin(),afeature.end(),feature.begin(),feature.end())){
                    to_add = false;
                    break;
                }
            }
            if (to_add){
                Phi.push_back(afeature);
            }
        }
        
    }
}

// [[Rcpp::export]]
Rcpp::List enumerateDF(Rcpp::List tree_mats){
    int nrTrees = tree_mats.length();
    Rcpp::NumericMatrix mat = tree_mats[0];
    int nrMuts = mat.cols();
    StringIntMap mutToId; // maps the column names to Id
    IntStringMap idToMut; // maps Id to column names
    int m = nrMuts;
    
    // Step 1. map mutations to ids and vice versa
    std::vector<std::string> mutNames = Rcpp::as<std::vector<std::string>>(Rcpp::colnames(tree_mats[0]));
    int idx = 0;
    for (auto mut : mutNames){
        mutToId[mut] = idx;
        idToMut[idx] = mut;
        idx++;
    }
    
    // Step 2. Convert each matrix the following way:
    //          - Each row becomes an intset (the mut id that corresponds to the column will be present in the set of m(row,col) = 1)
    //          - Each row is added to an intsetset (intsetset represents each matrix)
    std::vector<IntSetSet> treesFeatures(nrTrees);
    for (int i = 0; i < nrTrees; i++){ // for each tree
       Rcpp::checkUserInterrupt();
        std::vector<std::string> mutsList = Rcpp::as<std::vector<std::string>>(Rcpp::colnames(tree_mats[i])); // each tree may have different order of colnames
        for (int j = 0; j < nrMuts; j++){ // for each row (clone/feature)
            IntSet feature;
            for (int c = 0; c < nrMuts; c++){ // for each column (mutation)
                Rcpp::NumericMatrix mat = tree_mats[i];
                double elem = mat(j,c);
                if (elem == 1){
                    int mutid = mutToId[mutsList[c]];
                    feature.insert(mutid);
                }
            }
            treesFeatures[i].insert(feature);
        }
    }
    
    // Step 3. Enumerate all sets of size 1<= i <= m containing combinations of integers from 1 to m. Each combination of integers is called "Feature"
    FeatureFamily allFeatures;
    for (int i = 1; i <= m; ++i){
      Rcpp::checkUserInterrupt();
        Feature new_feature;
        for (int j = 1; j <= i; ++j){
            new_feature.push_back(j);
        }
        Feature add_feature2(new_feature);
        allFeatures.push_back(add_feature2);
        while(true){
            bool nothing_new = true;
            for (int idx = i-1; idx >= 0; --idx){
                if (new_feature[idx] < m && (idx == i-1 || new_feature[idx] < new_feature[idx+1]-1)){
                    new_feature[idx]+= 1;
                    int inc = new_feature[idx]+1;
                    // Reset all entries to the right of the incremented number
                    if(idx < i-1){
                        for(int idx2 = idx+1; idx2 < i; ++idx2){
                            new_feature[idx2] = inc++;
                        }
                    }
                    nothing_new = false;
                    break;
                }
            }
            if (nothing_new){
                break;
            }
            Feature add_feature2(new_feature);
            allFeatures.push_back(add_feature2);
        }
        
    }
    
    // Step 4. Convert features from step 2 to bitsets
    std::vector<IntSetVector> idxToFeatureMap2;
    std::vector<Family> FamilyVector2;
 
    for (int i = 0; i < nrTrees; ++i){
      Rcpp::checkUserInterrupt();
        Family fam; // Feature family of tree i (a collection of features)
        IntSetVector itfMap2(m);
        int idx = 0;
        for (auto it = treesFeatures[i].begin(); it != treesFeatures[i].end(); ++it){
            // For every mutation set (feature) of tree i, check if other trees have that feature
            Subset feat; // One feature of tree i
            
            itfMap2[idx++] = *it;
            for (int j = 0; j < nrTrees; ++j){
                if (i == j){
                    feat[j] = 1;
                }
                else if (treesFeatures[j].count(*it) == 0){
                    feat[j] = 1;
                }
            }
            fam.push_back(feat);
        }
        FamilyVector2.push_back(fam);
        idxToFeatureMap2.push_back(itfMap2);
    }
    
    // Step 5. Define the universe
    Subset Universe;
    for (int i = 0; i < nrTrees; ++i){
        Universe[i] = 1;
    }
    
    // Step 6. Enumerate distingushing features for each tree
    std::vector<FeatureFamily> PhiVector;
    for (auto fam : FamilyVector2){
        Rcpp::checkUserInterrupt();
        FeatureFamily Phi;
        getDFF(fam, Universe, m, Phi, allFeatures);
        PhiVector.push_back(Phi);
    }
    
    // Step 5. Convert distinguishing features to R List
    Rcpp::List ret = Rcpp::List::create();
    for (int i = 0; i < nrTrees; ++i){
        Rcpp::List treeDFF = Rcpp::List::create();
        for (auto Feature : PhiVector[i]){
            std::string featureStr = "";
            for (int mutIdx : Feature){
                for (auto mut : idxToFeatureMap2[i][mutIdx-1]){
                    featureStr += idToMut[mut];
                    featureStr += " ";
                }
                featureStr += ",";
            }
            featureStr.pop_back();
            featureStr.pop_back();
            treeDFF.push_back(featureStr);
        }
        ret.push_back(treeDFF);
    }
    return ret;
}