#include "Material.h"
//
// created by rob on 01/12/2020. :)
//

Material::Material() = default;
Material::Material(std::string Name) : name(Name) {opticalDensity=0; usingTexture=0; usingBump=0; reflectivity = 0;}
