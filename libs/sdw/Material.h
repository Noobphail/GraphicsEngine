#pragma once
#include <glm/glm.hpp>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Utils.h"
#include "TextureMap.h"
#include "Colour.h"

struct Material {
	std::string name{};
    Colour ambientColour; //Ka
    Colour diffuseColour; //Kd
    Colour specularColour; //Ks
    float specularExponent{}; //Ns
	float reflectivity{}; // maybe use this
    float opticalDensity{}; //refractive ndex, Ni
    float transparency{}; //Tr
    int illuminationModel{}; //illum
    bool usingTexture{};
    bool usingBump{};
    TextureMap textureMap{}; //map_Kd

    Material();
    Material(std::string Name);
};