#pragma once
#include <glm/glm.hpp>
#include <string>
#include <array>
#include "Material.h"
#include <vector>
#include "Colour.h"
#include "TexturePoint.h"

struct ModelTriangle {
	std::array<glm::vec3, 3> vertices{};
	std::array<glm::vec3, 3> vertexNorms{};
	std::array<TexturePoint, 3> texturePoints{};
	Colour colour{};
	glm::vec3 normal{};
	Material mat;
	//float reflectivity;
	bool isPhongShaded;
	std::vector<ModelTriangle*> likelyShadowCasters;
	unsigned long frameLastChanged;

	ModelTriangle();
	ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
