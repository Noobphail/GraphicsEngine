#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include "Material.h"
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <list>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <RayTriangleIntersection.h>
#include <ModelTriangle.h>
#include <map>
#include <glm/gtx/string_cast.hpp>
#include <algorithm>
#include "config.cpp"
#include "Helpers.cpp"
#include "Lines.cpp"
#include "Triangles.cpp"
#include "3DLines.cpp"
#include "3DTriangles.cpp"
#include "Raycast.cpp"
#include "ObjectParser.cpp"
# define M_PI           3.14159265358979323846

glm::vec4 cameraPos = glm::vec4(0, 0, 8, 1);
glm::mat4 cameraOrientation = glm::mat4x4(1, 0, 0, 0,
										  0, 1, 0, 0,
										  0, 0, 1, 0,
										  0, 0, 0, 1);
glm::mat3 globalOrientation = glm::mat3x3(1, 0, 0,
										  0, 1, 0,
										  0, 0, 1);
glm::vec3 sphereLightPos = glm::vec3(0, 1.5, 3);
std::vector<glm::vec3> lightPoints = {glm::vec3(0,-2.0,0),glm::vec3(0.1,-2.0,0),glm::vec3(-0.1,-2.0,0),glm::vec3(0,-2.0,0.1),glm::vec3(0,-2.0,-0.1)};
glm::vec3 lightPos = glm::vec3(0, -2.0, 0);
float lightIntensity = 300;
double minBrightness = 0.3;
float lightThreshold = 12288;

bool shiftHeld = false; //doesnt work if you release shift while processing is happening
int rendererSetting = 2;
const int numSettings = 3;
int lightSetting = 1;
//const int numLSettings = 6;

bool perspectiveCorrect = true;
bool redraw = true;

bool gouraudShading = false;
bool phongShading = false;
glm::vec4 lockedXBasis = glm::vec4(1, 0, 0, 0);
bool lockedXAxis = false;
bool softShadowsEnabled = true;

unsigned long renderCounter = 0;
std::vector<ModelTriangle *> trianglesChangedThisFrame;

std::vector<glm::vec3> createPlanarLight(float maxRadius, float ringSpacing, float multiplierPerRing, int numInFirstRing, glm::vec3 point, glm::vec3 normal ){

    std::vector<glm::vec3> lights;
    lights.push_back(point);
    if (ringSpacing == 0){
        return lights;
    }
    float currentRadius = ringSpacing;
    float currentAngle = 0;
    float currentIncrement = 2*M_PI / numInFirstRing;
    float currentNumLights = numInFirstRing;
    glm::vec3 temp = glm::vec3(0,1,0); //cross product doesn't work if vectors parallel
    if (normal.x == 0 && normal.z == 0){
        temp = glm::vec3(1,0,0);
    }
    glm::vec3 basis1 = glm::normalize(glm::cross(temp, normal));
    glm::vec3 basis2 = glm::normalize(glm::cross(basis1,normal));
    while (currentRadius < maxRadius){
        for (int i = 0; i < currentNumLights; ++i) {
            glm::vec3 newLight = point + currentRadius * basis1 * cos(currentAngle) +  currentRadius * basis2 * sin(currentAngle);
            lights.push_back(newLight);
            currentAngle += currentIncrement;
        }
        currentAngle = (currentAngle + currentIncrement/2) - (2*M_PI);
        currentRadius += ringSpacing;
        currentNumLights *= multiplierPerRing;
        currentIncrement = 2*M_PI / currentNumLights;
    }
    std::cout << "Light size: " << lights.size() << std::endl;
    return lights;
}

bool listContainsTriangle(std::vector<ModelTriangle *> list, ModelTriangle *ptr)
{
	for (int i = 0; i < list.size(); ++i)
	{
		if (list[i] == ptr)
		{
			return true;
		}
	}
	return false;
}

bool isShadowed(RayTriangleIntersection currentIntersection, std::vector<ModelTriangle *> triangles, glm::vec3 normal, glm::vec3 light)
{
	RayTriangleIntersection fullShadowCheck = getClosestIntersection(currentIntersection.intersectionPoint + glm::normalize(normal) * shadowBias, light - currentIntersection.intersectionPoint, triangles, 9999);
	if (fullShadowCheck.distanceFromCamera != -1)
	{
		if (glm::distance(currentIntersection.intersectionPoint, fullShadowCheck.intersectionPoint) < glm::distance(currentIntersection.intersectionPoint, light))
		{
			return true;
		}
	}
	return false;
}


bool getShadowed(glm::vec3 point, ModelTriangle *t, std::vector<ModelTriangle *> triangles, glm::vec3 normal)
{
	std::vector<ModelTriangle *> likelyShadowers = t->likelyShadowCasters;
	bool previousShadowerChanged = false;
	float lightDist = glm::distance(lightPos, point);
	for (int i = 0; i < likelyShadowers.size(); ++i)
	{
		if (checkForShadowIntersection(point + glm::normalize(normal) * shadowBias, lightPos - point, likelyShadowers[i], lightDist))
		{
			return true;
		}
		else
		{
			if (renderCounter == likelyShadowers[i]->frameLastChanged)
			{
				previousShadowerChanged = true;
			}
		}
	}
	//std::cout<<"2"<<std::endl;
	for (int i = 0; i < trianglesChangedThisFrame.size(); ++i)
	{
		std::cout << "triangles changed" << std::endl;
		if (checkForShadowIntersection(point + glm::normalize(normal) * shadowBias, lightPos - point, trianglesChangedThisFrame[i], lightDist))
		{
			if (!listContainsTriangle(likelyShadowers, trianglesChangedThisFrame[i]))
			{
				t->likelyShadowCasters.push_back(trianglesChangedThisFrame[i]);
			}
			return true;
		}
	}
	if (previousShadowerChanged)
	{
		std::cout << "full check" << std::endl;
		RayTriangleIntersection fullShadowCheck = getClosestIntersection(point + glm::normalize(normal) * shadowBias, lightPos - point, triangles, 9999);
		if (fullShadowCheck.distanceFromCamera != -1)
		{
			if (glm::distance(point, fullShadowCheck.intersectionPoint) < glm::distance(point, lightPos))
			{
				if (!listContainsTriangle(likelyShadowers, fullShadowCheck.intersectedTriangle))
				{
					t->likelyShadowCasters.push_back(fullShadowCheck.intersectedTriangle);
				}
				return true;
			}
		}
	}
	return false;
}

//returns the fraction of lights that hit the point 1 is fully lit 0 is dark
float getShadeStrength(std::vector<glm::vec3> pointLights, glm::vec3 normal, RayTriangleIntersection currentIntersection, std::vector<ModelTriangle*>triangles){
    float lights = pointLights.size();
    int hitLights = 0;
    for (int i = 0; i < pointLights.size(); ++i) {
        if (!isShadowed(currentIntersection,triangles,normal,pointLights[i])){
            hitLights++;
        }
    }
    return hitLights / lights;
}

double getBrightness(RayTriangleIntersection currentIntersection, std::vector<ModelTriangle *> triangles, glm::vec3 lightPosition, glm::vec3 normal)
{
	double proxScalar = 0;
	double incScalar = 0;
	double specScalar = 0;
	float lightness = 1;

	//RayTriangleIntersection shadowIntersection = getClosestIntersection(currentIntersection.intersectionPoint + glm::normalize(normal) * shadowBias, lightPosition - currentIntersection.intersectionPoint, triangles, 9999);
	bool shaded;
    if (softShadowsEnabled){
        lightness = getShadeStrength(lightPoints,normal,currentIntersection,triangles);
        shaded = (lightness == 0);
    }
    else{
        shaded = isShadowed(currentIntersection, triangles, normal, lightPosition);
    }

	if (shaded)
	{ //if not lit up
	}
	else
	{
		//proximity
		float distanceToLight = glm::length2(currentIntersection.intersectionPoint - lightPosition); //vscode shut up glm does have member length2
		if (distanceToLight == 0)
		{
			std::cout << "light inside object oh no" << std::endl;
		}
		proxScalar = lightIntensity / (4 * glm::pi<float>() * distanceToLight);
		//angle of incidence
		glm::vec3 lightVector = currentIntersection.intersectionPoint - lightPosition;
		incScalar = glm::dot(glm::normalize(lightVector), glm::normalize(normal));
		incScalar = clamp(incScalar, 0.0, 1.0);
		//specular
		glm::vec3 normDP = glm::normalize(normal) * (glm::dot(glm::normalize(lightVector), glm::normalize(normal)));
		glm::vec3 reflectionVec = glm::normalize(lightVector) - 2 * normDP;
		specScalar = glm::dot(glm::normalize(glm::vec3(cameraPos.x, cameraPos.y, cameraPos.z) - currentIntersection.intersectionPoint), glm::normalize(reflectionVec));
		specScalar = clamp(specScalar, 0, 1);
		specScalar = std::pow(specScalar, specExponent);
	}
	double brightness = glm::max(specScalar, proxScalar * incScalar);
	brightness = clamp(brightness, 0, 1);
	brightness = glm::max(brightness, minBrightness); //one line ambient lighting *dabs*
	return brightness * lightness + minBrightness * (1-lightness);
}

glm::vec3 barycentricCoefNormals(glm::vec3 p, ModelTriangle *t, std::array<glm::vec3, 3> vertexNormals)
{ //https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
	std::vector<double> baryCoefs = barycentricCoefs(p, t);
	return -(baryCoefs[0] * glm::normalize(vertexNormals[0]) + baryCoefs[1] * glm::normalize(vertexNormals[1]) + baryCoefs[2] * glm::normalize(vertexNormals[2]));
}

glm::vec3 barycentricCoefNormals(RayTriangleIntersection intersection)
{
	std::vector<double> baryCoefs = barycentricCoefs(intersection.intersectionPoint, intersection.intersectedTriangle);
	return -(baryCoefs[0] * glm::normalize(intersection.intersectedTriangle->vertexNorms[0]) + baryCoefs[1] * glm::normalize(intersection.intersectedTriangle->vertexNorms[1]) + baryCoefs[2] * glm::normalize(intersection.intersectedTriangle->vertexNorms[2]));
}

glm::vec3 getIntersectionNormal(RayTriangleIntersection intersection)
{
	if (intersection.intersectedTriangle->isPhongShaded)
	{
		return barycentricCoefNormals(intersection.intersectionPoint, intersection.intersectedTriangle, intersection.intersectedTriangle->vertexNorms);
	}
	else
	{
		return intersection.intersectedTriangle->normal;
	}
}

std::vector<std::vector<double>> getVertexBrightness(std::vector<ModelTriangle *> triangles)
{
	std::vector<std::vector<double>> vertexBrightness;
	RayTriangleIntersection currentIntersection; //= RayTriangleIntersection();//RayTriangleIntersection(triangles[0].vertices[0],0,triangles[0],0);
	for (int i = 0; i < triangles.size(); i++)
	{
		std::vector<double> currentBrightnesses;
		currentIntersection = RayTriangleIntersection(triangles[i]->vertices[0], 1, triangles[i], i);
		currentBrightnesses.push_back(getBrightness(currentIntersection, triangles, lightPos, -triangles[i]->vertexNorms[0]));
		currentIntersection = RayTriangleIntersection(triangles[i]->vertices[1], 1, triangles[i], i);
		currentBrightnesses.push_back(getBrightness(currentIntersection, triangles, lightPos, -triangles[i]->vertexNorms[1]));
		currentIntersection = RayTriangleIntersection(triangles[i]->vertices[2], 1, triangles[i], i);
		currentBrightnesses.push_back(getBrightness(currentIntersection, triangles, lightPos, -triangles[i]->vertexNorms[2]));
		vertexBrightness.push_back(currentBrightnesses);
	}
	return vertexBrightness;
}

double barycentricCoefBrightness(glm::vec3 p, ModelTriangle *t, std::vector<double> vertexBrightnesses)
{
	std::vector<double> baryCoefs = barycentricCoefs(p, t);

	//std::cout<< alpha << " " << vertexBrightnesses[0] << " " << beta << " " << vertexBrightnesses[1] << " " << gamma  <<" " << vertexBrightnesses[2] << " " << alpha * vertexBrightnesses[0] + beta * vertexBrightnesses[1] + gamma * vertexBrightnesses[2] << std::endl;
	return std::max(baryCoefs[0] * vertexBrightnesses[0] + baryCoefs[1] * vertexBrightnesses[1] + baryCoefs[2] * vertexBrightnesses[2], minBrightness);
}

Colour getRefractionColour(glm::vec3 incidenctRay, glm::vec3 point, glm::vec3 normal, std::vector<ModelTriangle *> triangles, float refractiveIndex, int rayDepth);

Colour getReflectionColour(glm::vec3 incidenceAng, glm::vec3 p, glm::vec3 normal, std::vector<ModelTriangle *> triangles, int currentDepth)
{
	Colour outCol;
	glm::vec3 N = glm::normalize(normal);
	glm::vec3 rI = glm::normalize(incidenceAng);
	glm::vec3 reflectedVector = glm::normalize(rI - 2 * N * (glm::dot(N, rI)));
	RayTriangleIntersection reflectIntersection = getClosestIntersection(p, reflectedVector, triangles, 30);
	if (reflectIntersection.distanceFromCamera != -1)
	{
		glm::vec3 nextNormal = getIntersectionNormal(reflectIntersection);
		if (reflectIntersection.intersectedTriangle->mat.reflectivity > minReflectivity && currentDepth < maxReflections)
		{
			return getReflectionColour(reflectIntersection.intersectionPoint - p, reflectIntersection.intersectionPoint, nextNormal, triangles, currentDepth + 1);
		}
		else if (reflectIntersection.intersectedTriangle->mat.opticalDensity > minOpticalDensity && currentDepth < maxReflections)
		{
			return getRefractionColour(glm::normalize(reflectedVector), reflectIntersection.intersectionPoint, reflectIntersection.intersectedTriangle->normal, triangles, reflectIntersection.intersectedTriangle->mat.opticalDensity, currentDepth + 1);
		}
		outCol = reflectIntersection.intersectedTriangle->colour;
		float brightness = getBrightness(reflectIntersection, triangles, lightPos, nextNormal);
		outCol.red = clamp(outCol.red * brightness, 0, 255);
		outCol.green = clamp(outCol.green * brightness, 0, 255);
		outCol.blue = clamp(outCol.blue * brightness, 0, 255);
	}
	else
	{
		outCol = Colour(0, 0, 0);
	}
	return outCol;
}

glm::vec3 calculateRefraction(glm::vec3 incidenctRay, glm::vec3 point, glm::vec3 normal, float refractiveIndex, int rayDepth)
{ //broken code from https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
	float cosi = clamp(glm::dot(incidenctRay, normal), -1, 1);
	float etai = 1;				 //incidence IOR
	float etat = refractiveIndex; //medium IOR
	bool outside = false;
	glm::vec3 tempNormal = normal;
	if (cosi < 0)
	{
		outside = true;
		cosi = -cosi;
	}
	else
	{
		std::swap(etai, etat);
		tempNormal = -normal;
	}
	float eta = etai / etat;
	float k = 1 - eta * eta * (1 - cosi * cosi);
	glm::vec3 refractedRay;
	if (k < 0)
	{
		k = fabsf(k);
	}
    refractedRay = eta * incidenctRay + (eta * cosi - glm::sqrt(k)) * tempNormal;
	return refractedRay;
}


Colour getRefractionColour(glm::vec3 incidenctRay, glm::vec3 point, glm::vec3 normal, std::vector<ModelTriangle *> triangles, float opticalDensity, int rayDepth)
{

	Colour outCol;
	incidenctRay = glm::normalize(incidenctRay);
	Colour reflectionCol, refractionCol;

	glm::vec3 refractedRay = calculateRefraction(incidenctRay, point, normal, opticalDensity, rayDepth);

	RayTriangleIntersection nextIntersection = getClosestIntersection(point, refractedRay, triangles, 999);
	if (nextIntersection.distanceFromCamera != -1)
	{
		glm::vec3 intersectedNormal = getIntersectionNormal(nextIntersection);

		if (nextIntersection.intersectedTriangle->mat.reflectivity > minReflectivity && rayDepth < maxReflections)
		{
			return getReflectionColour(refractedRay, nextIntersection.intersectionPoint, intersectedNormal, triangles, rayDepth + 1);
		}
		else if (nextIntersection.intersectedTriangle->mat.opticalDensity > minOpticalDensity && rayDepth < maxReflections)
		{
			return getRefractionColour(refractedRay, nextIntersection.intersectionPoint, nextIntersection.intersectedTriangle->normal, triangles, nextIntersection.intersectedTriangle->mat.opticalDensity, rayDepth + 1);
		}
		else
		{
			glm::vec3 intersectedNormal = getIntersectionNormal(nextIntersection);
			outCol = nextIntersection.intersectedTriangle->colour;
			float brightness = getBrightness(nextIntersection, triangles, lightPos, intersectedNormal);
			outCol.red = clamp(outCol.red * brightness, 0, 255);
			outCol.green = clamp(outCol.green * brightness, 0, 255);
			outCol.blue = clamp(outCol.blue * brightness, 0, 255);
		}
	}
	else
	{
		outCol = Colour(0, 0, 0);
	}
	return outCol;
}

Colour getTexturePointColour(ModelTriangle *t, glm::vec3 p)
{
	std::vector<double> baryCoefs = barycentricCoefs(p, t);
	float texturePointX = baryCoefs[0] * t->texturePoints[0].x + baryCoefs[1] * t->texturePoints[1].x + baryCoefs[2] * t->texturePoints[2].x;
	float texturePointY = baryCoefs[0] * t->texturePoints[0].y + baryCoefs[1] * t->texturePoints[1].y + baryCoefs[2] * t->texturePoints[2].y;

	return int32ToCol(t->mat.textureMap.pixels[texturePointToIndex(TexturePoint(texturePointX, texturePointY), t->mat.textureMap)]);
}

void raycastDraw(DrawingWindow &window, std::vector<ModelTriangle *> triangles, float scaleFactor)
{
	/*trianglesChangedThisFrame.clear();
	for (int i = 0; i < triangles.size(); ++i)
	{
		if (triangles[i]->frameLastChanged == renderCounter)
		{
			trianglesChangedThisFrame.push_back(triangles[i]);
		}
	}*/
	glm::mat3 cameraRot = glm::mat3(cameraOrientation);
	glm::vec3 right = -cameraRot[0];
	glm::vec3 up = cameraRot[1];
	glm::vec3 forward = -cameraRot[2];
	glm::vec3 camPos = glm::vec3(cameraPos);
	glm::vec3 imagePlaneCentre = camPos + forward * di;
	glm::vec3 topLeftImagePlane = imagePlaneCentre - right * ((WIDTH / scaleFactor) / 2) + up * ((HEIGHT / scaleFactor) / 2);
	glm::vec3 rightStep = right * 1 / scaleFactor;
	glm::vec3 upStep = up * 1 / scaleFactor;
    std::vector<std::vector<double>> vertexBrightness;
    if (gouraudShading) {
        vertexBrightness = getVertexBrightness(triangles);
    }
	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{

			RayTriangleIntersection currentIntersection = getClosestIntersection(camPos, (topLeftImagePlane + rightStep * x - upStep * y) - camPos, triangles, 99999);
			if (currentIntersection.distanceFromCamera != -1)
			{
				double brightness;
				Colour outCol;
				glm::vec3 intersectedNormal = getIntersectionNormal(currentIntersection);
				if (currentIntersection.intersectedTriangle->mat.reflectivity > minReflectivity)
				{
					outCol = getReflectionColour(currentIntersection.intersectionPoint - camPos, currentIntersection.intersectionPoint, intersectedNormal, triangles, 0);
				}
				else if (currentIntersection.intersectedTriangle->mat.opticalDensity > minOpticalDensity)
				{
					outCol = getRefractionColour(currentIntersection.intersectionPoint - camPos, currentIntersection.intersectionPoint, currentIntersection.intersectedTriangle->normal, triangles, currentIntersection.intersectedTriangle->mat.opticalDensity, 0);
				}
				else
				{
					if (currentIntersection.intersectedTriangle->mat.usingTexture)
					{
						outCol = getTexturePointColour(currentIntersection.intersectedTriangle, currentIntersection.intersectionPoint);
					}
					else
					{
						outCol = Colour(currentIntersection.intersectedTriangle->colour);
					}
					if (gouraudShading)
					{
						/*RayTriangleIntersection shadowIntersection = getClosestIntersection(currentIntersection.intersectionPoint, lightPos - currentIntersection.intersectionPoint, triangles, 9999);
                        if (shadowIntersection.distanceFromCamera == -1 || glm::distance(currentIntersection.intersectionPoint, shadowIntersection.intersectionPoint) > glm::distance(currentIntersection.intersectionPoint, lightPos))
                        {//if lit up
                            brightness = barycentricCoefBrightness(currentIntersection.intersectionPoint,currentIntersection.intersectedTriangle,vertexBrightness[currentIntersection.triangleIndex]);
                        }
                        else{
                            brightness = minBrightness;
                        }*/
						brightness = barycentricCoefBrightness(currentIntersection.intersectionPoint, currentIntersection.intersectedTriangle, vertexBrightness[currentIntersection.triangleIndex]);
						brightness = std::max(brightness, minBrightness);
					}
					else
					{
						glm::vec3 currentNormal = getIntersectionNormal(currentIntersection);
						brightness = getBrightness(currentIntersection, triangles, lightPos, currentNormal);
					}
					outCol.red = clamp(outCol.red * brightness, 0, 255);
					outCol.green = clamp(outCol.green * brightness, 0, 255);
					outCol.blue = clamp(outCol.blue * brightness, 0, 255);
				}

				window.setPixelColour(x, y, colToInt32(outCol));
			}
			else
			{
				window.setPixelColour(x, y, colToInt32(Colour(0, 0, 0)));
			}
		}
	}
}

//returns vector corresponding to a 3D image projected to the 2D image plane, scaled by a given factor
std::vector<CanvasTriangle> projectToImagePlane(DrawingWindow &window, std::vector<ModelTriangle *> triangles, float scaleFactor)
{
	glm::mat3 camOrientationMatrix = glm::mat3(cameraOrientation);
	std::vector<CanvasTriangle> ret;
	float zBuffer[WIDTH][HEIGHT];
	for (int i = 0; i < WIDTH; i++)
	{
		for (int j = 0; j < HEIGHT; j++)
		{
			zBuffer[i][j] = 0.0;
		}
	}
	glm::vec3 camZ = glm::normalize(glm::vec3(camOrientationMatrix[2]));
	glm::vec3 camPos = glm::vec3(cameraPos);
	for (int i = 0; i < triangles.size(); i++)
	{
		CanvasPoint vertices[3] = {};
		Material mat = triangles[i]->mat;
		for (int j = 0; j < 3; j++)
		{
			glm::vec3 currentVertex = glm::vec3(triangles[i]->vertices[j].x - cameraPos.x, triangles[i]->vertices[j].y - cameraPos.y, triangles[i]->vertices[j].z - cameraPos.z);
			glm::vec3 res = currentVertex * camOrientationMatrix;
			float x = res.x;
			float y = res.y;
			float z = res.z;
			float ui = di * (x / z) + window.width / 2;
			float vi = di * (y / z) + window.height / 2;
			ui = scaleFactor * (ui - window.width / 2) + window.width / 2;
			vi = scaleFactor * (vi - window.height / 2) + window.height / 2;
			vertices[j] = CanvasPoint(round(ui), round(vi));
			vertices[j].texturePoint = triangles[i]->texturePoints[j];
			vertices[j].depth = glm::dot(camPos - currentVertex, camZ) - di;
			vertices[j].inverseDepth = 1 / vertices[j].depth;
		}
		CanvasTriangle t = CanvasTriangle(vertices[0], vertices[1], vertices[2]);
		if (rendererSetting == 0)
		{
			drawStrokedTriangle(window, t, Colour(255, 255, 255));
		}
		else
		{
			if (triangles[i]->mat.usingTexture)
			{
				if (perspectiveCorrect)
				{
					texturePerspetiveTriangle3D(window, t, triangles[i]->mat.textureMap, zBuffer);
				}
				else
				{
					textureTriangle3D(window, t, triangles[i]->mat.textureMap, zBuffer);
				}
			}
			else
			{
				fillTriangle3D(window, t, triangles[i]->colour, zBuffer);
			}
		}
		ret.push_back(t);
	}
	return ret;
}

//now used to draw the model, perhaps make a global struct for currently loaded models and textures?
void draw(DrawingWindow &window, std::vector<ModelTriangle *> triangles)
{
	window.clearPixels();
	if (rendererSetting == 0 || rendererSetting == 1)
	{
		projectToImagePlane(window, triangles, 100);
	}
	else if (rendererSetting == 2)
	{
		raycastDraw(window, triangles, 100);
	}
	redraw = false;
	renderCounter++;
}

void update(DrawingWindow &window)
{
	// Function for performing animation (shifting artifacts or moving the camera)
}

void lookAt(glm::vec3 vector)
{
	glm::mat3 camMatrix = glm::mat3(cameraOrientation);
	camMatrix[2] = glm::normalize(glm::vec3(cameraPos.x, cameraPos.y, cameraPos.z) - vector);
	camMatrix[0] = glm::normalize(glm::cross(glm::vec3(0, -1, 0), camMatrix[2]));
	camMatrix[1] = glm::normalize(glm::cross(camMatrix[2], camMatrix[0]));
	cameraOrientation = glm::mat4(camMatrix);
	redraw = true;
	std::cout << "New position: " << glm::to_string(cameraPos) << std::endl;
	std::cout << "New oreintation: " << glm::to_string(cameraOrientation) << std::endl;
}

glm::mat4 getRelativeAxisRotationMatrix(float angle, glm::vec4 basis)
{
	float uX = basis.x;
	float uY = basis.y;
	float uZ = basis.z;
	float cosTheta = cos(angle);
	float sinTheta = sin(angle);
	return glm::mat4(
		cosTheta + uX * uX * (1 - cosTheta), uX * uY * (1 - cosTheta) + uZ * sinTheta, uX * uZ * (1 - cosTheta) - uY * sinTheta, 0.0, // first column (not row!)
		uX * uY * (1 - cosTheta) - uZ * sinTheta, cosTheta + uY * uY * (1 - cosTheta), uZ * uY * (1 - cosTheta) + uX * sinTheta, 0.0, // second column
		uX * uZ * (1 - cosTheta) + uY * sinTheta, uY * uZ * (1 - cosTheta) - uX * sinTheta, cosTheta + uZ * uZ * (1 - cosTheta), 0.0, // third column
		0.0, 0.0, 0.0, 1);
}

//ORDER DEPENDENT
void orientCameraX(float angle)
{
	//X-rotation
	cameraOrientation *= getRelativeAxisRotationMatrix(angle, cameraOrientation[0]);
	redraw = true;
}

void orientCameraY(float angle)
{
	//Y-rotation
	cameraOrientation *= glm::mat4(cos(angle), 0, -1 * sin(angle), 0,
								   0, 1, 0, 0,
								   sin(angle), 0, cos(angle), 0,
								   0, 0, 0, 1);

	redraw = true;
}

void moveCamera(glm::vec3 translation)
{
	glm::mat4 translationMatrix = glm::mat4(1, 0, 0, translation[0],
											0, 1, 0, translation[1],
											0, 0, 1, translation[2],
											0, 0, 0, 1);
	cameraPos = cameraPos * translationMatrix;
	redraw = true;
	std::cout << "Cam pos: " << cameraPos[0] <<","<< cameraPos[1] <<","<< cameraPos[2] << std::endl;
}

//rotates camera theta degrees around the x-axis
void rotateCameraAroundX(float thetaRadians)
{
	if (!lockedXAxis)
	{
		lockedXBasis = cameraOrientation[0];
		lockedXAxis = true;
	}
	glm::mat4 rotMat = getRelativeAxisRotationMatrix(thetaRadians, lockedXBasis);
	cameraPos = cameraPos * rotMat;
	redraw = true;
	std::cout << "campos: " << cameraPos[0] << cameraPos[1] << cameraPos[2] << std::endl;
}

//rotates camera theta degrees around the x-axis
void rotateCameraAroundY(float thetaRadians)
{
	glm::mat4 rotMat = glm::mat4(
		cos(thetaRadians), 0.0, -sin(thetaRadians), 0.0, // first column (not row!)
		0.0, 1.0, 0.0, 0.0,								 // second column
		sin(thetaRadians), 0.0, cos(thetaRadians), 0.0,	 // third column
		0.0, 0.0, 0.0, 1);
	cameraPos = cameraPos * rotMat;
	redraw = true;
	std::cout << "qq " << cameraPos[0] << ", " << cameraPos[1] << ", " << cameraPos[2] << std::endl;
}

void handleEvent(SDL_Event event, DrawingWindow &window)
{
	if (event.type == SDL_KEYDOWN)
	{
		if (event.key.keysym.sym == SDLK_LEFT)
			moveCamera(glm::vec3(cameraOrientation[0]));
		else if (event.key.keysym.sym == SDLK_RIGHT)
			moveCamera(-glm::vec3(cameraOrientation[0]));
		else if (event.key.keysym.sym == SDLK_UP) //into screen
			moveCamera(-glm::vec3(cameraOrientation[2]));
		else if (event.key.keysym.sym == SDLK_DOWN) //out of screen
			moveCamera(glm::vec3(cameraOrientation[2]));
		else if (event.key.keysym.sym == SDLK_w) //move up
			moveCamera(glm::vec3(cameraOrientation[1]));
		else if (event.key.keysym.sym == SDLK_s) //move down
			moveCamera(-glm::vec3(cameraOrientation[1]));
		else if (event.key.keysym.sym == SDLK_i)
		{ //rotate camera down(kinda)
			rotateCameraAroundX(0.261799);
			lookAt(glm::vec3(0, 0, 0));
		}
		else if (event.key.keysym.sym == SDLK_k)
		{ //rotate camera down(kinda)
			rotateCameraAroundX(-0.261799);
			lookAt(glm::vec3(0, 0, 0));
		}
        else if (event.key.keysym.sym == SDLK_q)
        { //rotate camera down(kinda)
            softShadowsEnabled = !softShadowsEnabled;
            redraw = true;
        }
		else if (event.key.keysym.sym == SDLK_j)
		{ //rotate camera left(kinda)
			rotateCameraAroundY(0.261799);
			lockedXAxis = false;
			lookAt(glm::vec3(0, 0, 0));
		}
		else if (event.key.keysym.sym == SDLK_l)
		{ //rotate camera right(kinda)
			rotateCameraAroundY(-0.261799);
			lockedXAxis = false;
			lookAt(glm::vec3(0, 0, 0));
		}
		else if (event.key.keysym.sym == SDLK_t) //pitch camera up
			orientCameraX(-0.261799);
		else if (event.key.keysym.sym == SDLK_g) //pitch camera down
			orientCameraX(0.261799);
		else if (event.key.keysym.sym == SDLK_f) //yaw camera left
			orientCameraY(-0.261799);
		else if (event.key.keysym.sym == SDLK_h) //yaw camera right
			orientCameraY(0.261799);
		else if (event.key.keysym.sym == SDLK_u)
		{
			CanvasPoint v0 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			Colour col = Colour(rand() % 256, rand() % 256, rand() % 256);
			drawStrokedTriangle(window, v0, v1, v2, col);
		}
		else if (event.key.keysym.sym == SDLK_COLON) //t
		{
			CanvasPoint v0 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			Colour col = Colour(rand() % 256, rand() % 256, rand() % 256);
			CanvasTriangle t = CanvasTriangle(v0, v1, v2);
			fillTriangle(window, t, col);
		}
		else if (event.key.keysym.sym == SDLK_AT) //g
		{
			CanvasPoint v0 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			Colour col = Colour(rand() % 256, rand() % 256, rand() % 256);
			drawStrokedTriangle(window, v0, v1, v2, Colour(255, 255, 255));
			CanvasTriangle t = CanvasTriangle(v0, v1, v2);
			fillTriangle(window, t, col);
		}
		else if (event.key.keysym.sym == SDLK_HASH) //t
		{
			TexturePoint t1 = TexturePoint(195, 5);
			TexturePoint t2 = TexturePoint(395, 380);
			TexturePoint t3 = TexturePoint(65, 330);
			CanvasPoint p1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			p1.texturePoint = t1;
			CanvasPoint p2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			p2.texturePoint = t2;
			CanvasPoint p3 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
			p3.texturePoint = t3;
			CanvasTriangle t = CanvasTriangle(p1, p2, p3);
			textureTriangle(window, t);
		}
		else if (event.key.keysym.sym == SDLK_1)
		{
			rendererSetting = 0;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_2)
		{
			rendererSetting = 1;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_3)
		{
			rendererSetting = 2;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_LEFTBRACKET)
		{
			rendererSetting = (rendererSetting - 1) % numSettings;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_e)
		{
            std::cout << "perspective : " << perspectiveCorrect << std::endl;
			perspectiveCorrect = !perspectiveCorrect;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_KP_2)
		{
			lightPos.y += 0.1;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_KP_8)
		{
			lightPos.y -= 0.1;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_KP_4)
		{
			lightPos.x -= 0.1;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_KP_6)
		{
			lightPos.x += 0.1;
			redraw = true;
		}
		else if (event.key.keysym.sym == SDLK_LSHIFT)
		{
			shiftHeld = true;
			std::cout << "shift down" << std::endl;
		}
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN)
	{
		window.savePPM("output.ppm");
	}

	else if (event.type == SDL_KEYUP)
	{
		if (event.key.keysym.sym == SDLK_LSHIFT)
		{
			shiftHeld = false;
			std::cout << "shift up" << std::endl;
		}
	}
}

int main(int argc, char *argv[])
{
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<glm::vec3> result;
	Colour black = Colour(0, 0, 0);
	glm::vec3 from(1, 4, 9.2);
	glm::vec3 to(4, 1, 9.8);
	Colour white = Colour(255, 255, 255);
	std::map<std::string, Material> matMap;
	std::map<std::string, TextureMap> textMap;
	std::vector<ModelTriangle> triangles;
	std::vector<ModelTriangle *> renderedTriangles;

	if (loadSphere)
	{
		loadMTL("src/cornell-box.mtl", &matMap);		  //, &textMap);
		triangles = loadOBJ("src/sphere.obj", 1, matMap); //, textMap);
		sphereLightPos = glm::vec3(4,1.5,8);
		lightPos = sphereLightPos;
		lightIntensity = 600;
		rendererSetting = 2;
		//gouraudShading = true;
		//phongShading = false;
	}
	else if (loadTextured)
	{
		loadMTL("src/textured-cornell-box.mtl", &matMap);				//, &textMap);
		triangles = loadOBJ("src/textured-cornell-box.obj", 1, matMap); //, textMap);
		rendererSetting = 2;
	}
	else
	{
		loadMTL("src/cornell-box.mtl", &matMap);			   //, &textMap);
		triangles = loadOBJ("src/cornell-box.obj", 1, matMap); //, textMap);
		rendererSetting = 2;
        lightPoints = createPlanarLight(0.11,0.1,2,7,lightPos,glm::vec3(0,1,0));
	}
	for (int i = 0; i < triangles.size(); ++i)
	{
		renderedTriangles.push_back(&triangles[i]);
	}
	redraw = true;
	lookAt(glm::vec3(0, 0, 0));
	Colour test = Colour(255, 123, 222);

	std::cout << "COLOUR " << int32ToCol(colToInt32(test)) << std::endl;
	while (true)
	{
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event))
			handleEvent(event, window);
		//update(window);
		if (redraw)
			draw(window, renderedTriangles);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}