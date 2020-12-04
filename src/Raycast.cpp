//
// Created by ralph on 26/11/2020.
//

//casts a ray from cPos in rayDirection (with distance limit maxDistance), returning intersection information. Returns distance = -1 if no intersection found
RayTriangleIntersection getClosestIntersection(glm::vec3 cPos, glm::vec3 rayDirection, std::vector<ModelTriangle*> triangles, float maxDistance)
{
    float closestDepth = 999999;
    RayTriangleIntersection currentIntersection = RayTriangleIntersection();
    currentIntersection.distanceFromCamera = -1;
    for (int i = 0; i < triangles.size(); ++i)
    {
        glm::vec3 e0 = triangles[i]->vertices[1] - triangles[i]->vertices[0];
        glm::vec3 e1 = triangles[i]->vertices[2] - triangles[i]->vertices[0];
        glm::vec3 SPVector = cPos - triangles[i]->vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
        if ((possibleSolution[1] >= 0) && (possibleSolution[1] <= 1.0) && (possibleSolution[2] >= 0.0) && (possibleSolution[2] <= 1.0) && (possibleSolution[1] + possibleSolution[2]) <= 1.0)
        {
            if (possibleSolution[0] > 0.00001 && possibleSolution[0] < closestDepth && possibleSolution[0] <= maxDistance)
            {
                currentIntersection = RayTriangleIntersection(triangles[i]->vertices[0] + e0 * possibleSolution[1] + e1 * possibleSolution[2], possibleSolution[0], triangles[i], i);
                closestDepth = possibleSolution[0];
            }
        }
    }
    return currentIntersection;
}

//Returns an intersection between the triangle and a point if there is one, does not populate RayTriangleIntersection.triangleIndex
RayTriangleIntersection getIntersection(glm::vec3 cPos, glm::vec3 rayDirection, ModelTriangle* t, float maxDistance){
    glm::vec3 e0 = t->vertices[1] - t->vertices[0];
    glm::vec3 e1 = t->vertices[2] - t->vertices[0];
    glm::vec3 SPVector = cPos - t->vertices[0];
    glm::mat3 DEMatrix(-rayDirection, e0, e1);
    glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
    RayTriangleIntersection emptyIntersection = RayTriangleIntersection();
    emptyIntersection.distanceFromCamera = -1;
    if ((possibleSolution[1] >= 0) && (possibleSolution[1] <= 1.0) && (possibleSolution[2] >= 0.0) && (possibleSolution[2] <= 1.0) && (possibleSolution[1] + possibleSolution[2]) <= 1.0)
    {
        if (possibleSolution[0] > 0.00001 && possibleSolution[0] <= maxDistance)
        {
            return RayTriangleIntersection(t->vertices[0] + e0 * possibleSolution[1] + e1 * possibleSolution[2], possibleSolution[0], t, 0);
        }
    }
    return emptyIntersection;
}

bool checkForShadowIntersection(glm::vec3 cPos, glm::vec3 rayDirection, ModelTriangle* t, float lightDistance){
    glm::vec3 e0 = t->vertices[1] - t->vertices[0];
    glm::vec3 e1 = t->vertices[2] - t->vertices[0];
    glm::vec3 SPVector = cPos - t->vertices[0];
    glm::mat3 DEMatrix(-rayDirection, e0, e1);
    glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
    if ((possibleSolution[1] >= 0) && (possibleSolution[1] <= 1.0) && (possibleSolution[2] >= 0.0) && (possibleSolution[2] <= 1.0) && (possibleSolution[1] + possibleSolution[2]) <= 1.0)
    {
        if (possibleSolution[0] > 0.00001 && possibleSolution[0] <= lightDistance)
        {
            return true;
        }
    }
    return false;
}