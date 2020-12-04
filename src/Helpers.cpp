//
// Created by ralph on 26/11/2020. but why?
//

template <typename T>
//get sign of a T
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

void print(std::string str){
    std::cout<<str<<std::endl;
}

Colour int32ToCol(int col)
{
    return Colour((col>>16)&255,(col>>8)&255,col&255);
}

//convert from Colour type to packed uint32_t
uint32_t colToInt32(Colour col)
{
    return (255 << 24) + ((int)(col.red) << 16) + (int(col.green) << 8) + int(col.blue);
}


float clamp(float val, float low, float high)
{
    return (val < low) ? low : (high < val) ? high : val;
}


std::string printVector(glm::vec3 vec) {
    return "x: " + std::to_string(vec.x) + " y: " + std::to_string(vec.y) + " z: " + std::to_string(vec.z);
}

//return vector of numberOfValues interpolated values between from and to
std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues)
{
    std::vector<float> ret;
    ret.assign(0, numberOfValues);
    if (numberOfValues == 1)
    {
        ret.push_back(from);
    }
    else
    {
        if (from == to)
        {
            for (int i = 0; i < numberOfValues; i++)
            {
                ret.push_back(from);
            }
        }
        else
        {
            for (int i = 0; i < numberOfValues; i++)
            {
                ret.push_back(from + (i) * (to - from) / (numberOfValues - 1));
            }
        }
    }
    return ret;
}

//curse you vscode >:(
void sortPointsByY(CanvasPoint *points, int length)
{
    for (int j = 0; j < length - 1; j++)
    {
        int sortedNum = 0;
        for (int i = 0; i < length - 1; i++)
        {
            if (points[i].y > points[i + 1].y)
            {
                std::swap(points[i], points[i + 1]);
            }
            else
            {
                sortedNum++;
            }
        }
        if (sortedNum == length - 1)
        {
            return;
        }
    }
}

//converts texture point to index in texturemap
int texturePointToIndex(TexturePoint p, TextureMap m)
{
    return round(p.y) * m.width + round(p.x);
}


std::vector<double> barycentricCoefs(glm::vec3 p, ModelTriangle *t){
    glm::vec3 p1 = t->vertices[0];
    glm::vec3 p2 = t->vertices[1];
    glm::vec3 p3 = t->vertices[2];
    // calculate vectors from point p to vertices p1, p2 and p3:
    glm::vec3 f1 = p1-p;
    glm::vec3 f2 = p2-p;
    glm::vec3 f3 = p3-p;
    // calculate the areas and factors (order of parameters doesn't matter):
    double a = glm::length(glm::cross(p1-p2, p1-p3)); // main triangle area a
    double alpha = glm::length(glm::cross(f2, f3)) / a; // p1's triangle area / a
    double beta = glm::length(glm::cross(f3, f1)) / a; // p2's triangle area / a
    double gamma = glm::length(glm::cross(f1, f2)) / a; // p3's triangle area / a
    std::vector<double> baryCoefs;
    baryCoefs.push_back(alpha);
    baryCoefs.push_back(beta);
    baryCoefs.push_back(gamma);
    return baryCoefs;
}

std::vector<double> barycentricCoefsNew(glm::vec3 P, ModelTriangle t){
    glm::vec3 v0 = t.vertices[0];
    glm::vec3 v1 = t.vertices[1];
    glm::vec3 v2 = t.vertices[2];
    glm::vec3 v0v1 = t.vertices[1] - t.vertices[0];
    glm::vec3 v0v2 = t.vertices[2] - t.vertices[0];
    // no need to normalize
    glm::vec3 N = glm::cross(v0v1,v0v2); // N
    float denom = glm::dot(N,N);


    glm::vec3 edge1 = v2 - v1;
    glm::vec3 vp1 = P - v1;

    glm::vec3 C = glm::cross(edge1,vp1);
    double u = glm::dot(N,C);

    // edge 2
    glm::vec3 edge2 = v0 - v2;
    glm::vec3 vp2 = P - v2;
    C = glm::cross(edge2,vp2);
    double v = glm::dot(N,C);

    u /= denom;
    v /= denom;

    std::vector<double> baryCoefs;
    baryCoefs.push_back(u);
    baryCoefs.push_back(v);
    baryCoefs.push_back(1-u-v);
    return baryCoefs;
}

//interpolates numberOfValues points between two vec3s
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues)
{
    std::vector<glm::vec3> ret;
    for (int i = 0; i < numberOfValues; i++)
    {
        glm::vec3 temp(0);
        for (int j = 0; j < 3; j++)
        {
            temp[j] = (from[j] + (i) * (to[j] - from[j]) / (numberOfValues - 1));
        }
        ret.push_back(temp);
    }
    return ret;
}

//does grayscale interpolation hting
void greyscale(DrawingWindow &window)
{
    std::vector<float> r = interpolateSingleFloats(255, 0, window.width);
    std::vector<float> g = interpolateSingleFloats(255, 0, window.width);
    std::vector<float> b = interpolateSingleFloats(255, 0, window.width);
    for (size_t y = 0; y < window.height; y++)
    {
        for (size_t x = 0; x < window.width; x++)
        {
            Colour col = Colour(r[x], g[x], b[x]);
            uint32_t packed = colToInt32(col);
            window.setPixelColour(x, y, packed);
        }
    }
}

//rgb interpolation
void rgbInterpolate(DrawingWindow &window)
{
    std::vector<glm::vec3> left;
    std::vector<glm::vec3> right;
    glm::vec3 topLeft(255, 0, 0);	   //red
    glm::vec3 topRight(0, 0, 255);	   //blu
    glm::vec3 bottomRight(0, 255, 0);  //grn
    glm::vec3 bottomLeft(255, 255, 0); //ylw
    left = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
    right = interpolateThreeElementValues(topRight, bottomRight, window.height);
    for (int i = 0; i < window.height; i++)
    {
        std::vector<glm::vec3> row = interpolateThreeElementValues(left[i], right[i], window.width);
        for (int j = 0; j < window.width; j++)
        {
            uint32_t packed = (255 << 24) + ((int)(row[j][0]) << 16) + ((int)(row[j][1]) << 8) + (int)(row[j][2]);
            window.setPixelColour(j, i, packed);
        }
    }
}